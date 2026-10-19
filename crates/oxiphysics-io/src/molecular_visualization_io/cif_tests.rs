// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the CIF reader.

use super::cif::{parse_cif, parse_cif_number};
use super::types::{CifBlock, MolVizError};

const ATOM_SITE_HEADER: &str = "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n";

fn sample_cif() -> &'static str {
    "data_test\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\nFe1 Fe 0.0 0.0 0.0 1.0\nFe2 Fe 0.5 0.5 0.5 1.0\n"
}

fn only_block(text: &str) -> CifBlock {
    let mut blocks = parse_cif(text).unwrap();
    assert_eq!(blocks.len(), 1, "expected exactly one data block");
    blocks.remove(0)
}

fn parse_message(text: &str) -> String {
    match parse_cif(text) {
        Err(MolVizError::ParseError { message }) => message,
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn test_cif_parse_block() {
    let blocks = parse_cif(sample_cif()).unwrap();
    assert_eq!(blocks.len(), 1);
    assert_eq!(blocks[0].name, "test");
}

#[test]
fn test_cif_parse_atom_sites() {
    let blocks = parse_cif(sample_cif()).unwrap();
    assert_eq!(blocks[0].atom_sites.len(), 2);
}

#[test]
fn test_cif_atom_site_values() {
    let blocks = parse_cif(sample_cif()).unwrap();
    let site0 = &blocks[0].atom_sites[0];
    assert_eq!(site0.element, "Fe");
    assert!((site0.fx - 0.0).abs() < 1e-12);
    assert!((site0.occupancy - 1.0).abs() < 1e-12);
}

#[test]
fn test_cif_parse_tags() {
    let blocks = parse_cif(sample_cif()).unwrap();
    assert!(blocks[0].tags.contains_key("_cell_length_a"));
}

#[test]
fn standard_uncertainties_are_dropped_from_coordinates_and_occupancy() {
    let text = format!("data_x\n{ATOM_SITE_HEADER}Fe1 Fe 0.2500(3) -0.1250(12) 0.5 0.95(2)\n");
    let blk = only_block(&text);
    assert_eq!(blk.atom_sites.len(), 1);
    let site = &blk.atom_sites[0];
    assert_eq!(site.label, "Fe1");
    assert_eq!(site.fx, 0.25);
    assert_eq!(site.fy, -0.125);
    assert_eq!(site.fz, 0.5);
    assert_eq!(site.occupancy, 0.95);
}

#[test]
fn cell_tags_after_the_atom_loop_are_kept() {
    let text = format!(
        "data_x\n{ATOM_SITE_HEADER}Fe1 Fe 0.0 0.0 0.0 1.0\nFe2 Fe 0.5 0.5 0.5 1.0\n_cell_length_a 4.0\n_cell_length_b 5.0(1)\n_cell_angle_gamma 90\n"
    );
    let blk = only_block(&text);
    assert_eq!(blk.atom_sites.len(), 2, "tags must not become table rows");
    assert_eq!(blk.tags["_cell_length_a"], "4.0");
    assert_eq!(blk.tags["_cell_length_b"], "5.0(1)");
    assert_eq!(blk.tags["_cell_angle_gamma"], "90");
    assert_eq!(
        parse_cif_number(&blk.tags["_cell_length_b"]).unwrap(),
        Some(5.0)
    );
}

#[test]
fn a_second_table_after_the_atom_loop_is_not_mixed_into_the_sites() {
    let text = format!(
        "data_x\n{ATOM_SITE_HEADER}Fe1 Fe 0.1 0.2 0.3 1.0\nloop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n'-x, -y, -z'\n_cell_volume 64.0\n"
    );
    let blk = only_block(&text);
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.tags["_cell_volume"], "64.0");
    assert!(!blk.tags.contains_key("_symmetry_equiv_pos_as_xyz"));
}

#[test]
fn anisotropic_displacement_tables_do_not_create_atom_sites() {
    let text = format!(
        "data_x\n{ATOM_SITE_HEADER}Fe1 Fe 0.1 0.2 0.3 1.0\nloop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\n_atom_site_aniso_U_22\n_atom_site_aniso_U_33\nFe1 0.01 0.02 0.03\n"
    );
    let blk = only_block(&text);
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].fx, 0.1);
}

#[test]
fn quoted_values_keep_their_blanks_and_quotes_end_before_whitespace() {
    let text = format!(
        "data_x\n_symmetry_space_group_name_H-M 'P 21/c'\n_chemical_name \"O'Brien salt\"\n_other 'it's here'\n{ATOM_SITE_HEADER}'Fe 1' 'Fe 3+' 0.0 0.0 0.0 1.0\n"
    );
    let blk = only_block(&text);
    assert_eq!(blk.tags["_symmetry_space_group_name_H-M"], "P 21/c");
    assert_eq!(blk.tags["_chemical_name"], "O'Brien salt");
    assert_eq!(blk.tags["_other"], "it's here");
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].label, "Fe 1");
    assert_eq!(blk.atom_sites[0].element, "Fe 3+");
}

#[test]
fn table_rows_may_span_lines_and_share_a_line() {
    let text = "data_x\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nA1 0.1\n0.2 0.3 A2 0.4 0.5\n0.6\nA3 0.7 0.8 0.9\n";
    let blk = only_block(text);
    let sites: Vec<(&str, f64, f64, f64)> = blk
        .atom_sites
        .iter()
        .map(|s| (s.label.as_str(), s.fx, s.fy, s.fz))
        .collect();
    assert_eq!(
        sites,
        vec![
            ("A1", 0.1, 0.2, 0.3),
            ("A2", 0.4, 0.5, 0.6),
            ("A3", 0.7, 0.8, 0.9)
        ]
    );
}

#[test]
fn values_on_the_line_after_their_data_name_are_found() {
    let blk = only_block("data_x\n_cell_length_a\n  4.25(2)\n_cell_length_b   \t 5.5\n");
    assert_eq!(blk.tags["_cell_length_a"], "4.25(2)");
    assert_eq!(blk.tags["_cell_length_b"], "5.5");
}

#[test]
fn text_fields_are_values() {
    let text = "data_x\n_publ_section_title\n;\n Crystal structure\n of iron\n;\n_cell_length_a 2.87\n_note\n;trailing text\n;\n";
    let blk = only_block(text);
    assert_eq!(
        blk.tags["_publ_section_title"],
        " Crystal structure\n of iron"
    );
    assert_eq!(blk.tags["_cell_length_a"], "2.87");
    assert_eq!(blk.tags["_note"], "trailing text");
}

#[test]
fn text_fields_can_hold_table_cells() {
    let text = "data_x\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n;site\nwith newline\n;\n0.1 0.2 0.3\n";
    let blk = only_block(text);
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].label, "site\nwith newline");
    assert_eq!(blk.atom_sites[0].fx, 0.1);
}

#[test]
fn comments_are_ignored_everywhere() {
    let text = "# leading\ndata_x # block\n_cell_length_a 4.0 # trailing\nloop_ # table\n_atom_site_label\n_atom_site_fract_x # x\n_atom_site_fract_y\n_atom_site_fract_z\nFe 0.1 # first\n 0.2 0.3\n# done\n";
    let blk = only_block(text);
    assert_eq!(blk.name, "x");
    assert_eq!(blk.tags["_cell_length_a"], "4.0");
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].fz, 0.3);
}

#[test]
fn hash_inside_a_word_or_quotes_is_not_a_comment() {
    let blk = only_block("data_x\n_a abc#def\n_b 'x # y'\n");
    assert_eq!(blk.tags["_a"], "abc#def");
    assert_eq!(blk.tags["_b"], "x # y");
}

#[test]
fn keywords_and_data_names_are_case_insensitive() {
    let text = "DATA_Upper\n_Cell_Length_A 3.0\nLOOP_\n_ATOM_SITE_LABEL\n_Atom_Site_Fract_X\n_atom_site_fract_Y\n_atom_site_fract_z\nC1 0.1 0.2 0.3\n";
    let blk = only_block(text);
    assert_eq!(blk.name, "Upper");
    assert_eq!(blk.tags["_Cell_Length_A"], "3.0");
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].label, "C1");
}

#[test]
fn dotted_category_separators_are_accepted() {
    let text = "data_x\nloop_\n_atom_site.label\n_atom_site.type_symbol\n_atom_site.fract_x\n_atom_site.fract_y\n_atom_site.fract_z\nSi1 Si 0.125 0.25 0.375\n";
    let blk = only_block(text);
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].element, "Si");
    assert_eq!(blk.atom_sites[0].fz, 0.375);
}

#[test]
fn windows_line_endings_and_a_byte_order_mark_are_handled() {
    let text = format!("\u{feff}{}", sample_cif().replace('\n', "\r\n"));
    let blk = only_block(&text);
    assert_eq!(blk.name, "test");
    assert_eq!(blk.tags["_cell_length_c"], "4.0");
    assert_eq!(blk.atom_sites.len(), 2);
    assert_eq!(blk.atom_sites[1].label, "Fe2");
    assert_eq!(blk.atom_sites[1].fz, 0.5);

    let field = "data_x\r\n_t\r\n;\r\nline one\r\nline two\r\n;\r\n_u 1\r\n";
    let blk = only_block(field);
    assert_eq!(blk.tags["_t"], "line one\r\nline two");
    assert_eq!(blk.tags["_u"], "1");
}

#[test]
fn unknown_and_inapplicable_occupancy_defaults_to_one() {
    let text = format!(
        "data_x\n{ATOM_SITE_HEADER}A 1 0.1 0.2 0.3 ?\nB 1 0.1 0.2 0.3 .\nC 1 0.1 0.2 0.3 0.5(1)\n"
    );
    let blk = only_block(&text);
    let occupancies: Vec<f64> = blk.atom_sites.iter().map(|s| s.occupancy).collect();
    assert_eq!(occupancies, vec![1.0, 1.0, 0.5]);
}

#[test]
fn missing_label_type_and_occupancy_columns_use_defaults() {
    let text =
        "data_x\nloop_\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n0.1 0.2 0.3\n";
    let blk = only_block(text);
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].label, "X");
    assert_eq!(blk.atom_sites[0].element, "X");
    assert_eq!(blk.atom_sites[0].occupancy, 1.0);
}

#[test]
fn the_element_is_taken_from_the_label_without_a_type_column() {
    let text = "data_x\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nFe1 0 0 0\nO2 0 0 0\nCa 0 0 0\nH12A 0 0 0\n";
    let blk = only_block(text);
    let elements: Vec<&str> = blk.atom_sites.iter().map(|s| s.element.as_str()).collect();
    assert_eq!(elements, vec!["Fe", "O", "Ca", "H"]);
}

#[test]
fn unknown_coordinates_are_an_error_naming_the_site() {
    let text = format!("data_x\n{ATOM_SITE_HEADER}Fe1 Fe ? 0.0 0.0 1.0\n");
    let message = parse_message(&text);
    assert!(message.contains("Fe1"), "{message}");
    assert!(message.contains("_atom_site_fract_x"), "{message}");
}

#[test]
fn garbage_coordinates_and_occupancies_are_errors_not_zero() {
    for row in [
        "Fe1 Fe abc 0.0 0.0 1.0",
        "Fe1 Fe 0.0 0.0 0.5(x) 1.0",
        "Fe1 Fe 0.0 nan 0.0 1.0",
        "Fe1 Fe 0.0 0.0 0.0 high",
    ] {
        let text = format!("data_x\n{ATOM_SITE_HEADER}{row}\n");
        let message = parse_message(&text);
        assert!(message.contains("Fe1"), "{row}: {message}");
    }
}

#[test]
fn an_incomplete_last_row_is_an_error() {
    let text = format!("data_x\n{ATOM_SITE_HEADER}Fe1 Fe 0.0 0.0 0.0 1.0\nFe2 Fe 0.5 0.5\n");
    let message = parse_message(&text);
    assert!(message.contains("incomplete"), "{message}");
}

#[test]
fn unterminated_quotes_and_text_fields_are_errors() {
    let message = parse_message("data_x\n_a 'never closed\n_b 1\n");
    assert!(message.contains("line 2"), "{message}");
    let message = parse_message("data_x\n_a\n;\nnever closed\n_b 1\n");
    assert!(message.contains("line 3"), "{message}");
}

#[test]
fn malformed_structure_is_reported_with_the_line() {
    let message = parse_message("data_x\n_a\n_b 1\n");
    assert!(
        message.contains("_a") && message.contains("line 2"),
        "{message}"
    );
    let message = parse_message("data_x\n_a 1 2\n");
    assert!(
        message.contains("'2'") && message.contains("line 2"),
        "{message}"
    );
    let message = parse_message("data_x\nloop_\nvalue\n");
    assert!(message.contains("no data names"), "{message}");
    let message = parse_message("data_x\nloop_\nloop_\n_a\n1\n");
    assert!(message.contains("loop_"), "{message}");
    let message = parse_message("data_x\n_a\n");
    assert!(message.contains("_a"), "{message}");
}

#[test]
fn several_blocks_keep_their_own_tags_and_sites() {
    let text = format!(
        "ignored preamble _x 1\ndata_first\n_a 1\n{ATOM_SITE_HEADER}A 1 0.1 0.2 0.3 1\ndata_second\n_a 2\n{ATOM_SITE_HEADER}B 1 0.4 0.5 0.6 1\nC 1 0.7 0.8 0.9 1\ndata_third\n"
    );
    let blocks = parse_cif(&text).unwrap();
    assert_eq!(blocks.len(), 3);
    assert_eq!(blocks[0].name, "first");
    assert_eq!(blocks[0].tags["_a"], "1");
    assert_eq!(blocks[0].atom_sites.len(), 1);
    assert_eq!(blocks[1].name, "second");
    assert_eq!(blocks[1].tags["_a"], "2");
    assert_eq!(blocks[1].atom_sites.len(), 2);
    assert_eq!(blocks[1].atom_sites[1].fx, 0.7);
    assert_eq!(blocks[2].name, "third");
    assert!(blocks[2].tags.is_empty() && blocks[2].atom_sites.is_empty());
}

#[test]
fn empty_input_and_input_without_blocks_give_no_blocks() {
    assert!(parse_cif("").unwrap().is_empty());
    assert!(parse_cif("# only a comment\n\n").unwrap().is_empty());
    assert!(parse_cif("_cell_length_a 4.0\n").unwrap().is_empty());
}

#[test]
fn a_realistic_database_entry_round_trips_its_numbers() {
    let text = "\
#------------------------------------------------------------------------------
data_9008460
_chemical_name_mineral 'Iron'
_cell_length_a 2.8665(3)
_cell_length_b 2.8665(3)
_cell_length_c 2.8665(3)
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'I m -3 m'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-x,-y,z
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
Fe 0.0000(0) 0.00000 0.000(0) 0.0123(4)
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
Fe 0.0123(4)
";
    let blk = only_block(text);
    assert_eq!(blk.name, "9008460");
    assert_eq!(blk.tags["_symmetry_space_group_name_H-M"], "I m -3 m");
    assert_eq!(blk.atom_sites.len(), 1);
    assert_eq!(blk.atom_sites[0].label, "Fe");
    assert_eq!(blk.atom_sites[0].element, "Fe");
    assert_eq!(
        (
            blk.atom_sites[0].fx,
            blk.atom_sites[0].fy,
            blk.atom_sites[0].fz
        ),
        (0.0, 0.0, 0.0)
    );
    assert_eq!(
        parse_cif_number(&blk.tags["_cell_length_a"]).unwrap(),
        Some(2.8665)
    );
}

#[test]
fn cif_numbers_drop_their_standard_uncertainty() {
    let cases = [
        ("0.2500(3)", 0.25),
        ("-0.0012(4)", -0.0012),
        ("12(3)", 12.0),
        ("1.5e-3(2)", 0.0015),
        ("  7.5 ", 7.5),
        ("+.5", 0.5),
        ("5.", 5.0),
        ("-3", -3.0),
        ("1.0E2", 100.0),
        ("0.25(12345)", 0.25),
    ];
    for (text, want) in cases {
        assert_eq!(parse_cif_number(text).unwrap(), Some(want), "{text}");
    }
}

#[test]
fn cif_number_markers_for_unknown_and_inapplicable_are_none() {
    assert_eq!(parse_cif_number("?").unwrap(), None);
    assert_eq!(parse_cif_number(".").unwrap(), None);
    assert_eq!(parse_cif_number(" ? ").unwrap(), None);
}

#[test]
fn malformed_cif_numbers_are_errors() {
    for text in [
        "", "abc", "(3)", "0.25(", "0.25()", "0.25(a)", "0.25(1", "nan", "inf", "-inf", "1e999",
        "0.25(3)x", "1,5", "0x10",
    ] {
        assert!(
            matches!(parse_cif_number(text), Err(MolVizError::ParseError { .. })),
            "{text:?} must be rejected"
        );
    }
}
