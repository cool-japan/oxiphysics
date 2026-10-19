// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Interpretation of the atom-site table and the symmetry-operation list of a
//! parsed CIF data block.

use super::cif_syntax::{CifBlock, CifLoop};
use super::lattice_math;
use super::symmetry::{SymmetryOp, parse_symmetry_op};
use super::{CifParser, CrystalAtom};
use crate::{Error, Result};
use std::f64::consts::PI;

/// Data names that can hold the fractional x coordinate of an atom site.
pub(super) const FRACT_X_TAGS: [&str; 2] = ["_atom_site_fract_x", "_atom_site_x_fract"];
/// Data name of the Cartesian x coordinate of an atom site.
pub(super) const CARTN_X_TAGS: [&str; 1] = ["_atom_site_cartn_x"];

/// Chemical element symbols, used to recognise the element in atom labels.
const ELEMENT_SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Reduce an atom type symbol or site label to a chemical element symbol.
///
/// Oxidation states and site numbering are dropped: `Fe3+`, `O2-`, `Fe1` and
/// `O2A` become `Fe`, `O`, `Fe` and `O`.  The leading letters are read as a
/// two-letter symbol when that is a real element (`Ca`, `FE` -> `Fe`) and as a
/// one-letter symbol otherwise (`OW` -> `O`).  A name that matches no element,
/// such as `D`, is returned capitalised so that it is not lost.
pub(super) fn element_symbol(raw: &str) -> String {
    let letters: Vec<char> = raw
        .chars()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    let Some(&first) = letters.first() else {
        return raw.trim().to_string();
    };
    let one = first.to_ascii_uppercase().to_string();
    if let Some(&second) = letters.get(1) {
        let two = format!("{one}{}", second.to_ascii_lowercase());
        if ELEMENT_SYMBOLS.contains(&two.as_str()) {
            return two;
        }
    }
    if ELEMENT_SYMBOLS.contains(&one.as_str()) {
        return one;
    }
    let mut capitalised = one;
    capitalised.extend(letters.iter().skip(1).map(|c| c.to_ascii_lowercase()));
    capitalised
}

/// Whether a CIF value is the "unknown" (`?`) or "not applicable" (`.`)
/// placeholder.
fn is_placeholder(token: &str) -> bool {
    matches!(token.trim(), "?" | ".")
}

/// Numeric value of a CIF token, ignoring a standard uncertainty `1.234(5)`.
fn cif_number(token: &str) -> Option<f64> {
    CifParser::strip_su(token.trim())
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
}

fn cell(row: &[String], column: Option<usize>) -> Option<&str> {
    column.and_then(|c| row.get(c)).map(String::as_str)
}

/// Column indices of a coordinate triple, if the table has all three columns.
fn coordinate_columns(table: &CifLoop, tags: [&[&str]; 3]) -> Option<[usize; 3]> {
    Some([
        table.column(tags[0])?,
        table.column(tags[1])?,
        table.column(tags[2])?,
    ])
}

/// Read the atoms of an `_atom_site` table.
///
/// Positions are taken from the `_atom_site_fract_*` columns, or converted
/// from `_atom_site_Cartn_*` through `lattice` when only those exist.  The
/// displacement parameter is the isotropic B factor in square Angstrom; when
/// only `U_iso` is given it is converted with `B = 8 pi^2 U`.  A missing or
/// `?` occupancy means 1.
///
/// # Errors
///
/// Fails when the table has neither atom type nor label column, lacks a full
/// coordinate triple, contains a coordinate or occupancy that is not a
/// number, or needs a Cartesian conversion with a singular lattice.
pub(super) fn extract_atoms(table: &CifLoop, lattice: &[[f64; 3]; 3]) -> Result<Vec<CrystalAtom>> {
    let i_symbol = table.column(&["_atom_site_type_symbol"]);
    let i_label = table.column(&["_atom_site_label"]);
    if i_symbol.is_none() && i_label.is_none() {
        return Err(Error::Parse(
            "CIF: atom-site table has neither _atom_site_type_symbol nor _atom_site_label"
                .to_string(),
        ));
    }
    let fractional = coordinate_columns(
        table,
        [
            &FRACT_X_TAGS,
            &["_atom_site_fract_y", "_atom_site_y_fract"],
            &["_atom_site_fract_z", "_atom_site_z_fract"],
        ],
    );
    let cartesian = coordinate_columns(
        table,
        [
            &CARTN_X_TAGS,
            &["_atom_site_cartn_y"],
            &["_atom_site_cartn_z"],
        ],
    );
    let (columns, is_cartesian) = match (fractional, cartesian) {
        (Some(columns), _) => (columns, false),
        (None, Some(columns)) => (columns, true),
        (None, None) => {
            return Err(Error::Parse(
                "CIF: atom-site table lacks a complete set of fractional or Cartesian coordinates"
                    .to_string(),
            ));
        }
    };
    let lattice_inverse = if is_cartesian {
        Some(lattice_math::inverse(lattice).ok_or_else(|| {
            Error::Parse(
                "CIF: unit cell is singular, cannot convert Cartesian atom coordinates".to_string(),
            )
        })?)
    } else {
        None
    };
    let i_occupancy = table.column(&["_atom_site_occupancy"]);
    let i_b = table.column(&["_atom_site_b_iso_or_equiv", "_atom_site_b_iso"]);
    let i_u = table.column(&["_atom_site_u_iso_or_equiv", "_atom_site_u_iso"]);

    let mut atoms = Vec::new();
    for row in table.rows()? {
        let name = cell(row, i_symbol)
            .filter(|token| !is_placeholder(token))
            .or_else(|| cell(row, i_label))
            .unwrap_or_default();
        let label = cell(row, i_label).unwrap_or(name);
        let element = element_symbol(name);
        if element.is_empty() {
            return Err(Error::Parse(format!(
                "CIF: atom '{label}' has an empty element symbol"
            )));
        }
        let mut position = [0.0_f64; 3];
        for (axis, (slot, column)) in position.iter_mut().zip(columns).enumerate() {
            let token = cell(row, Some(column)).unwrap_or_default();
            *slot = cif_number(token).ok_or_else(|| {
                Error::Parse(format!(
                    "CIF: atom '{label}' has an invalid coordinate '{token}' on axis {axis}"
                ))
            })?;
        }
        if let Some(inverse) = &lattice_inverse {
            position = lattice_math::cartesian_to_fractional(position, inverse);
        }
        let occupancy = match cell(row, i_occupancy) {
            None => 1.0,
            Some(token) if is_placeholder(token) => 1.0,
            Some(token) => cif_number(token).ok_or_else(|| {
                Error::Parse(format!(
                    "CIF: atom '{label}' has an invalid occupancy '{token}'"
                ))
            })?,
        };
        let b_iso = cell(row, i_b).and_then(cif_number);
        let u_iso = cell(row, i_u).and_then(cif_number);
        let b_factor = match (b_iso, u_iso) {
            (Some(b), _) => b,
            (None, Some(u)) => 8.0 * PI * PI * u,
            (None, None) => 0.0,
        };
        atoms.push(CrystalAtom {
            element,
            fractional_pos: position,
            occupancy,
            b_factor,
        });
    }
    Ok(atoms)
}

/// Collect the space-group operations listed in a data block.
///
/// The newer `_space_group_symop_operation_xyz` is preferred over the legacy
/// `_symmetry_equiv_pos_as_xyz`; either may appear as a table or, for a single
/// operation, as a plain item.  An empty result means the file lists none.
///
/// # Errors
///
/// Fails if a listed operation cannot be parsed.
pub(super) fn collect_symmetry_ops(block: &CifBlock) -> Result<Vec<SymmetryOp>> {
    const TAGS: [&str; 2] = [
        "_space_group_symop_operation_xyz",
        "_symmetry_equiv_pos_as_xyz",
    ];
    for tag in TAGS {
        let mut ops = Vec::new();
        for table in &block.loops {
            let Some(column) = table.column(&[tag]) else {
                continue;
            };
            for row in table.rows()? {
                let text = row.get(column).map(String::as_str).unwrap_or_default();
                ops.push(parse_symmetry_op(text)?);
            }
        }
        if let Some(text) = block.items.get(tag)
            && !text.is_empty()
        {
            ops.push(parse_symmetry_op(text)?);
        }
        if !ops.is_empty() {
            return Ok(ops);
        }
    }
    Ok(Vec::new())
}

#[cfg(test)]
mod tests {
    use super::super::cif_syntax::parse_block;
    use super::*;

    const CUBIC_5: [[f64; 3]; 3] = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]];

    fn table_of(cif: &str) -> CifLoop {
        parse_block(cif)
            .expect("valid CIF")
            .loops
            .into_iter()
            .next()
            .expect("one loop")
    }

    const SITE_HEADER: &str = "\
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
";

    #[test]
    fn element_symbols_drop_charges_and_site_numbers() {
        for (raw, expected) in [
            ("Fe", "Fe"),
            ("Fe3+", "Fe"),
            ("O2-", "O"),
            ("Fe1", "Fe"),
            ("O2A", "O"),
            ("Si_3", "Si"),
            ("Ca1", "Ca"),
            ("FE", "Fe"),
            ("cu", "Cu"),
            ("OW", "O"),
            ("H", "H"),
            ("Og", "Og"),
            ("D", "D"),
            ("Xx", "Xx"),
        ] {
            assert_eq!(element_symbol(raw), expected, "{raw}");
        }
        assert_eq!(element_symbol("1"), "1");
        assert_eq!(element_symbol(""), "");
    }

    #[test]
    fn every_element_symbol_is_a_fixed_point() {
        for symbol in ELEMENT_SYMBOLS {
            assert_eq!(element_symbol(symbol), symbol);
            assert_eq!(element_symbol(&format!("{symbol}12")), symbol);
        }
    }

    #[test]
    fn uiso_is_converted_to_a_b_factor() {
        let table = table_of(&format!(
            "{SITE_HEADER}_atom_site_U_iso_or_equiv\nFe1 0.0 0.0 0.0 0.01\n"
        ));
        let atoms = extract_atoms(&table, &CUBIC_5).unwrap();
        assert_eq!(atoms.len(), 1);
        let expected = 8.0 * PI * PI * 0.01;
        assert!((atoms[0].b_factor - expected).abs() < 1e-12);
        assert!((atoms[0].b_factor - 0.789_568_352).abs() < 1e-8);
        assert_eq!(atoms[0].element, "Fe");
    }

    #[test]
    fn b_iso_wins_over_u_iso_and_placeholders_fall_back() {
        let table = table_of(
            "\
loop_
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_B_iso_or_equiv
_atom_site_U_iso_or_equiv
Na 0.1(1) 0.2 0.3 0.5 1.5 0.9
Cl 0 0 0 . ? 0.02
Ar 0 0 0 ? . ?
",
        );
        let atoms = extract_atoms(&table, &CUBIC_5).unwrap();
        assert_eq!(atoms.len(), 3);
        assert_eq!(atoms[0].element, "Na");
        assert_eq!(atoms[0].fractional_pos, [0.1, 0.2, 0.3]);
        assert_eq!(atoms[0].occupancy, 0.5);
        assert_eq!(atoms[0].b_factor, 1.5);
        assert_eq!(atoms[1].occupancy, 1.0);
        assert!((atoms[1].b_factor - 8.0 * PI * PI * 0.02).abs() < 1e-12);
        assert_eq!(atoms[2].occupancy, 1.0);
        assert_eq!(atoms[2].b_factor, 0.0);
    }

    #[test]
    fn the_label_names_the_element_when_there_is_no_type_symbol() {
        let table = table_of(&format!("{SITE_HEADER}O2A 0 0 0\nFe3 0.5 0.5 0.5\n"));
        let atoms = extract_atoms(&table, &CUBIC_5).unwrap();
        assert_eq!(atoms[0].element, "O");
        assert_eq!(atoms[1].element, "Fe");
    }

    #[test]
    fn an_anisotropic_displacement_table_is_not_an_atom_table() {
        let aniso = table_of(
            "\
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
Fe1 0.01 0.02
",
        );
        let err = extract_atoms(&aniso, &CUBIC_5).unwrap_err();
        assert!(err.to_string().contains("neither"), "{err}");
    }

    #[test]
    fn tables_without_a_full_coordinate_triple_are_rejected() {
        let table = table_of(
            "\
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
Fe1 0.0 0.0
",
        );
        let err = extract_atoms(&table, &CUBIC_5).unwrap_err();
        assert!(err.to_string().contains("coordinates"), "{err}");
    }

    #[test]
    fn unreadable_numbers_are_errors_not_zeros() {
        for bad in ["?", ".", "abc", "1.2.3"] {
            let table = table_of(&format!("{SITE_HEADER}Fe1 {bad} 0.0 0.0\n"));
            let err = extract_atoms(&table, &CUBIC_5).unwrap_err();
            assert!(err.to_string().contains("Fe1"), "{bad}: {err}");
            assert!(err.to_string().contains("coordinate"), "{bad}: {err}");
        }
        let table = table_of(&format!(
            "{SITE_HEADER}_atom_site_occupancy\nFe1 0 0 0 half\n"
        ));
        let err = extract_atoms(&table, &CUBIC_5).unwrap_err();
        assert!(err.to_string().contains("occupancy"), "{err}");
    }

    #[test]
    fn a_ragged_atom_table_is_an_error() {
        let table = table_of(&format!("{SITE_HEADER}Fe1 0 0 0\nFe2 0.5 0.5\n"));
        assert!(extract_atoms(&table, &CUBIC_5).is_err());
    }

    #[test]
    fn cartesian_sites_are_converted_with_the_cell() {
        let hexagonal = CifParser::build_lattice(
            3.0,
            3.0,
            5.0,
            90.0_f64.to_radians(),
            90.0_f64.to_radians(),
            120.0_f64.to_radians(),
        );
        let table = table_of(
            "\
loop_
_atom_site_label
_atom_site_Cartn_x
_atom_site_Cartn_y
_atom_site_Cartn_z
Mg1 0.0 1.7320508075688772 2.5
",
        );
        let atoms = extract_atoms(&table, &hexagonal).unwrap();
        let pos = atoms[0].fractional_pos;
        assert!((pos[0] - 1.0 / 3.0).abs() < 1e-12, "{pos:?}");
        assert!((pos[1] - 2.0 / 3.0).abs() < 1e-12, "{pos:?}");
        assert!((pos[2] - 0.5).abs() < 1e-12, "{pos:?}");
        let singular = [[0.0; 3]; 3];
        assert!(extract_atoms(&table, &singular).is_err());
    }

    fn block_of(cif: &str) -> CifBlock {
        parse_block(cif).expect("valid CIF")
    }

    #[test]
    fn symmetry_operations_are_read_from_a_two_column_loop() {
        let block = block_of(
            "\
data_x
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x, y, z'
2 '-x, -y, -z'
3 'x+1/2, y+1/2, z'
",
        );
        let ops = collect_symmetry_ops(&block).unwrap();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[2].translation, [0.5, 0.5, 0.0]);
    }

    #[test]
    fn the_newer_tag_is_preferred_and_single_items_work() {
        let both = block_of(
            "\
data_x
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-x,-y,-z
loop_
_space_group_symop_operation_xyz
x,y,z
",
        );
        assert_eq!(collect_symmetry_ops(&both).unwrap().len(), 1);
        let single = block_of("data_x\n_symmetry_equiv_pos_as_xyz 'x, y, z'\n");
        assert_eq!(collect_symmetry_ops(&single).unwrap().len(), 1);
        assert!(
            collect_symmetry_ops(&block_of("data_x\n_cell_length_a 4\n"))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_malformed_operation_fails_the_whole_read() {
        let block = block_of(
            "\
data_x
loop_
_space_group_symop_operation_xyz
x,y,z
x,y
",
        );
        assert!(collect_symmetry_ops(&block).is_err());
    }
}
