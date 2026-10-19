// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the CIF parser: table handling, displacement
//! parameters and space-group symmetry expansion.

use super::*;
use std::f64::consts::PI;

fn parse(cif: &str) -> CrystalStructure {
    CifParser::new().parse(cif).expect("valid CIF")
}

fn assert_close(got: [f64; 3], want: [f64; 3]) {
    for (g, w) in got.iter().zip(want.iter()) {
        assert!((g - w).abs() < 1e-9, "{got:?} vs {want:?}");
    }
}

fn cubic_cell(a: f64) -> String {
    format!(
        "_cell_length_a {a}\n_cell_length_b {a}\n_cell_length_c {a}\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
    )
}

const SITE_COLUMNS: &str = "\
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
";

#[test]
fn the_atom_site_table_survives_a_trailing_aniso_table() {
    let cif = format!(
        "data_test\n{cell}{SITE_COLUMNS}_atom_site_U_iso_or_equiv\n\
         Fe1 Fe 0.0 0.0 0.0 0.010\n\
         O1 O 0.25 0.5 0.75 0.020\n\
         loop_\n_atom_site_aniso_label\n_atom_site_aniso_type_symbol\n\
         _atom_site_aniso_U_11\n_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n\
         Fe1 Fe 0.01 0.01 0.01\n\
         O1 O 0.02 0.02 0.02\n",
        cell = cubic_cell(5.0)
    );
    let crystal = parse(&cif);
    assert_eq!(crystal.atoms.len(), 2);
    assert_eq!(crystal.atoms[0].element, "Fe");
    assert_eq!(crystal.atoms[1].element, "O");
    assert_close(crystal.atoms[0].fractional_pos, [0.0, 0.0, 0.0]);
    assert_close(crystal.atoms[1].fractional_pos, [0.25, 0.5, 0.75]);
    assert!((crystal.volume() - 125.0).abs() < 1e-9);
}

#[test]
fn the_atom_site_table_may_come_before_other_tables_and_items() {
    // Aniso table first, atom sites second, publication table last.
    let cif = format!(
        "data_test\n{cell}\
         loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\nFe1 0.01\n\
         {SITE_COLUMNS}Fe1 Fe 0.5 0.5 0.5\n\
         loop_\n_publ_author_name\n'Doe, J.'\n'Roe, R.'\n",
        cell = cubic_cell(4.0)
    );
    let crystal = parse(&cif);
    assert_eq!(crystal.atoms.len(), 1);
    assert_close(crystal.atoms[0].fractional_pos, [0.5, 0.5, 0.5]);
}

#[test]
fn items_after_a_loop_are_read_as_items_again() {
    let cif = format!(
        "data_test\n{SITE_COLUMNS}Si1 Si 0 0 0\n\
         _cell_length_a 5.43\n_cell_length_b 5.43\n_cell_length_c 5.43\n\
         _symmetry_Int_Tables_number 227\n"
    );
    let crystal = parse(&cif);
    assert!((crystal.lattice[0][0] - 5.43).abs() < 1e-12);
    assert!((crystal.volume() - 5.43_f64.powi(3)).abs() < 1e-9);
    assert_eq!(crystal.space_group, 227);
    assert_eq!(crystal.atoms.len(), 1);
}

#[test]
fn u_iso_becomes_b_iso_with_8_pi_squared() {
    let u_only = format!(
        "data_u\n{cell}{SITE_COLUMNS}_atom_site_U_iso_or_equiv\nFe1 Fe 0 0 0 0.0125(4)\n",
        cell = cubic_cell(3.0)
    );
    let atom = &parse(&u_only).atoms[0];
    assert!((atom.b_factor - 8.0 * PI * PI * 0.0125).abs() < 1e-12);

    let b_only = format!(
        "data_b\n{cell}{SITE_COLUMNS}_atom_site_B_iso_or_equiv\nFe1 Fe 0 0 0 0.987\n",
        cell = cubic_cell(3.0)
    );
    assert!((parse(&b_only).atoms[0].b_factor - 0.987).abs() < 1e-12);

    let neither = format!(
        "data_n\n{cell}{SITE_COLUMNS}Fe1 Fe 0 0 0\n",
        cell = cubic_cell(3.0)
    );
    assert_eq!(parse(&neither).atoms[0].b_factor, 0.0);
}

#[test]
fn centering_operations_expand_a_single_atom_to_the_conventional_cell() {
    let cif = format!(
        "data_Cu\n{cell}_symmetry_Int_Tables_number 225\n\
         loop_\n_symmetry_equiv_pos_as_xyz\n\
         'x, y, z'\n'x, y+1/2, z+1/2'\n'x+1/2, y, z+1/2'\n'x+1/2, y+1/2, z'\n\
         {SITE_COLUMNS}Cu1 Cu 0.0 0.0 0.0\n",
        cell = cubic_cell(3.615)
    );
    let crystal = parse(&cif);
    assert_eq!(crystal.atoms.len(), 4);
    let expected = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    for (atom, want) in crystal.atoms.iter().zip(expected) {
        assert_eq!(atom.element, "Cu");
        assert_close(atom.fractional_pos, want);
    }
    assert_eq!(crystal.space_group, 225);
}

#[test]
fn rock_salt_is_completed_from_its_asymmetric_unit() {
    let cif = format!(
        "data_NaCl\n{cell}\
         loop_\n_space_group_symop_id\n_space_group_symop_operation_xyz\n\
         1 x,y,z\n2 x,y+1/2,z+1/2\n3 x+1/2,y,z+1/2\n4 x+1/2,y+1/2,z\n\
         {SITE_COLUMNS}Na1 Na 0 0 0\nCl1 Cl 0.5 0.5 0.5\n",
        cell = cubic_cell(5.64)
    );
    let crystal = parse(&cif);
    let count = |element: &str| {
        crystal
            .atoms
            .iter()
            .filter(|a| a.element == element)
            .count()
    };
    assert_eq!(count("Na"), 4);
    assert_eq!(count("Cl"), 4);
    for atom in &crystal.atoms {
        for c in atom.fractional_pos {
            assert!((0.0..1.0).contains(&c), "{:?}", atom.fractional_pos);
        }
    }
    // Every Na has six Cl neighbours at a / 2 once periodic images are counted.
    let cart = crystal.to_cartesian();
    let a = 5.64;
    for (i, _) in crystal
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, x)| x.element == "Na")
    {
        let mut neighbours = 0;
        for (j, _) in crystal
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, x)| x.element == "Cl")
        {
            for sx in -1..=1 {
                for sy in -1..=1 {
                    for sz in -1..=1 {
                        let dx = cart[j][0] + f64::from(sx) * a - cart[i][0];
                        let dy = cart[j][1] + f64::from(sy) * a - cart[i][1];
                        let dz = cart[j][2] + f64::from(sz) * a - cart[i][2];
                        let distance = (dx * dx + dy * dy + dz * dz).sqrt();
                        if (distance - a / 2.0).abs() < 1e-9 {
                            neighbours += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(neighbours, 6, "Na {i}");
    }
}

#[test]
fn atoms_on_special_positions_are_not_duplicated() {
    let ops = "loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n-x,-y,-z\n";
    let at_inversion_centre = format!(
        "data_p1bar\n{cell}{ops}{SITE_COLUMNS}Pt1 Pt 0 0 0\n",
        cell = cubic_cell(4.0)
    );
    assert_eq!(parse(&at_inversion_centre).atoms.len(), 1);
    let general = format!(
        "data_p1bar\n{cell}{ops}{SITE_COLUMNS}Pt1 Pt 0.1 0.2 0.3\n",
        cell = cubic_cell(4.0)
    );
    let crystal = parse(&general);
    assert_eq!(crystal.atoms.len(), 2);
    assert_close(crystal.atoms[1].fractional_pos, [0.9, 0.8, 0.7]);
}

#[test]
fn monoclinic_screw_and_glide_operations_wrap_into_the_cell() {
    let cif = format!(
        "data_p21c\n{cell}\
         loop_\n_symmetry_equiv_pos_as_xyz\n\
         'x, y, z'\n'-x, y+1/2, -z+1/2'\n'-x, -y, -z'\n'x, -y+1/2, z+1/2'\n\
         {SITE_COLUMNS}Fe1 Fe 0.1 0.2 0.3\n",
        cell = cubic_cell(6.0)
    );
    let crystal = parse(&cif);
    let expected = [
        [0.1, 0.2, 0.3],
        [0.9, 0.7, 0.2],
        [0.9, 0.8, 0.7],
        [0.1, 0.3, 0.8],
    ];
    assert_eq!(crystal.atoms.len(), 4);
    for (atom, want) in crystal.atoms.iter().zip(expected) {
        assert_close(atom.fractional_pos, want);
    }
}

#[test]
fn hexagonal_operations_with_mixed_axes_and_a_special_position() {
    let cell = "_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 6.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\n";
    let ops = "loop_\n_space_group_symop_operation_xyz\nx,y,z\n-y,x-y,z\n-x+y,-x,z\n";
    let special = format!("data_p3\n{cell}{ops}{SITE_COLUMNS}Mg1 Mg 0.333333333 0.666666667 0.0\n");
    assert_eq!(parse(&special).atoms.len(), 1);
    let general = format!("data_p3\n{cell}{ops}{SITE_COLUMNS}Mg1 Mg 0.1 0.2 0.3\n");
    let crystal = parse(&general);
    let expected = [[0.1, 0.2, 0.3], [0.8, 0.9, 0.3], [0.1, 0.9, 0.3]];
    assert_eq!(crystal.atoms.len(), 3);
    for (atom, want) in crystal.atoms.iter().zip(expected) {
        assert_close(atom.fractional_pos, want);
    }
}

#[test]
fn rutile_expands_to_two_titanium_and_four_oxygen_atoms() {
    // P4_2/mnm (No. 136): 16 operations, Ti on 2a (0,0,0), O on 4f (x,x,0).
    let ops = [
        "x,y,z",
        "-x,-y,z",
        "-y+1/2,x+1/2,z+1/2",
        "y+1/2,-x+1/2,z+1/2",
        "-x+1/2,y+1/2,-z+1/2",
        "x+1/2,-y+1/2,-z+1/2",
        "y,x,-z",
        "-y,-x,-z",
        "-x,-y,-z",
        "x,y,-z",
        "y+1/2,-x+1/2,-z+1/2",
        "-y+1/2,x+1/2,-z+1/2",
        "x+1/2,-y+1/2,z+1/2",
        "-x+1/2,y+1/2,z+1/2",
        "-y,-x,z",
        "y,x,z",
    ];
    let mut cif = String::from(
        "data_rutile\n_cell_length_a 4.5937\n_cell_length_b 4.5937\n_cell_length_c 2.9587\n\
         _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
         _symmetry_Int_Tables_number 136\nloop_\n_symmetry_equiv_pos_as_xyz\n",
    );
    for op in ops {
        cif.push_str(&format!("'{op}'\n"));
    }
    cif.push_str(SITE_COLUMNS);
    cif.push_str("Ti1 Ti4+ 0 0 0\nO1 O2- 0.3053 0.3053 0\n");
    let crystal = parse(&cif);
    assert_eq!(crystal.space_group, 136);
    let positions_of = |element: &str| -> Vec<[f64; 3]> {
        crystal
            .atoms
            .iter()
            .filter(|a| a.element == element)
            .map(|a| a.fractional_pos)
            .collect()
    };
    let ti = positions_of("Ti");
    assert_eq!(ti.len(), 2, "{ti:?}");
    assert_close(ti[0], [0.0, 0.0, 0.0]);
    assert_close(ti[1], [0.5, 0.5, 0.5]);
    let o = positions_of("O");
    assert_eq!(o.len(), 4, "{o:?}");
    let expected = [
        [0.3053, 0.3053, 0.0],
        [0.6947, 0.6947, 0.0],
        [0.1947, 0.8053, 0.5],
        [0.8053, 0.1947, 0.5],
    ];
    for want in expected {
        assert!(
            o.iter().any(|got| got
                .iter()
                .zip(want.iter())
                .all(|(g, w)| (g - w).abs() < 1e-9)),
            "missing O at {want:?} in {o:?}"
        );
    }
    assert_eq!(crystal.atoms.len(), 6);
}

#[test]
fn files_without_symmetry_operations_keep_their_atoms_as_listed() {
    let cif = format!(
        "data_p1\n{cell}{SITE_COLUMNS}Fe1 Fe 1.2 -0.3 0.5\nFe2 Fe 0.5 0.5 0.5\n",
        cell = cubic_cell(4.0)
    );
    let crystal = parse(&cif);
    assert_eq!(crystal.atoms.len(), 2);
    assert_close(crystal.atoms[0].fractional_pos, [1.2, -0.3, 0.5]);
}

#[test]
fn an_invalid_symmetry_operation_fails_the_parse() {
    let cif = format!(
        "data_bad\n{cell}loop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n'x, y'\n{SITE_COLUMNS}Fe1 Fe 0 0 0\n",
        cell = cubic_cell(4.0)
    );
    let err = CifParser::new().parse(&cif).unwrap_err();
    assert!(err.to_string().contains("symmetry operation"), "{err}");
}

#[test]
fn unreadable_atom_data_fails_the_parse_instead_of_moving_atoms_to_the_origin() {
    for row in ["Fe1 Fe ? 0 0", "Fe1 Fe 0.1 abc 0", "Fe1 Fe 0.1 0.2"] {
        let cif = format!(
            "data_bad\n{cell}{SITE_COLUMNS}{row}\n",
            cell = cubic_cell(4.0)
        );
        assert!(CifParser::new().parse(&cif).is_err(), "{row}");
    }
}

#[test]
fn only_the_first_data_block_is_read() {
    let cif = format!(
        "data_first\n{cell}{SITE_COLUMNS}Fe1 Fe 0 0 0\n\
         data_second\n_cell_length_a 99\n{SITE_COLUMNS}Ni1 Ni 0.5 0.5 0.5\n",
        cell = cubic_cell(4.0)
    );
    let crystal = parse(&cif);
    assert_eq!(crystal.atoms.len(), 1);
    assert_eq!(crystal.atoms[0].element, "Fe");
    assert!((crystal.lattice[0][0] - 4.0).abs() < 1e-12);
}

#[test]
fn comments_and_text_fields_do_not_disturb_the_tables() {
    let cif = format!(
        "# a comment\ndata_text\n_publ_section_title\n;\nloop_\n_atom_site_label\nNot a table\n;\n\
         {cell}# cell done\n{SITE_COLUMNS}\
         Fe1 Fe 0.0 0.0 0.0 # first site\n\
         Fe2 Fe 0.5 0.5 0.5\n",
        cell = cubic_cell(4.0)
    );
    let crystal = parse(&cif);
    assert_eq!(crystal.atoms.len(), 2);
    assert_close(crystal.atoms[1].fractional_pos, [0.5, 0.5, 0.5]);
}

#[test]
fn cartesian_atom_sites_are_converted_to_fractional_coordinates() {
    let cif = "data_hex\n_cell_length_a 3.0\n_cell_length_b 3.0\n_cell_length_c 5.0\n\
               _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\n\
               loop_\n_atom_site_label\n_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\n\
               Mg1 0.0 1.7320508075688772 2.5\n";
    let crystal = parse(cif);
    assert_close(crystal.atoms[0].fractional_pos, [1.0 / 3.0, 2.0 / 3.0, 0.5]);
}

#[test]
fn write_cif_output_reads_back_to_the_same_structure() {
    let mut crystal = CrystalStructure::new_cubic(5.5);
    crystal.space_group = 225;
    crystal.atoms.push(CrystalAtom {
        element: "Na".into(),
        fractional_pos: [0.0, 0.0, 0.0],
        occupancy: 1.0,
        b_factor: 0.5,
    });
    crystal.atoms.push(CrystalAtom {
        element: "Cl".into(),
        fractional_pos: [0.5, 0.25, 0.75],
        occupancy: 0.9,
        b_factor: 1.25,
    });
    let back = parse(&write_cif(&crystal));
    assert_eq!(back.space_group, 225);
    assert!((back.volume() - crystal.volume()).abs() < 1e-6);
    assert_eq!(back.atoms.len(), 2);
    for (got, want) in back.atoms.iter().zip(&crystal.atoms) {
        assert_eq!(got.element, want.element);
        assert_close(got.fractional_pos, want.fractional_pos);
        assert!((got.occupancy - want.occupancy).abs() < 1e-9);
        assert!((got.b_factor - want.b_factor).abs() < 1e-3);
    }
}
