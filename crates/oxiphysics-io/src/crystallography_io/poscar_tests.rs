// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the VASP POSCAR reader and writer: Cartesian
//! coordinates on skewed cells, the universal scale factor (including the
//! negative cell-volume convention) and the coordinate-mode label.

use super::*;

const HEX_A: f64 = 3.2;
const HEX_C: f64 = 5.2;

/// A one-atom hexagonal POSCAR (gamma = 120 degrees) with the given scale line,
/// coordinate-mode line and atom position.
fn hexagonal_poscar(scale: &str, mode: &str, position: [f64; 3]) -> String {
    format!(
        "Hexagonal Mg\n{scale}\n  {a:.12} 0.000000000000 0.000000000000\n  {bx:.12} {by:.12} 0.000000000000\n  0.000000000000 0.000000000000 {c:.12}\nMg\n1\n{mode}\n  {x:.12} {y:.12} {z:.12}\n",
        a = HEX_A,
        bx = -HEX_A / 2.0,
        by = HEX_A * 3.0_f64.sqrt() / 2.0,
        c = HEX_C,
        x = position[0],
        y = position[1],
        z = position[2],
    )
}

fn assert_close(got: [f64; 3], want: [f64; 3], tol: f64) {
    for (g, w) in got.iter().zip(want.iter()) {
        assert!((g - w).abs() < tol, "{got:?} vs {want:?}");
    }
}

fn volume(lattice: &[[f64; 3]; 3]) -> f64 {
    lattice_math::determinant(lattice).abs()
}

#[test]
fn cartesian_input_on_a_hexagonal_cell_uses_the_full_lattice() {
    let reader = VaspReader::new();
    // r = 1/3 a1 + 2/3 a2 + 1/2 a3 = (0, a / sqrt(3), c / 2).
    let cartesian = hexagonal_poscar(
        "1.0",
        "Cartesian",
        [0.0, HEX_A / 3.0_f64.sqrt(), HEX_C / 2.0],
    );
    let direct = hexagonal_poscar("1.0", "Direct", [1.0 / 3.0, 2.0 / 3.0, 0.5]);
    let from_cartesian = reader.parse_poscar_full(&cartesian).unwrap();
    let from_direct = reader.parse_poscar_full(&direct).unwrap();
    assert!(from_cartesian.is_direct);
    assert!(from_direct.is_direct);
    assert_close(
        from_cartesian.positions[0],
        [1.0 / 3.0, 2.0 / 3.0, 0.5],
        1e-10,
    );
    assert_close(from_cartesian.positions[0], from_direct.positions[0], 1e-10);
    // The diagonal-only division of the old reader gave x = 0, not 1/3.
    assert!((from_cartesian.positions[0][0] - 1.0 / 3.0).abs() < 1e-10);
}

#[test]
fn parse_poscar_returns_fractional_atoms_for_cartesian_input() {
    let cartesian = hexagonal_poscar(
        "1.0",
        "cartesian",
        [0.0, HEX_A / 3.0_f64.sqrt(), HEX_C / 2.0],
    );
    let crystal = VaspReader::new().parse_poscar(&cartesian).unwrap();
    assert_close(
        crystal.atoms[0].fractional_pos,
        [1.0 / 3.0, 2.0 / 3.0, 0.5],
        1e-10,
    );
    // Converting back reproduces the Cartesian position.
    let cart = crystal.to_cartesian()[0];
    assert_close(cart, [0.0, HEX_A / 3.0_f64.sqrt(), HEX_C / 2.0], 1e-10);
}

#[test]
fn kartesian_spelling_is_also_cartesian() {
    let text = hexagonal_poscar(
        "1.0",
        "Kartesian",
        [0.0, HEX_A / 3.0_f64.sqrt(), HEX_C / 2.0],
    );
    let ps = VaspReader::new().parse_poscar_full(&text).unwrap();
    assert_close(ps.positions[0], [1.0 / 3.0, 2.0 / 3.0, 0.5], 1e-10);
}

#[test]
fn cartesian_coordinates_are_multiplied_by_the_universal_scale() {
    // Unit cube scaled by 2: the cell edge is 2 Å and a Cartesian position
    // (0.5, 0.25, 0.75) in file units is (1.0, 0.5, 1.5) Å, i.e. fractional
    // (0.5, 0.25, 0.75).
    let text = "Scaled\n2.0\n1 0 0\n0 1 0\n0 0 1\nH\n1\nCartesian\n0.5 0.25 0.75\n";
    let ps = VaspReader::new().parse_poscar_full(text).unwrap();
    assert!((ps.scale - 2.0).abs() < 1e-15);
    assert!((ps.lattice[0][0] - 2.0).abs() < 1e-15);
    assert_close(ps.positions[0], [0.5, 0.25, 0.75], 1e-14);

    // Direct coordinates are not touched by the scale.
    let direct = "Scaled\n2.0\n1 0 0\n0 1 0\n0 0 1\nH\n1\nDirect\n0.5 0.25 0.75\n";
    let ps = VaspReader::new().parse_poscar_full(direct).unwrap();
    assert_close(ps.positions[0], [0.5, 0.25, 0.75], 1e-15);
}

#[test]
fn a_negative_scale_is_the_cell_volume() {
    // Unit cube with scale -27: the cell is stretched to 27 Å^3, an edge of 3 Å.
    let text = "Volume\n-27.0\n1 0 0\n0 1 0\n0 0 1\nH\n1\nCartesian\n0.5 0.5 0.5\n";
    let ps = VaspReader::new().parse_poscar_full(text).unwrap();
    assert!((ps.scale - 3.0).abs() < 1e-12, "scale {}", ps.scale);
    assert!((ps.lattice[0][0] - 3.0).abs() < 1e-12);
    assert!((volume(&ps.lattice) - 27.0).abs() < 1e-10);
    // 0.5 file units * 3 = 1.5 Å = half of the 3 Å edge.
    assert_close(ps.positions[0], [0.5, 0.5, 0.5], 1e-12);
    // A negative factor must never invert the lattice.
    assert!(ps.lattice[0][0] > 0.0 && ps.lattice[1][1] > 0.0 && ps.lattice[2][2] > 0.0);

    // Non-cubic cell: diag(2, 3, 4) has volume 24; -48 doubles it.
    let ortho = "Volume\n-48\n2 0 0\n0 3 0\n0 0 4\nH\n1\nDirect\n0 0 0\n";
    let ps = VaspReader::new().parse_poscar_full(ortho).unwrap();
    assert!((ps.scale - 2.0_f64.cbrt()).abs() < 1e-12);
    assert!((volume(&ps.lattice) - 48.0).abs() < 1e-9);
    let crystal = VaspReader::new().parse_poscar(ortho).unwrap();
    assert!((crystal.volume() - 48.0).abs() < 1e-9);
}

#[test]
fn degenerate_scales_and_lattices_are_errors() {
    let reader = VaspReader::new();
    for scale in ["0.0", "-0.0", "nan", "inf", "-inf", "abc"] {
        let text = format!("T\n{scale}\n1 0 0\n0 1 0\n0 0 1\nH\n1\nDirect\n0 0 0\n");
        assert!(reader.parse_poscar_full(&text).is_err(), "scale {scale}");
    }
    let singular_volume = "T\n-10\n1 0 0\n2 0 0\n0 0 1\nH\n1\nDirect\n0 0 0\n";
    let err = reader.parse_poscar_full(singular_volume).unwrap_err();
    assert!(err.to_string().contains("singular"), "{err}");
    let singular_cartesian = "T\n1.0\n1 0 0\n2 0 0\n0 0 1\nH\n1\nCartesian\n0 0 0\n";
    let err = reader.parse_poscar_full(singular_cartesian).unwrap_err();
    assert!(err.to_string().contains("singular"), "{err}");
}

#[test]
fn written_poscar_is_labelled_direct_and_round_trips() {
    let reader = VaspReader::new();
    let cartesian = hexagonal_poscar(
        "1.0",
        "Cartesian",
        [0.0, HEX_A / 3.0_f64.sqrt(), HEX_C / 2.0],
    );
    let parsed = reader.parse_poscar_full(&cartesian).unwrap();
    let written = parsed.to_poscar_string();
    // Fractional values must not be announced as Cartesian.
    assert!(written.contains("Direct"), "{written}");
    assert!(!written.contains("Cartesian"), "{written}");
    let again = reader.parse_poscar_full(&written).unwrap();
    assert!(again.is_direct);
    assert_close(again.positions[0], parsed.positions[0], 1e-8);
    for (row_a, row_b) in again.lattice.iter().zip(parsed.lattice.iter()) {
        assert_close(*row_a, *row_b, 1e-8);
    }
}

#[test]
fn a_hand_built_cartesian_structure_is_converted_and_written_as_cartesian() {
    let lattice = [[4.0, 0.0, 0.0], [1.5, 3.5, 0.0], [0.5, 0.75, 5.0]];
    let fractional = [0.2, 0.4, 0.6];
    let cartesian = lattice_math::row_times_matrix(fractional, &lattice);
    let structure = PoscarStructure {
        comment: "hand built".into(),
        scale: 1.0,
        lattice,
        species: vec!["Fe".into()],
        counts: vec![1],
        selective_dynamics: false,
        is_direct: false,
        positions: vec![cartesian],
        sd_flags: None,
        magmom: None,
        forces: None,
    };
    let frac = structure.fractional_positions().unwrap();
    assert_close(frac[0], fractional, 1e-12);
    let crystal = structure.to_crystal_structure();
    assert_close(crystal.atoms[0].fractional_pos, fractional, 1e-12);

    let written = structure.to_poscar_string();
    assert!(written.contains("Cartesian"));
    let reparsed = VaspReader::new().parse_poscar_full(&written).unwrap();
    assert_close(reparsed.positions[0], fractional, 1e-8);

    // A direct-mode structure passes its positions through untouched.
    let direct = PoscarStructure {
        is_direct: true,
        positions: vec![fractional],
        ..structure.clone()
    };
    assert_eq!(direct.fractional_positions().unwrap(), vec![fractional]);

    // A singular lattice has no Cartesian conversion: the strict accessor
    // reports it and the lenient conversion keeps the stored numbers.
    let flat = PoscarStructure {
        lattice: [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        ..structure
    };
    assert!(flat.fractional_positions().is_err());
    assert_eq!(
        flat.to_crystal_structure().atoms[0].fractional_pos,
        cartesian
    );
}

#[test]
fn selective_dynamics_flags_survive_a_cartesian_conversion() {
    // Skewed cell: the second atom sits at fractional (0.5, 0.5, 0.5), which is
    // Cartesian (3.0, 2.125, 2.5) here.
    let text = "SD\n1.0\n4 0 0\n1.5 3.5 0\n0.5 0.75 5\nFe\n2\nSelective dynamics\nCartesian\n0 0 0 T T T\n3.0 2.125 2.5 F F T\n";
    let ps = VaspReader::new().parse_poscar_full(text).unwrap();
    assert!(ps.selective_dynamics);
    let flags = ps.sd_flags.as_ref().expect("flags");
    assert_eq!(flags[0], [true, true, true]);
    assert_eq!(flags[1], [false, false, true]);
    assert_close(ps.positions[0], [0.0, 0.0, 0.0], 1e-14);
    assert_close(ps.positions[1], [0.5, 0.5, 0.5], 1e-12);
}
