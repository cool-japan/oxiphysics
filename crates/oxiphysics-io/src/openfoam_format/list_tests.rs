// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The simplified OpenFOAM readers must take the cell values from the
//! `internalField` entry only: parentheses in the `boundaryField` section may
//! neither add cells nor break the parse.

use super::{FoamField, parse_foam_scalar, parse_foam_vector, read_foam_mesh, write_foam_field};

const FOAM_HEADER: &str = "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       volVectorField;\n    object      U;\n}\n";

/// A `U` file with `values` as a nonuniform internal field and `n_patches`
/// boundary patches that each carry `value uniform (...)`.
fn vector_file(values: &[[f64; 3]], n_patches: usize) -> String {
    let mut s = String::from(FOAM_HEADER);
    s.push_str("dimensions      [0 1 -1 0 0 0 0];\n\n");
    s.push_str("internalField   nonuniform List<vector>\n");
    s.push_str(&format!("{}\n(\n", values.len()));
    for v in values {
        s.push_str(&format!("({} {} {})\n", v[0], v[1], v[2]));
    }
    s.push_str(")\n;\n\nboundaryField\n{\n");
    for p in 0..n_patches {
        s.push_str(&format!(
            "    patch{p}\n    {{\n        type            fixedValue;\n        value           uniform ({p} 0 0);\n    }}\n"
        ));
    }
    s.push_str("}\n");
    s
}

fn scalar_file(internal: &str, boundary_value: &str) -> String {
    format!(
        "dimensions      [0 0 0 0 0 0 0];\n\ninternalField   {internal};\n\nboundaryField\n{{\n    inlet\n    {{\n        type            fixedValue;\n        value           {boundary_value};\n    }}\n    outlet\n    {{\n        type            zeroGradient;\n    }}\n}}\n"
    )
}

#[test]
fn two_boundary_patches_with_uniform_vectors_add_no_cells() {
    let cells = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
    let field = parse_foam_vector(&vector_file(&cells, 2), "U").unwrap();
    assert_eq!(field.n_cells(), 3);
    assert_eq!(
        field.internal_field,
        vec![1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]
    );
}

#[test]
fn any_number_of_boundary_patches_leaves_the_cell_count_alone() {
    let cells: Vec<[f64; 3]> = (0..7).map(|i| [i as f64, -(i as f64), 0.5]).collect();
    for n_patches in [0usize, 1, 2, 3, 12, 60] {
        let field = parse_foam_vector(&vector_file(&cells, n_patches), "U").unwrap();
        assert_eq!(field.n_cells(), 7, "{n_patches} patches");
        assert_eq!(field.internal_field.len(), 21);
        assert_eq!(&field.internal_field[18..], &[6.0, -6.0, 0.5]);
    }
}

#[test]
fn boundary_patches_are_still_reported() {
    let cells = [[1.0, 2.0, 3.0]];
    let field = parse_foam_vector(&vector_file(&cells, 2), "U").unwrap();
    assert_eq!(field.patches.len(), 2);
    assert_eq!(field.patches[1].name, "patch1");
    assert_eq!(field.patches[1].patch_type, "fixedValue");
}

#[test]
fn parentheses_in_the_boundary_field_do_not_break_a_scalar_list() {
    let src = "dimensions [0 0 0 0 0 0 0];\ninternalField nonuniform List<scalar>\n3\n(\n1.5\n2.5\n3.5\n)\n;\nboundaryField\n{\n    wall\n    {\n        type codedFixedValue;\n        code #{ value = (a + (b * c)); #};\n        value uniform 0;\n    }\n}\n";
    let field = parse_foam_scalar(src, "p").unwrap();
    assert_eq!(field.internal_field, vec![1.5, 2.5, 3.5]);
}

#[test]
fn uniform_internal_field_is_not_hijacked_by_a_boundary_list() {
    let src = format!(
        "// nCells: 4\n{}",
        scalar_file("uniform 7.5", "nonuniform List<scalar> 2(1 2)")
    );
    let field = parse_foam_scalar(&src, "T").unwrap();
    assert_eq!(field.internal_field, vec![7.5; 4]);
}

#[test]
fn uniform_vector_internal_field_is_expanded() {
    let src = "// nCells: 3\ninternalField   uniform (1 2 3);\nboundaryField\n{\n    wall\n    {\n        type fixedValue;\n        value uniform (0 0 0);\n    }\n}\n";
    let field = parse_foam_vector(src, "U").unwrap();
    assert_eq!(field.n_cells(), 3);
    assert_eq!(
        field.internal_field,
        vec![1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
    );
    assert_eq!(field.n_components, 3);
}

#[test]
fn uniform_detection_tolerates_any_whitespace() {
    for internal in ["uniform 5", "uniform\t5", "uniform\n5"] {
        let field = parse_foam_scalar(&scalar_file(internal, "uniform 0"), "p").unwrap();
        assert_eq!(field.internal_field, vec![5.0], "{internal:?}");
    }
    let src = "internalField\tuniform\t(4 5 6);\n";
    assert_eq!(
        parse_foam_vector(src, "U").unwrap().internal_field,
        vec![4.0, 5.0, 6.0]
    );
}

#[test]
fn compact_size_glued_to_the_parenthesis_is_read() {
    let scalar = scalar_file("nonuniform List<scalar> 3(1 2 3)", "uniform 0");
    assert_eq!(
        parse_foam_scalar(&scalar, "p").unwrap().internal_field,
        vec![1.0, 2.0, 3.0]
    );
    let vector = "internalField nonuniform List<vector> 2((1 2 3) (4 5 6));\nboundaryField\n{\n    a\n    {\n        value uniform (9 9 9);\n    }\n    b\n    {\n        value uniform (8 8 8);\n    }\n}\n";
    let field = parse_foam_vector(vector, "U").unwrap();
    assert_eq!(field.n_cells(), 2);
    assert_eq!(field.internal_field, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
}

#[test]
fn declared_size_must_match_the_entries() {
    let scalar = scalar_file("nonuniform List<scalar> 4(1 2 3)", "uniform 0");
    let err = parse_foam_scalar(&scalar, "p").unwrap_err().to_string();
    assert!(
        err.contains("declares 4") && err.contains("holds 3"),
        "{err}"
    );
    let too_many = scalar_file("nonuniform List<scalar> 2(1 2 3)", "uniform 0");
    assert!(parse_foam_scalar(&too_many, "p").is_err());

    let vector = "internalField nonuniform List<vector> 2((1 2 3) (4 5 6) (7 8 9));\n";
    let err = parse_foam_vector(vector, "U").unwrap_err().to_string();
    assert!(
        err.contains("declares 2") && err.contains("holds 3"),
        "{err}"
    );
}

#[test]
fn empty_lists_are_valid() {
    let scalar = scalar_file("nonuniform List<scalar> 0()", "uniform 0");
    assert!(
        parse_foam_scalar(&scalar, "p")
            .unwrap()
            .internal_field
            .is_empty()
    );
    let vector = "internalField nonuniform 0();\n";
    assert_eq!(parse_foam_vector(vector, "U").unwrap().n_cells(), 0);
}

#[test]
fn mentions_of_the_keyword_elsewhere_are_ignored() {
    let src = "// internalField (a comment with (parentheses))\n/* internalField uniform 99; */\ninternalField   nonuniform List<scalar> 2(1 2);\nboundaryField\n{\n    a\n    {\n        value $internalField;\n    }\n}\n";
    let field = parse_foam_scalar(src, "p").unwrap();
    assert_eq!(field.internal_field, vec![1.0, 2.0]);
}

#[test]
fn nonuniform_without_a_list_does_not_borrow_one_from_the_boundary() {
    let src = scalar_file("nonuniform List<scalar>", "nonuniform List<scalar> 1(5)");
    assert!(parse_foam_scalar(&src, "p").is_err());
    let src = "internalField nonuniform List<vector>;\nboundaryField\n{\n    a\n    {\n        value uniform (1 2 3);\n    }\n}\n";
    assert!(parse_foam_vector(src, "U").is_err());
}

#[test]
fn wrong_list_type_and_unknown_forms_are_errors() {
    let src = "internalField nonuniform List<vector> 1((1 2 3));\n";
    assert!(parse_foam_scalar(src, "p").is_err());
    let src = "internalField nonuniform List<scalar> 1(1);\n";
    assert!(parse_foam_vector(src, "U").is_err());
    assert!(parse_foam_scalar("internalField $internal;\n", "p").is_err());
    assert!(parse_foam_scalar("internalField uniform;\n", "p").is_err());
    assert!(parse_foam_vector("internalField uniform 5;\n", "U").is_err());
    assert!(parse_foam_vector("internalField uniform (1 2);\n", "U").is_err());
}

#[test]
fn malformed_vectors_are_errors() {
    for body in ["(1 2)", "(1 2 3 4)", "(1 two 3)", "(1 2 3"] {
        let src = format!("internalField nonuniform List<vector> 1({body});\n");
        assert!(parse_foam_vector(&src, "U").is_err(), "{body}");
    }
    let src = "internalField nonuniform List<vector> 2((1 2 3) 4 5 6);\n";
    assert!(parse_foam_vector(src, "U").is_err());
}

#[test]
fn unbalanced_parentheses_are_errors() {
    assert!(parse_foam_scalar("internalField nonuniform List<scalar> 2(1 2;\n", "p").is_err());
    assert!(parse_foam_vector("internalField nonuniform List<vector> 1((1 2 3);\n", "U").is_err());
}

#[test]
fn absurd_cell_count_comment_is_rejected_instead_of_allocated() {
    let src = format!(
        "// nCells: 99999999999999\n{}",
        scalar_file("uniform 1", "uniform 0")
    );
    assert!(parse_foam_scalar(&src, "p").is_err());
    let src = "// nCells: 99999999999999\ninternalField uniform (1 2 3);\n";
    assert!(parse_foam_vector(src, "U").is_err());
}

#[test]
fn bare_list_snippets_without_the_keyword_still_parse() {
    let field = parse_foam_scalar("List<scalar> 2(4 5)", "p").unwrap();
    assert_eq!(field.internal_field, vec![4.0, 5.0]);
    let field = parse_foam_vector("List<vector> 1((1 2 3))", "U").unwrap();
    assert_eq!(field.internal_field, vec![1.0, 2.0, 3.0]);
}

#[test]
fn written_fields_read_back_with_many_patches() {
    let mut field = FoamField::new_vector("U", vec![1.5, -2.5, 0.125, 4.0, 5.0, 6.0]);
    field.patches = (0..5)
        .map(|i| super::foam_boundary_patch(format!("p{i}"), "fixedValue", Some(i as f64)))
        .collect();
    let text = write_foam_field(&field);
    let back = parse_foam_vector(&text, "U").unwrap();
    assert_eq!(back.internal_field, field.internal_field);
    assert_eq!(back.patches.len(), 5);
    assert_eq!(back.patches[3].uniform_value, Some(3.0));
}

const MESH_HEADER: &str = "/*--------------------------------*- C++ -*----------------------------------*\\\n| =========                 |                                                 |\n| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\\*---------------------------------------------------------------------------*/\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       vectorField;\n    object      points;\n}\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

#[test]
fn polymesh_files_with_banner_headers_are_read() {
    let points = format!("{MESH_HEADER}3\n(\n(0 0 0)\n(1 0 0)\n(0.5 1 0)\n)\n\n// ***** //\n");
    let faces = format!("{MESH_HEADER}1\n(\n3(0 1 2)\n)\n\n// ***** //\n");
    let owner = format!("{MESH_HEADER}1\n(\n0\n)\n\n// ***** //\n");
    let neighbour = format!("{MESH_HEADER}0\n(\n)\n\n// ***** //\n");
    let mesh = read_foam_mesh(&points, &faces, &owner, &neighbour).unwrap();
    assert_eq!(mesh.n_points(), 3);
    assert_eq!(mesh.faces, vec![vec![0, 1, 2]]);
    assert_eq!(mesh.owner, vec![0]);
    assert_eq!(mesh.n_internal_faces(), 0);
}

#[test]
fn polymesh_declared_sizes_are_checked() {
    let faces = "1\n(\n3(0 1 2)\n)\n";
    let owner = "1\n(\n0\n)\n";
    let neighbour = "0\n(\n)\n";
    let two_points = "3\n(\n(0 0 0)\n(1 0 0)\n)\n";
    assert!(read_foam_mesh(two_points, faces, owner, neighbour).is_err());
    let points = "3\n(\n(0 0 0)\n(1 0 0)\n(0.5 1 0)\n)\n";
    let short_face = "1\n(\n4(0 1 2)\n)\n";
    assert!(read_foam_mesh(points, short_face, owner, neighbour).is_err());
    let two_owners = "2\n(\n0\n)\n";
    assert!(read_foam_mesh(points, faces, two_owners, neighbour).is_err());
    let junk_in_faces = "1\n(\nface 3(0 1 2)\n)\n";
    assert!(read_foam_mesh(points, junk_in_faces, owner, neighbour).is_err());
    let negative = "1\n(\n-1\n)\n";
    assert!(read_foam_mesh(points, faces, owner, negative).is_err());
}

#[test]
fn trailing_text_after_the_list_is_not_swallowed() {
    // The old rfind-based lookup extended the list up to the last ')' of the file.
    let points = "2\n(\n(0 0 0)\n(1 1 1)\n)\n// note (see the manual)\n";
    let faces = "0\n(\n)\n";
    let owner = "0\n(\n)\n";
    let neighbour = "0\n(\n)\n";
    let mesh = read_foam_mesh(points, faces, owner, neighbour).unwrap();
    assert_eq!(mesh.n_points(), 2);
}
