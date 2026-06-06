//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    FieldValues, FoamBc, FoamDict, FoamField, FoamPatch, FoamResidual, FoamTimeDir, FoamValue,
};

/// Returns the standard FoamFile header dict for an OpenFOAM file.
pub fn foam_header(class: &str, object: &str) -> String {
    format!(
        "/*--------------------------------*- C++ -*----------------------------------*\\\n\
         | =========                 |                                                 |\n\
         | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
         |  \\\\    /   O peration     | Version:  v2312                                 |\n\
         |   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n\
         |    \\\\/     M anipulation  |                                                 |\n\
         \\*---------------------------------------------------------------------------*/\n\
         FoamFile\n\
         {{\n\
             version     2.0;\n\
             format      ascii;\n\
             class       {class};\n\
             object      {object};\n\
         }}\n\
         // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )
}
/// Strip C and C++ style comments from input.
pub fn strip_foam_comments(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    let chars: Vec<char> = input.chars().collect();
    let n = chars.len();
    let mut i = 0;
    while i < n {
        if i + 1 < n && chars[i] == '/' && chars[i + 1] == '/' {
            while i < n && chars[i] != '\n' {
                i += 1;
            }
        } else if i + 1 < n && chars[i] == '/' && chars[i + 1] == '*' {
            i += 2;
            while i + 1 < n && !(chars[i] == '*' && chars[i + 1] == '/') {
                i += 1;
            }
            if i + 1 < n {
                i += 2;
            }
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}
/// Tokenise OpenFOAM dictionary text into simple tokens.
pub fn tokenise_foam(input: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let n = chars.len();
    let mut i = 0;
    while i < n {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        if c == '{' || c == '}' || c == '(' || c == ')' || c == ';' {
            tokens.push(c.to_string());
            i += 1;
            continue;
        }
        if c == '"' {
            let mut tok = String::new();
            i += 1;
            while i < n && chars[i] != '"' {
                tok.push(chars[i]);
                i += 1;
            }
            if i < n {
                i += 1;
            }
            tokens.push(tok);
            continue;
        }
        let mut tok = String::new();
        while i < n
            && !chars[i].is_whitespace()
            && chars[i] != '{'
            && chars[i] != '}'
            && chars[i] != '('
            && chars[i] != ')'
            && chars[i] != ';'
        {
            tok.push(chars[i]);
            i += 1;
        }
        if !tok.is_empty() {
            tokens.push(tok);
        }
    }
    tokens
}
/// Parse tokens into a `FoamDict`, starting at position `pos`.
/// Returns the dict and the next token position.
pub fn parse_dict_tokens(tokens: &[String], mut pos: usize) -> (FoamDict, usize) {
    let mut dict = FoamDict::new();
    while pos < tokens.len() {
        let tok = &tokens[pos];
        if tok == "}" {
            pos += 1;
            break;
        }
        if tok == ";" || tok == "{" || tok == "(" || tok == ")" {
            pos += 1;
            continue;
        }
        let key = tok.clone();
        pos += 1;
        if pos >= tokens.len() {
            break;
        }
        if tokens[pos] == "{" {
            pos += 1;
            let (sub, next) = parse_dict_tokens(tokens, pos);
            dict.insert(key, FoamValue::Dict(sub));
            pos = next;
            continue;
        }
        if tokens[pos] == "(" {
            pos += 1;
            let mut items = Vec::new();
            while pos < tokens.len() && tokens[pos] != ")" {
                if let Ok(v) = tokens[pos].parse::<f64>() {
                    items.push(FoamValue::Scalar(v));
                } else {
                    items.push(FoamValue::Word(tokens[pos].clone()));
                }
                pos += 1;
            }
            if pos < tokens.len() {
                pos += 1;
            }
            if pos < tokens.len() && tokens[pos] == ";" {
                pos += 1;
            }
            if items.len() == 3 && items.iter().all(|x| matches!(x, FoamValue::Scalar(_))) {
                let mut v = [0.0; 3];
                for (idx, item) in items.iter().enumerate() {
                    if let FoamValue::Scalar(s) = item {
                        v[idx] = *s;
                    }
                }
                dict.insert(key, FoamValue::Vector(v));
            } else {
                dict.insert(key, FoamValue::List(items));
            }
            continue;
        }
        let val_tok = tokens[pos].clone();
        pos += 1;
        if pos < tokens.len() && tokens[pos] == ";" {
            pos += 1;
        }
        if let Ok(v) = val_tok.parse::<f64>() {
            dict.insert(key, FoamValue::Scalar(v));
        } else {
            dict.insert(key, FoamValue::Word(val_tok));
        }
    }
    (dict, pos)
}
/// Parse an OpenFOAM `points` file content into a vector of 3D coordinates.
pub fn parse_foam_points(input: &str) -> Vec<[f64; 3]> {
    let cleaned = strip_foam_comments(input);
    let mut points = Vec::new();
    let mut in_list = false;
    for line in cleaned.lines() {
        let trimmed = line.trim();
        if trimmed == "(" {
            in_list = true;
            continue;
        }
        if trimmed == ")" {
            break;
        }
        if !in_list {
            continue;
        }
        let inner = trimmed.trim_start_matches('(').trim_end_matches(')');
        let parts: Vec<&str> = inner.split_whitespace().collect();
        if parts.len() == 3
            && let (Ok(x), Ok(y), Ok(z)) = (
                parts[0].parse::<f64>(),
                parts[1].parse::<f64>(),
                parts[2].parse::<f64>(),
            )
        {
            points.push([x, y, z]);
        }
    }
    points
}
/// Parse an OpenFOAM `owner` or `neighbour` file content into a list of integers.
pub fn parse_foam_label_list(input: &str) -> Vec<i64> {
    let cleaned = strip_foam_comments(input);
    let mut labels = Vec::new();
    let mut in_list = false;
    for line in cleaned.lines() {
        let trimmed = line.trim();
        if trimmed == "(" {
            in_list = true;
            continue;
        }
        if trimmed == ")" {
            break;
        }
        if !in_list {
            continue;
        }
        if let Ok(v) = trimmed.parse::<i64>() {
            labels.push(v);
        }
    }
    labels
}
/// Parse an OpenFOAM `faces` file content into face vertex lists.
pub fn parse_foam_faces(input: &str) -> Vec<Vec<usize>> {
    let cleaned = strip_foam_comments(input);
    let mut faces = Vec::new();
    let mut in_list = false;
    for line in cleaned.lines() {
        let trimmed = line.trim();
        if trimmed == "(" {
            in_list = true;
            continue;
        }
        if trimmed == ")" {
            break;
        }
        if !in_list {
            continue;
        }
        if let Some(paren_pos) = trimmed.find('(') {
            let inner = trimmed[paren_pos + 1..].trim_end_matches(')');
            let verts: Vec<usize> = inner
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if !verts.is_empty() {
                faces.push(verts);
            }
        }
    }
    faces
}
/// Parse an OpenFOAM boundary file into a list of `FoamPatch`.
///
/// Handles both the top-level list wrapper `N ( ... )` and bare dictionary entries.
pub fn parse_foam_boundary(input: &str) -> Vec<FoamPatch> {
    let cleaned = strip_foam_comments(input);
    let tokens = tokenise_foam(&cleaned);
    let mut pos = 0;
    while pos < tokens.len() && tokens[pos] != "(" {
        pos += 1;
    }
    if pos < tokens.len() {
        pos += 1;
    }
    let mut patches = Vec::new();
    while pos < tokens.len() {
        if tokens[pos] == ")" {
            break;
        }
        if tokens[pos] == ";" {
            pos += 1;
            continue;
        }
        let name = tokens[pos].clone();
        pos += 1;
        if pos < tokens.len() && tokens[pos] == "{" {
            pos += 1;
            let (sub, next) = parse_dict_tokens(&tokens, pos);
            pos = next;
            let patch_type = sub.get_word("type").unwrap_or("patch").to_string();
            let n_faces = sub.get_scalar("nFaces").unwrap_or(0.0) as usize;
            let start_face = sub.get_scalar("startFace").unwrap_or(0.0) as usize;
            patches.push(FoamPatch {
                name,
                patch_type,
                start_face,
                n_faces,
            });
        }
    }
    patches
}
/// Parse a uniform or non-uniform scalar field from OpenFOAM field file text.
pub fn parse_foam_scalar_field(input: &str) -> FieldValues {
    let cleaned = strip_foam_comments(input);
    if let Some(pos) = cleaned.find("internalField") {
        let rest = &cleaned[pos..];
        if let Some(u_pos) = rest.find("uniform") {
            let after_uniform = &rest[u_pos + 7..];
            let after_uniform = after_uniform.trim();
            if after_uniform.starts_with('(')
                && let Some(end) = after_uniform.find(')')
            {
                let inner = &after_uniform[1..end];
                let parts: Vec<f64> = inner
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if parts.len() == 3 {
                    return FieldValues::UniformVec([parts[0], parts[1], parts[2]]);
                }
            }
            let val_str: String = after_uniform
                .chars()
                .take_while(|c| !c.is_whitespace() && *c != ';')
                .collect();
            if let Ok(v) = val_str.parse::<f64>() {
                return FieldValues::Uniform(v);
            }
        }
        if rest.contains("nonuniform") {
            if rest.contains("List<vector>") {
                return parse_nonuniform_vector_list(&cleaned);
            }
            return parse_nonuniform_scalar_list(&cleaned);
        }
    }
    FieldValues::Uniform(0.0)
}
pub(super) fn parse_nonuniform_scalar_list(input: &str) -> FieldValues {
    let mut vals = Vec::new();
    let mut in_list = false;
    let mut found_marker = false;
    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.contains("nonuniform") {
            found_marker = true;
            continue;
        }
        if found_marker && !in_list {
            if trimmed == "(" {
                in_list = true;
                continue;
            }
            if trimmed.parse::<usize>().is_ok() {
                continue;
            }
        }
        if in_list {
            if trimmed == ")" || trimmed == ");" {
                break;
            }
            if let Ok(v) = trimmed.parse::<f64>() {
                vals.push(v);
            }
        }
    }
    FieldValues::NonUniform(vals)
}
pub(super) fn parse_nonuniform_vector_list(input: &str) -> FieldValues {
    let mut vals = Vec::new();
    let mut in_list = false;
    let mut found_marker = false;
    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.contains("nonuniform") {
            found_marker = true;
            continue;
        }
        if found_marker && !in_list {
            if trimmed == "(" {
                in_list = true;
                continue;
            }
            if trimmed.parse::<usize>().is_ok() {
                continue;
            }
        }
        if in_list {
            if trimmed == ")" || trimmed == ");" {
                break;
            }
            let inner = trimmed.trim_start_matches('(').trim_end_matches(')');
            let parts: Vec<f64> = inner
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if parts.len() == 3 {
                vals.push([parts[0], parts[1], parts[2]]);
            }
        }
    }
    FieldValues::NonUniformVec(vals)
}
/// Sort time directory names numerically.
pub fn sort_time_dirs(dirs: &[String]) -> Vec<String> {
    let mut numeric: Vec<(f64, String)> = dirs
        .iter()
        .filter_map(|d| d.parse::<f64>().ok().map(|v| (v, d.clone())))
        .collect();
    numeric.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    numeric.into_iter().map(|(_, s)| s).collect()
}
/// Check if a directory name is a valid OpenFOAM time directory.
pub fn is_time_dir(name: &str) -> bool {
    name.parse::<f64>().is_ok()
}
/// Find the latest time from a list of directory names.
pub fn latest_time(dirs: &[String]) -> Option<f64> {
    dirs.iter()
        .filter_map(|d| d.parse::<f64>().ok())
        .fold(None, |acc, v| match acc {
            None => Some(v),
            Some(prev) => Some(if v > prev { v } else { prev }),
        })
}
/// Generate the path string for a field file at a given time.
pub fn field_path(case_dir: &str, time: f64, field_name: &str) -> String {
    if (time - time.round()).abs() < 1e-12 {
        format!("{}/{}/{}", case_dir, time as i64, field_name)
    } else {
        format!("{}/{}/{}", case_dir, time, field_name)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::openfoam::ControlDict;

    use crate::openfoam::FoamMesh;
    use crate::openfoam::FvSchemes;
    use crate::openfoam::FvSolution;
    use crate::openfoam::TransportProperties;
    use crate::openfoam::types::*;
    #[test]
    fn test_foam_header() {
        let h = foam_header("vectorField", "points");
        assert!(h.contains("FoamFile"), "header must contain FoamFile");
        assert!(h.contains("class"), "header must contain 'class'");
        assert!(h.contains("vectorField"));
        assert!(h.contains("points"));
    }
    #[test]
    fn test_box_mesh_point_count() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.points.len(), 3 * 3 * 3, "2x2x2 mesh → 27 points");
    }
    #[test]
    fn test_box_mesh_cell_count() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.n_cells, 8, "2x2x2 mesh → 8 cells");
    }
    #[test]
    fn test_write_points_format() {
        let mesh = FoamMesh::box_mesh(2.0, 3.0, 4.0, 1, 1, 1);
        let out = mesh.write_points();
        assert!(
            out.contains("(0 0 0)") || out.contains("(0.0 0.0 0.0)") || out.contains("(0 0 0)")
        );
        assert!(out.contains("2"), "lx=2.0 should appear");
        assert!(out.contains("3"), "ly=3.0 should appear");
        assert!(out.contains("4"), "lz=4.0 should appear");
        assert!(out.contains("8\n("));
    }
    #[test]
    fn test_write_faces_format() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let out = mesh.write_faces();
        assert!(out.contains("6\n("), "should list 6 faces");
        assert!(out.contains("4("), "faces should be quads");
    }
    #[test]
    fn test_control_dict_format() {
        let cd = ControlDict::new("icoFoam", 1.0, 0.01);
        let out = cd.to_string();
        assert!(out.contains("application"), "must contain 'application'");
        assert!(out.contains("endTime"), "must contain 'endTime'");
        assert!(out.contains("icoFoam"));
        assert!(out.contains("1"));
        assert!(out.contains("0.01"));
    }
    #[test]
    fn test_foam_scalar_field_uniform() {
        let field = FoamField {
            n_cells: 8,
            field_name: "p".to_string(),
            field_class: "volScalarField".to_string(),
            dimensions: "[0 2 -2 0 0 0 0]".to_string(),
            internal_values: FieldValues::Uniform(0.0),
            boundary_conditions: vec![FoamBc {
                patch_name: "inlet".to_string(),
                bc_type: "zeroGradient".to_string(),
                value: None,
            }],
        };
        let out = field.to_string();
        assert!(out.contains("uniform 0"), "uniform scalar zero");
        assert!(out.contains("volScalarField"));
        assert!(out.contains("zeroGradient"));
    }
    #[test]
    fn test_foam_vector_field() {
        let field = FoamField {
            n_cells: 4,
            field_name: "U".to_string(),
            field_class: "volVectorField".to_string(),
            dimensions: "[0 1 -1 0 0 0 0]".to_string(),
            internal_values: FieldValues::UniformVec([1.0, 0.0, 0.0]),
            boundary_conditions: vec![
                FoamBc {
                    patch_name: "walls".to_string(),
                    bc_type: "noSlip".to_string(),
                    value: None,
                },
                FoamBc {
                    patch_name: "inlet".to_string(),
                    bc_type: "fixedValue".to_string(),
                    value: Some("uniform (1 0 0)".to_string()),
                },
            ],
        };
        let out = field.to_string();
        assert!(out.contains("uniform (1 0 0)"), "uniform vec (1 0 0)");
        assert!(out.contains("volVectorField"));
        assert!(out.contains("noSlip"));
        assert!(out.contains("fixedValue"));
        assert!(out.contains("(1 0 0)"));
    }
    #[test]
    fn test_foam_dict_empty() {
        let d = FoamDict::new();
        assert!(d.is_empty());
        assert_eq!(d.len(), 0);
    }
    #[test]
    fn test_foam_dict_insert_and_get() {
        let mut d = FoamDict::new();
        d.insert("alpha", FoamValue::Scalar(0.5));
        d.insert("solver", FoamValue::Word("PCG".to_string()));
        assert_eq!(d.len(), 2);
        assert_eq!(d.get_scalar("alpha"), Some(0.5));
        assert_eq!(d.get_word("solver"), Some("PCG"));
        assert!(d.get_scalar("missing").is_none());
    }
    #[test]
    fn test_foam_dict_keys() {
        let mut d = FoamDict::new();
        d.insert("a", FoamValue::Scalar(1.0));
        d.insert("b", FoamValue::Scalar(2.0));
        assert_eq!(d.keys(), vec!["a", "b"]);
    }
    #[test]
    fn test_foam_dict_vector() {
        let mut d = FoamDict::new();
        d.insert("gravity", FoamValue::Vector([0.0, -9.81, 0.0]));
        let v = d.get_vector("gravity").unwrap();
        assert!((v[1] + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_foam_dict_sub_dict() {
        let mut sub = FoamDict::new();
        sub.insert("tolerance", FoamValue::Scalar(1e-6));
        let mut d = FoamDict::new();
        d.insert("solvers", FoamValue::Dict(sub));
        let s = d.get_dict("solvers").unwrap();
        assert_eq!(s.get_scalar("tolerance"), Some(1e-6));
    }
    #[test]
    fn test_foam_dict_parse_simple() {
        let input = r#"
            application icoFoam;
            endTime 1.0;
            deltaT 0.01;
        "#;
        let d = FoamDict::parse(input);
        assert_eq!(d.get_word("application"), Some("icoFoam"));
        assert_eq!(d.get_scalar("endTime"), Some(1.0));
        assert_eq!(d.get_scalar("deltaT"), Some(0.01));
    }
    #[test]
    fn test_foam_dict_parse_with_comments() {
        let input = r#"
            // This is a comment
            value 42; /* inline comment */
            name test;
        "#;
        let d = FoamDict::parse(input);
        assert_eq!(d.get_scalar("value"), Some(42.0));
        assert_eq!(d.get_word("name"), Some("test"));
    }
    #[test]
    fn test_foam_dict_parse_nested() {
        let input = r#"
            outer
            {
                inner
                {
                    x 10;
                }
            }
        "#;
        let d = FoamDict::parse(input);
        let outer = d.get_dict("outer").unwrap();
        let inner = outer.get_dict("inner").unwrap();
        assert_eq!(inner.get_scalar("x"), Some(10.0));
    }
    #[test]
    fn test_foam_dict_parse_vector() {
        let input = "gravity (0 -9.81 0);";
        let d = FoamDict::parse(input);
        let v = d.get_vector("gravity").unwrap();
        assert!((v[1] + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_foam_dict_to_foam_string() {
        let mut d = FoamDict::new();
        d.insert("solver", FoamValue::Word("PCG".to_string()));
        d.insert("tolerance", FoamValue::Scalar(1e-6));
        let out = d.to_foam_string(0);
        assert!(out.contains("solver"));
        assert!(out.contains("PCG"));
        assert!(out.contains("tolerance"));
    }
    #[test]
    fn test_strip_foam_comments() {
        let input = "hello // world\nfoo /* bar */ baz";
        let out = strip_foam_comments(input);
        assert!(out.contains("hello"));
        assert!(!out.contains("world"));
        assert!(out.contains("foo"));
        assert!(!out.contains("bar"));
        assert!(out.contains("baz"));
    }
    #[test]
    fn test_parse_foam_points() {
        let input = r#"
FoamFile { version 2.0; class vectorField; object points; }
8
(
(0 0 0)
(1 0 0)
(1 1 0)
(0 1 0)
(0 0 1)
(1 0 1)
(1 1 1)
(0 1 1)
)
"#;
        let pts = parse_foam_points(input);
        assert_eq!(pts.len(), 8);
        assert!((pts[0][0]).abs() < 1e-10);
        assert!((pts[1][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_parse_foam_label_list() {
        let input = r#"
FoamFile { version 2.0; class labelList; object owner; }
4
(
0
1
2
3
)
"#;
        let labels = parse_foam_label_list(input);
        assert_eq!(labels, vec![0, 1, 2, 3]);
    }
    #[test]
    fn test_parse_foam_faces() {
        let input = r#"
FoamFile { version 2.0; class faceList; object faces; }
2
(
4(0 1 2 3)
3(4 5 6)
)
"#;
        let faces = parse_foam_faces(input);
        assert_eq!(faces.len(), 2);
        assert_eq!(faces[0], vec![0, 1, 2, 3]);
        assert_eq!(faces[1], vec![4, 5, 6]);
    }
    #[test]
    fn test_parse_foam_boundary() {
        let input = r#"
2
(
    inlet
    {
        type patch;
        nFaces 10;
        startFace 100;
    }
    walls
    {
        type wall;
        nFaces 20;
        startFace 110;
    }
)
"#;
        let patches = parse_foam_boundary(input);
        assert_eq!(patches.len(), 2);
        assert_eq!(patches[0].name, "inlet");
        assert_eq!(patches[0].patch_type, "patch");
        assert_eq!(patches[0].n_faces, 10);
        assert_eq!(patches[0].start_face, 100);
        assert_eq!(patches[1].name, "walls");
        assert_eq!(patches[1].patch_type, "wall");
    }
    #[test]
    fn test_parse_uniform_scalar() {
        let input = r#"
FoamFile { version 2.0; class volScalarField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField { }
"#;
        match parse_foam_scalar_field(input) {
            FieldValues::Uniform(v) => assert!((v).abs() < 1e-10),
            _ => panic!("expected Uniform"),
        }
    }
    #[test]
    fn test_parse_uniform_vector() {
        let input = r#"
FoamFile { version 2.0; class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField   uniform (1 0 0);
boundaryField { }
"#;
        match parse_foam_scalar_field(input) {
            FieldValues::UniformVec(v) => {
                assert!((v[0] - 1.0).abs() < 1e-10);
                assert!((v[1]).abs() < 1e-10);
            }
            _ => panic!("expected UniformVec"),
        }
    }
    #[test]
    fn test_parse_nonuniform_scalar() {
        let input = r#"
internalField   nonuniform List<scalar>
3
(
1.5
2.5
3.5
);
"#;
        match parse_foam_scalar_field(input) {
            FieldValues::NonUniform(v) => {
                assert_eq!(v.len(), 3);
                assert!((v[0] - 1.5).abs() < 1e-10);
            }
            _ => panic!("expected NonUniform"),
        }
    }
    #[test]
    fn test_parse_nonuniform_vector() {
        let input = r#"
internalField   nonuniform List<vector>
2
(
(1 0 0)
(0 1 0)
);
"#;
        match parse_foam_scalar_field(input) {
            FieldValues::NonUniformVec(v) => {
                assert_eq!(v.len(), 2);
                assert!((v[0][0] - 1.0).abs() < 1e-10);
                assert!((v[1][1] - 1.0).abs() < 1e-10);
            }
            _ => panic!("expected NonUniformVec"),
        }
    }
    #[test]
    fn test_bc_zero_gradient() {
        let bc = FoamBc::zero_gradient("outlet");
        assert_eq!(bc.bc_type, "zeroGradient");
        assert!(bc.value.is_none());
    }
    #[test]
    fn test_bc_fixed_scalar() {
        let bc = FoamBc::fixed_scalar("inlet", 1.0);
        assert_eq!(bc.bc_type, "fixedValue");
        assert!(bc.value.as_ref().unwrap().contains("uniform 1"));
    }
    #[test]
    fn test_bc_fixed_vector() {
        let bc = FoamBc::fixed_vector("inlet", [1.0, 0.0, 0.0]);
        assert!(bc.value.as_ref().unwrap().contains("(1 0 0)"));
    }
    #[test]
    fn test_bc_no_slip() {
        let bc = FoamBc::no_slip("walls");
        assert_eq!(bc.bc_type, "noSlip");
    }
    #[test]
    fn test_bc_symmetry() {
        let bc = FoamBc::symmetry("sym");
        assert_eq!(bc.bc_type, "symmetry");
    }
    #[test]
    fn test_bc_empty() {
        let bc = FoamBc::empty("frontAndBack");
        assert_eq!(bc.bc_type, "empty");
    }
    #[test]
    fn test_is_time_dir() {
        assert!(is_time_dir("0"));
        assert!(is_time_dir("0.5"));
        assert!(is_time_dir("100"));
        assert!(!is_time_dir("constant"));
        assert!(!is_time_dir("system"));
    }
    #[test]
    fn test_sort_time_dirs() {
        let dirs: Vec<String> = vec!["1", "0.5", "0", "10", "2"]
            .into_iter()
            .map(String::from)
            .collect();
        let sorted = sort_time_dirs(&dirs);
        assert_eq!(sorted, vec!["0", "0.5", "1", "2", "10"]);
    }
    #[test]
    fn test_latest_time() {
        let dirs: Vec<String> = vec!["0", "0.5", "1", "2"]
            .into_iter()
            .map(String::from)
            .collect();
        assert_eq!(latest_time(&dirs), Some(2.0));
    }
    #[test]
    fn test_latest_time_empty() {
        let dirs: Vec<String> = Vec::new();
        assert_eq!(latest_time(&dirs), None);
    }
    #[test]
    fn test_field_path_integer_time() {
        let p = field_path("/case", 0.0, "U");
        assert_eq!(p, "/case/0/U");
    }
    #[test]
    fn test_field_path_fractional_time() {
        let p = field_path("/case", 0.5, "p");
        assert_eq!(p, "/case/0.5/p");
    }
    #[test]
    fn test_fv_schemes_default() {
        let s = FvSchemes::default_second_order();
        let out = s.to_string();
        assert!(out.contains("ddtSchemes"));
        assert!(out.contains("Euler"));
        assert!(out.contains("gradSchemes"));
        assert!(out.contains("Gauss linear"));
        assert!(out.contains("divSchemes"));
        assert!(out.contains("laplacianSchemes"));
        assert!(out.contains("snGradSchemes"));
    }
    #[test]
    fn test_fv_solution_default_piso() {
        let s = FvSolution::default_piso();
        let out = s.to_string();
        assert!(out.contains("solvers"));
        assert!(out.contains("PISO"));
        assert!(out.contains("nCorrectors"));
        assert!(out.contains("PCG"));
        assert!(out.contains("smoothSolver"));
    }
    #[test]
    fn test_transport_properties_newtonian() {
        let tp = TransportProperties::newtonian(1e-6);
        let out = tp.to_string();
        assert!(out.contains("Newtonian"));
        assert!(out.contains("nu"));
        assert!(out.contains("[0 2 -1 0 0 0 0]"));
    }
    #[test]
    fn test_mesh_internal_faces() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.n_internal_faces(), 12);
    }
    #[test]
    fn test_mesh_boundary_faces() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.n_boundary_faces(), 24);
    }
    #[test]
    fn test_mesh_total_faces() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.n_faces(), 36);
    }
    #[test]
    fn test_mesh_bounding_box() {
        let mesh = FoamMesh::box_mesh(2.0, 3.0, 4.0, 1, 1, 1);
        let (min, max) = mesh.bounding_box();
        assert!((min[0]).abs() < 1e-10);
        assert!((max[0] - 2.0).abs() < 1e-10);
        assert!((max[1] - 3.0).abs() < 1e-10);
        assert!((max[2] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_centre() {
        let mesh = FoamMesh::box_mesh(2.0, 4.0, 6.0, 1, 1, 1);
        let c = mesh.centre();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
        assert!((c[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_check_topology() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        assert!(mesh.check_topology());
    }
    #[test]
    fn test_mesh_find_patch() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        assert!(mesh.find_patch("xmin").is_some());
        assert!(mesh.find_patch("nonexistent").is_none());
    }
    #[test]
    fn test_mesh_patch_names() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let names = mesh.patch_names();
        assert_eq!(names.len(), 6);
        assert!(names.contains(&"xmin"));
        assert!(names.contains(&"zmax"));
    }
    #[test]
    fn test_write_then_parse_points_roundtrip() {
        let mesh = FoamMesh::box_mesh(1.0, 2.0, 3.0, 1, 1, 1);
        let written = mesh.write_points();
        let parsed = parse_foam_points(&written);
        assert_eq!(parsed.len(), mesh.points.len());
        for (orig, parsed_pt) in mesh.points.iter().zip(parsed.iter()) {
            for i in 0..3 {
                assert!((orig[i] - parsed_pt[i]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_write_then_parse_faces_roundtrip() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let written = mesh.write_faces();
        let parsed = parse_foam_faces(&written);
        assert_eq!(parsed.len(), mesh.faces.len());
    }
    #[test]
    fn test_write_then_parse_owner_roundtrip() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let written = mesh.write_owner();
        let parsed = parse_foam_label_list(&written);
        assert_eq!(parsed.len(), mesh.owner.len());
        for (orig, parsed_v) in mesh.owner.iter().zip(parsed.iter()) {
            assert_eq!(*orig as i64, *parsed_v);
        }
    }
    #[test]
    fn test_nonuniform_scalar_field_write_and_parse() {
        let field = FoamField {
            n_cells: 3,
            field_name: "T".to_string(),
            field_class: "volScalarField".to_string(),
            dimensions: "[0 0 0 1 0 0 0]".to_string(),
            internal_values: FieldValues::NonUniform(vec![300.0, 310.0, 320.0]),
            boundary_conditions: vec![],
        };
        let written = field.to_string();
        match parse_foam_scalar_field(&written) {
            FieldValues::NonUniform(v) => {
                assert_eq!(v.len(), 3);
                assert!((v[0] - 300.0).abs() < 1e-10);
                assert!((v[2] - 320.0).abs() < 1e-10);
            }
            _ => panic!("expected NonUniform"),
        }
    }
    #[test]
    fn test_nonuniform_vector_field_write_and_parse() {
        let field = FoamField {
            n_cells: 2,
            field_name: "U".to_string(),
            field_class: "volVectorField".to_string(),
            dimensions: "[0 1 -1 0 0 0 0]".to_string(),
            internal_values: FieldValues::NonUniformVec(vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            boundary_conditions: vec![],
        };
        let written = field.to_string();
        match parse_foam_scalar_field(&written) {
            FieldValues::NonUniformVec(v) => {
                assert_eq!(v.len(), 2);
                assert!((v[0][0] - 1.0).abs() < 1e-10);
                assert!((v[1][1] - 1.0).abs() < 1e-10);
            }
            _ => panic!("expected NonUniformVec"),
        }
    }
}
/// Discover time directories from a list of directory entry names.
///
/// Filters names that parse as floating-point numbers and returns
/// them sorted in ascending order along with their numeric values.
pub fn discover_time_dirs(all_entries: &[String]) -> Vec<FoamTimeDir> {
    let mut result: Vec<FoamTimeDir> = all_entries
        .iter()
        .filter_map(|name| {
            name.parse::<f64>().ok().map(|t| FoamTimeDir {
                time: t,
                dir_name: name.clone(),
                fields: Vec::new(),
            })
        })
        .collect();
    result.sort_by(|a, b| {
        a.time
            .partial_cmp(&b.time)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    result
}
/// Discover time directories and associate field names with each time step.
///
/// `field_map` is a list of `(dir_name, field_name)` pairs that describe
/// which fields exist in which time directory.
pub fn discover_time_dirs_with_fields(
    all_entries: &[String],
    field_map: &[(String, String)],
) -> Vec<FoamTimeDir> {
    let mut dirs = discover_time_dirs(all_entries);
    for dir in &mut dirs {
        for (d, f) in field_map {
            if *d == dir.dir_name {
                dir.fields.push(f.clone());
            }
        }
        dir.fields.sort();
    }
    dirs
}
/// Write a non-uniform scalar field to an OpenFOAM ASCII field file string.
///
/// This is a convenience wrapper that creates a `FoamField` with `NonUniform`
/// internal values and no boundary conditions, then renders it.
pub fn write_scalar_field_ascii(
    field_name: &str,
    dimensions: &str,
    values: &[f64],
    bcs: Vec<FoamBc>,
) -> String {
    let field = FoamField {
        n_cells: values.len(),
        field_name: field_name.to_string(),
        field_class: "volScalarField".to_string(),
        dimensions: dimensions.to_string(),
        internal_values: FieldValues::NonUniform(values.to_vec()),
        boundary_conditions: bcs,
    };
    field.to_string()
}
/// Write a uniform scalar field file.
pub fn write_uniform_scalar_field(
    field_name: &str,
    dimensions: &str,
    value: f64,
    bcs: Vec<FoamBc>,
) -> String {
    let field = FoamField {
        n_cells: 0,
        field_name: field_name.to_string(),
        field_class: "volScalarField".to_string(),
        dimensions: dimensions.to_string(),
        internal_values: FieldValues::Uniform(value),
        boundary_conditions: bcs,
    };
    field.to_string()
}
/// Write a uniform vector field file.
pub fn write_uniform_vector_field(
    field_name: &str,
    dimensions: &str,
    value: [f64; 3],
    bcs: Vec<FoamBc>,
) -> String {
    let field = FoamField {
        n_cells: 0,
        field_name: field_name.to_string(),
        field_class: "volVectorField".to_string(),
        dimensions: dimensions.to_string(),
        internal_values: FieldValues::UniformVec(value),
        boundary_conditions: bcs,
    };
    field.to_string()
}
/// Parse OpenFOAM solver residuals from a log file string.
///
/// Recognises lines of the form:
/// ```text
/// Time = 0.1
/// smoothSolver:  Solving for Ux, Initial residual = 1e-3, Final residual = 1e-6, No Iterations 5
/// ```
pub fn parse_foam_residuals(log: &str) -> Vec<FoamResidual> {
    let mut residuals = Vec::new();
    let mut current_time = 0.0_f64;
    for line in log.lines() {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("Time =") {
            if let Ok(t) = rest.trim().parse::<f64>() {
                current_time = t;
            }
            continue;
        }
        if trimmed.contains("Solving for") {
            let field = extract_foam_field_name(trimmed).unwrap_or_default();
            let init = extract_foam_value(trimmed, "Initial residual =").unwrap_or(0.0);
            let fin = extract_foam_value(trimmed, "Final residual =").unwrap_or(0.0);
            let n_iter = extract_foam_int(trimmed, "No Iterations").unwrap_or(0);
            residuals.push(FoamResidual {
                time: current_time,
                field,
                initial_residual: init,
                final_residual: fin,
                n_iterations: n_iter,
            });
        }
    }
    residuals
}
/// Extract the field name from a "Solving for X" line.
pub(super) fn extract_foam_field_name(line: &str) -> Option<String> {
    let pos = line.find("Solving for")?;
    let rest = &line[pos + "Solving for".len()..];
    let trimmed = rest.trim_start();
    let end = trimmed.find([',', ' ', '\t']).unwrap_or(trimmed.len());
    Some(trimmed[..end].to_string())
}
/// Extract a float value after a `key` in a log line.
pub(super) fn extract_foam_value(line: &str, key: &str) -> Option<f64> {
    let pos = line.find(key)?;
    let rest = &line[pos + key.len()..].trim_start();
    let end = rest
        .find(|c: char| c == ',' || c.is_whitespace())
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}
/// Extract an integer value after a `key` in a log line.
pub(super) fn extract_foam_int(line: &str, key: &str) -> Option<usize> {
    let pos = line.find(key)?;
    let rest = &line[pos + key.len()..].trim_start();
    let end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}
/// Encode a GROMACS-style dimension set for common physical quantities.
///
/// Returns an OpenFOAM dimension string `[kg m s K mol A cd]`.
pub fn foam_dimensions(kg: i32, m: i32, s: i32, k: i32, mol: i32, a: i32, cd: i32) -> String {
    format!("[{} {} {} {} {} {} {}]", kg, m, s, k, mol, a, cd)
}
/// Dimension set for pressure (Pa = kg m⁻¹ s⁻²).
pub fn dims_pressure() -> String {
    foam_dimensions(1, -1, -2, 0, 0, 0, 0)
}
/// Dimension set for velocity (m s⁻¹).
pub fn dims_velocity() -> String {
    foam_dimensions(0, 1, -1, 0, 0, 0, 0)
}
/// Dimension set for temperature (K).
pub fn dims_temperature() -> String {
    foam_dimensions(0, 0, 0, 1, 0, 0, 0)
}
/// Dimension set for kinematic viscosity (m² s⁻¹).
pub fn dims_kinematic_viscosity() -> String {
    foam_dimensions(0, 2, -1, 0, 0, 0, 0)
}
#[cfg(test)]
mod tests_openfoam_ext {
    use super::*;

    use crate::openfoam::FoamFileHeader;

    #[test]
    fn test_foam_file_header_parse_basic() {
        let input = r#"
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
"#;
        let hdr = FoamFileHeader::parse(input).unwrap();
        assert!((hdr.version - 2.0).abs() < 1e-10);
        assert_eq!(hdr.format, "ascii");
        assert_eq!(hdr.class, "volScalarField");
        assert_eq!(hdr.object, "p");
        assert!(hdr.note.is_none());
        assert!(hdr.location.is_none());
    }
    #[test]
    fn test_foam_file_header_with_location_and_note() {
        let input = r#"
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
    note        "velocity field";
}
"#;
        let hdr = FoamFileHeader::parse(input).unwrap();
        assert_eq!(hdr.class, "volVectorField");
        assert_eq!(hdr.object, "U");
        assert_eq!(hdr.location.as_deref(), Some("0"));
        assert_eq!(hdr.note.as_deref(), Some("velocity field"));
    }
    #[test]
    fn test_foam_file_header_roundtrip() {
        let input = r#"
FoamFile
{
    version     2.0;
    format      ascii;
    class       labelList;
    object      owner;
}
"#;
        let hdr = FoamFileHeader::parse(input).unwrap();
        let written = hdr.to_string();
        let hdr2 = FoamFileHeader::parse(&written).unwrap();
        assert_eq!(hdr2.class, hdr.class);
        assert_eq!(hdr2.object, hdr.object);
    }
    #[test]
    fn test_foam_file_header_missing() {
        let input = "just some text without a FoamFile block";
        assert!(FoamFileHeader::parse(input).is_none());
    }
    #[test]
    fn test_foam_file_header_write_format() {
        let hdr = FoamFileHeader {
            version: 2.0,
            format: "ascii".to_string(),
            class: "polyMesh".to_string(),
            object: "points".to_string(),
            note: None,
            location: Some("constant/polyMesh".to_string()),
        };
        let s = hdr.to_string();
        assert!(s.contains("FoamFile"));
        assert!(s.contains("polyMesh"));
        assert!(s.contains("points"));
        assert!(s.contains("constant/polyMesh"));
    }
    #[test]
    fn test_discover_time_dirs_basic() {
        let entries: Vec<String> = vec!["0", "0.5", "1", "constant", "system"]
            .into_iter()
            .map(String::from)
            .collect();
        let dirs = discover_time_dirs(&entries);
        assert_eq!(dirs.len(), 3);
        assert!((dirs[0].time).abs() < 1e-10);
        assert!((dirs[1].time - 0.5).abs() < 1e-10);
        assert!((dirs[2].time - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_discover_time_dirs_sorted() {
        let entries: Vec<String> = vec!["10", "2", "0.1", "1"]
            .into_iter()
            .map(String::from)
            .collect();
        let dirs = discover_time_dirs(&entries);
        assert!((dirs[0].time - 0.1).abs() < 1e-10);
        assert!((dirs[3].time - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_discover_time_dirs_with_fields() {
        let entries: Vec<String> = vec!["0", "1"].into_iter().map(String::from).collect();
        let field_map: Vec<(String, String)> = vec![
            ("0".to_string(), "p".to_string()),
            ("0".to_string(), "U".to_string()),
            ("1".to_string(), "p".to_string()),
        ];
        let dirs = discover_time_dirs_with_fields(&entries, &field_map);
        assert_eq!(dirs[0].fields.len(), 2);
        assert_eq!(dirs[1].fields.len(), 1);
    }
    #[test]
    fn test_discover_time_dirs_empty() {
        let entries: Vec<String> = vec!["constant", "system"]
            .into_iter()
            .map(String::from)
            .collect();
        let dirs = discover_time_dirs(&entries);
        assert!(dirs.is_empty());
    }
    #[test]
    fn test_write_scalar_field_ascii() {
        let vals = vec![1.0, 2.0, 3.0];
        let out = write_scalar_field_ascii("T", "[0 0 0 1 0 0 0]", &vals, vec![]);
        assert!(out.contains("volScalarField"));
        assert!(out.contains("T"));
        assert!(out.contains("[0 0 0 1 0 0 0]"));
        assert!(out.contains("nonuniform"));
    }
    #[test]
    fn test_write_uniform_scalar_field() {
        let out = write_uniform_scalar_field(
            "p",
            "[0 2 -2 0 0 0 0]",
            0.0,
            vec![FoamBc::zero_gradient("outlet")],
        );
        assert!(out.contains("uniform 0"));
        assert!(out.contains("zeroGradient"));
    }
    #[test]
    fn test_write_uniform_vector_field() {
        let out = write_uniform_vector_field(
            "U",
            "[0 1 -1 0 0 0 0]",
            [1.0, 0.0, 0.0],
            vec![FoamBc::no_slip("walls")],
        );
        assert!(out.contains("volVectorField"));
        assert!(out.contains("(1 0 0)"));
        assert!(out.contains("noSlip"));
    }
    #[test]
    fn test_parse_foam_residuals_basic() {
        let log = "Time = 0.1\n\
smoothSolver:  Solving for Ux, Initial residual = 1e-3, Final residual = 1e-6, No Iterations 5\n\
PCG:  Solving for p, Initial residual = 1e-4, Final residual = 1e-7, No Iterations 10\n\
Time = 0.2\n\
smoothSolver:  Solving for Ux, Initial residual = 9e-4, Final residual = 8e-7, No Iterations 4\n";
        let res = parse_foam_residuals(log);
        assert_eq!(res.len(), 3);
        assert!((res[0].time - 0.1).abs() < 1e-10);
        assert_eq!(res[0].field, "Ux");
        assert!((res[0].initial_residual - 1e-3).abs() < 1e-12);
        assert_eq!(res[0].n_iterations, 5);
        assert_eq!(res[1].field, "p");
        assert!((res[2].time - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_parse_foam_residuals_empty_log() {
        let res = parse_foam_residuals("");
        assert!(res.is_empty());
    }
    #[test]
    fn test_parse_foam_residuals_no_time_line() {
        let log = "smoothSolver:  Solving for Ux, Initial residual = 1e-3, Final residual = 1e-6, No Iterations 5\n";
        let res = parse_foam_residuals(log);
        assert_eq!(res.len(), 1);
        assert!((res[0].time).abs() < 1e-10);
    }
    #[test]
    fn test_foam_dimensions_pressure() {
        let d = dims_pressure();
        assert_eq!(d, "[1 -1 -2 0 0 0 0]");
    }
    #[test]
    fn test_foam_dimensions_velocity() {
        let d = dims_velocity();
        assert_eq!(d, "[0 1 -1 0 0 0 0]");
    }
    #[test]
    fn test_foam_dimensions_temperature() {
        let d = dims_temperature();
        assert_eq!(d, "[0 0 0 1 0 0 0]");
    }
    #[test]
    fn test_foam_dimensions_custom() {
        let d = foam_dimensions(1, 2, -3, 0, 0, 0, 0);
        assert_eq!(d, "[1 2 -3 0 0 0 0]");
    }
    #[test]
    fn test_foam_header_class_matches_field_write() {
        let out = write_scalar_field_ascii("T", "[0 0 0 1 0 0 0]", &[300.0, 310.0], vec![]);
        let hdr = FoamFileHeader::parse(&out).unwrap();
        assert_eq!(hdr.class, "volScalarField");
        assert_eq!(hdr.object, "T");
    }
    #[test]
    fn test_write_nonuniform_and_parse_back() {
        let vals = vec![100.0, 200.0, 300.0, 400.0, 500.0];
        let out = write_scalar_field_ascii("p", "[1 -1 -2 0 0 0 0]", &vals, vec![]);
        match parse_foam_scalar_field(&out) {
            FieldValues::NonUniform(v) => {
                assert_eq!(v.len(), 5);
                assert!((v[0] - 100.0).abs() < 1e-10);
                assert!((v[4] - 500.0).abs() < 1e-10);
            }
            _ => panic!("expected NonUniform"),
        }
    }
    #[test]
    fn test_foam_kinematic_viscosity_dim() {
        let d = dims_kinematic_viscosity();
        assert_eq!(d, "[0 2 -1 0 0 0 0]");
    }
}
/// Compute the face area of a planar polygon given its vertex coordinates.
///
/// Uses the cross-product shoelace formula projected onto the best-fit plane.
/// For convex quads and triangles this gives an exact result.
pub fn face_area(vertices: &[[f64; 3]]) -> f64 {
    if vertices.len() < 3 {
        return 0.0;
    }
    let mut area_vec = [0.0_f64; 3];
    let v0 = vertices[0];
    for i in 1..vertices.len() - 1 {
        let v1 = vertices[i];
        let v2 = vertices[i + 1];
        let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        area_vec[0] += e1[1] * e2[2] - e1[2] * e2[1];
        area_vec[1] += e1[2] * e2[0] - e1[0] * e2[2];
        area_vec[2] += e1[0] * e2[1] - e1[1] * e2[0];
    }
    0.5 * (area_vec[0] * area_vec[0] + area_vec[1] * area_vec[1] + area_vec[2] * area_vec[2]).sqrt()
}
/// Compute the face centroid (arithmetic mean of vertices).
pub fn face_centroid(vertices: &[[f64; 3]]) -> [f64; 3] {
    if vertices.is_empty() {
        return [0.0; 3];
    }
    let n = vertices.len() as f64;
    let mut c = [0.0_f64; 3];
    for v in vertices {
        c[0] += v[0];
        c[1] += v[1];
        c[2] += v[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}
/// Compute the face normal (unit vector) using the cross-product approach.
///
/// Returns a zero vector if the face is degenerate.
pub fn face_normal(vertices: &[[f64; 3]]) -> [f64; 3] {
    if vertices.len() < 3 {
        return [0.0; 3];
    }
    let v0 = vertices[0];
    let v1 = vertices[1];
    let v2 = vertices[2];
    let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if mag < 1e-30 {
        return [0.0; 3];
    }
    [n[0] / mag, n[1] / mag, n[2] / mag]
}
/// Compute the volume of a hexahedral cell given 8 corner vertices in
/// OpenFOAM ordering: bottom-face corners 0-3, top-face corners 4-7.
///
/// Uses the divergence theorem with 6 quad faces decomposed into 2 triangles
/// each. For non-orthogonal cells this is approximate.
pub fn hex_cell_volume(corners: &[[f64; 3]; 8]) -> f64 {
    let tets: [[usize; 4]; 5] = [
        [0, 1, 3, 4],
        [1, 2, 3, 6],
        [4, 5, 6, 1],
        [4, 6, 7, 3],
        [1, 3, 4, 6],
    ];
    let mut vol = 0.0_f64;
    for [a, b, c, d] in &tets {
        vol += tet_signed_volume(corners[*a], corners[*b], corners[*c], corners[*d]).abs();
    }
    vol
}
pub(super) fn tet_signed_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    (cross[0] * ad[0] + cross[1] * ad[1] + cross[2] * ad[2]) / 6.0
}
/// Non-orthogonality angle (degrees) between a face normal and the cell-centre
/// to face-centre vector.
///
/// Returns 0 if the vectors are aligned (ideal orthogonal mesh).
pub fn non_orthogonality(
    face_norm: [f64; 3],
    owner_centre: [f64; 3],
    face_centre: [f64; 3],
) -> f64 {
    let d = [
        face_centre[0] - owner_centre[0],
        face_centre[1] - owner_centre[1],
        face_centre[2] - owner_centre[2],
    ];
    let mag_d = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let mag_n =
        (face_norm[0] * face_norm[0] + face_norm[1] * face_norm[1] + face_norm[2] * face_norm[2])
            .sqrt();
    if mag_d < 1e-30 || mag_n < 1e-30 {
        return 0.0;
    }
    let cos_theta = ((d[0] * face_norm[0] + d[1] * face_norm[1] + d[2] * face_norm[2])
        / (mag_d * mag_n))
        .clamp(-1.0, 1.0);
    cos_theta.acos() * 180.0 / std::f64::consts::PI
}
/// Compute an approximate cell centroid as the mean of all face centroids.
pub fn cell_centroid_from_faces(
    face_indices: &[usize],
    all_faces: &[Vec<usize>],
    all_points: &[[f64; 3]],
) -> [f64; 3] {
    if face_indices.is_empty() {
        return [0.0; 3];
    }
    let mut sum = [0.0_f64; 3];
    let mut count = 0usize;
    for &fi in face_indices {
        if let Some(face) = all_faces.get(fi) {
            let verts: Vec<[f64; 3]> = face
                .iter()
                .filter_map(|&vi| all_points.get(vi).copied())
                .collect();
            if !verts.is_empty() {
                let c = face_centroid(&verts);
                sum[0] += c[0];
                sum[1] += c[1];
                sum[2] += c[2];
                count += 1;
            }
        }
    }
    if count == 0 {
        return [0.0; 3];
    }
    let n = count as f64;
    [sum[0] / n, sum[1] / n, sum[2] / n]
}
/// Compute the cell-centred gradient of a scalar field using a simple
/// Green-Gauss scheme on a structured 3D grid.
///
/// `values[i]` is the scalar at cell `i`, `cell_centres[i]` its position,
/// `faces`, `owner`, `neighbour` are the standard OpenFOAM topology arrays.
///
/// Returns one gradient vector `[dφ/dx, dφ/dy, dφ/dz]` per cell.
pub fn green_gauss_gradient(
    values: &[f64],
    cell_centres: &[[f64; 3]],
    faces: &[Vec<usize>],
    points: &[[f64; 3]],
    owner: &[usize],
    neighbour: &[i64],
    n_cells: usize,
) -> Vec<[f64; 3]> {
    let mut grad = vec![[0.0_f64; 3]; n_cells];
    let mut vol = vec![0.0_f64; n_cells];
    for face_idx in 0..faces.len() {
        let face_verts: Vec<[f64; 3]> = faces[face_idx]
            .iter()
            .filter_map(|&vi| points.get(vi).copied())
            .collect();
        if face_verts.len() < 3 {
            continue;
        }
        let area = face_area(&face_verts);
        let norm = face_normal(&face_verts);
        let sf = [norm[0] * area, norm[1] * area, norm[2] * area];
        let fc = face_centroid(&face_verts);
        let own = owner[face_idx];
        let nbr = neighbour[face_idx];
        let phi_face = if nbr >= 0 {
            let n_idx = nbr as usize;
            let phi_own = values.get(own).copied().unwrap_or(0.0);
            let phi_nbr = values.get(n_idx).copied().unwrap_or(0.0);
            0.5 * (phi_own + phi_nbr)
        } else {
            values.get(own).copied().unwrap_or(0.0)
        };
        if own < n_cells {
            grad[own][0] += phi_face * sf[0];
            grad[own][1] += phi_face * sf[1];
            grad[own][2] += phi_face * sf[2];
            let cc = cell_centres.get(own).copied().unwrap_or([0.0; 3]);
            vol[own] +=
                ((fc[0] - cc[0]) * sf[0] + (fc[1] - cc[1]) * sf[1] + (fc[2] - cc[2]) * sf[2]).abs();
        }
        if nbr >= 0 {
            let n_idx = nbr as usize;
            if n_idx < n_cells {
                grad[n_idx][0] -= phi_face * sf[0];
                grad[n_idx][1] -= phi_face * sf[1];
                grad[n_idx][2] -= phi_face * sf[2];
                let cc = cell_centres.get(n_idx).copied().unwrap_or([0.0; 3]);
                vol[n_idx] +=
                    ((fc[0] - cc[0]) * sf[0] + (fc[1] - cc[1]) * sf[1] + (fc[2] - cc[2]) * sf[2])
                        .abs();
            }
        }
    }
    for i in 0..n_cells {
        let v = vol[i];
        if v > 1e-30 {
            grad[i][0] /= v;
            grad[i][1] /= v;
            grad[i][2] /= v;
        }
    }
    grad
}
/// k-epsilon turbulence model initial conditions helper.
///
/// Estimates initial `k` and `epsilon` from turbulence intensity `I` and
/// a length scale `L` at a reference velocity `U_ref`.
pub fn k_epsilon_initial(u_ref: f64, turbulence_intensity: f64, length_scale: f64) -> (f64, f64) {
    pub(super) const C_MU: f64 = 0.09;
    let k = 1.5 * (u_ref * turbulence_intensity).powi(2);
    let epsilon = C_MU.powf(0.75) * k.powf(1.5) / length_scale;
    (k, epsilon)
}
/// k-omega SST turbulence model initial conditions helper.
///
/// Returns `(k, omega)` from turbulence intensity and length scale.
pub fn k_omega_initial(u_ref: f64, turbulence_intensity: f64, length_scale: f64) -> (f64, f64) {
    pub(super) const C_MU: f64 = 0.09;
    let k = 1.5 * (u_ref * turbulence_intensity).powi(2);
    let omega = k.sqrt() / (C_MU.powf(0.25) * length_scale);
    (k, omega)
}
/// Compute turbulent viscosity from k-epsilon model.
///
/// `nu_t = C_mu * k^2 / epsilon`
pub fn turbulent_viscosity_ke(k: f64, epsilon: f64) -> f64 {
    pub(super) const C_MU: f64 = 0.09;
    if epsilon < 1e-30 {
        return 0.0;
    }
    C_MU * k * k / epsilon
}
/// Compute turbulent viscosity from k-omega SST model.
///
/// `nu_t = k / omega`
pub fn turbulent_viscosity_ko(k: f64, omega: f64) -> f64 {
    if omega < 1e-30 {
        return 0.0;
    }
    k / omega
}
/// Find the cell index in a structured box mesh that contains point `p`.
///
/// Works for `FoamMesh::box_mesh`-style meshes with uniform spacing.
/// Returns `None` if the point is outside the domain.
pub fn find_cell_containing_point(
    p: [f64; 3],
    lx: f64,
    ly: f64,
    lz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Option<usize> {
    if p[0] < 0.0 || p[0] >= lx || p[1] < 0.0 || p[1] >= ly || p[2] < 0.0 || p[2] >= lz {
        return None;
    }
    let i = ((p[0] / lx) * nx as f64) as usize;
    let j = ((p[1] / ly) * ny as f64) as usize;
    let k = ((p[2] / lz) * nz as f64) as usize;
    let i = i.min(nx - 1);
    let j = j.min(ny - 1);
    let k = k.min(nz - 1);
    Some(k * ny * nx + j * nx + i)
}
/// Probe a scalar field at a given cell index.
///
/// Returns the field value at the cell, or `None` if `cell_idx` is out of range.
pub fn probe_scalar(values: &[f64], cell_idx: usize) -> Option<f64> {
    values.get(cell_idx).copied()
}
/// Probe a vector field at a given cell index.
pub fn probe_vector(values: &[[f64; 3]], cell_idx: usize) -> Option<[f64; 3]> {
    values.get(cell_idx).copied()
}
/// Compute the average cell volume for a structured box mesh.
///
/// For a uniform mesh, this is simply `(lx*ly*lz) / (nx*ny*nz)`.
pub fn average_cell_volume(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> f64 {
    if nx == 0 || ny == 0 || nz == 0 {
        return 0.0;
    }
    (lx * ly * lz) / (nx * ny * nz) as f64
}
/// Compute the maximum face aspect ratio (longest edge / shortest edge) for a quad face.
///
/// Returns 1.0 for a perfect square face.
pub fn quad_face_aspect_ratio(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
    let d01 = edge_len(v0, v1);
    let d12 = edge_len(v1, v2);
    let d23 = edge_len(v2, v3);
    let d30 = edge_len(v3, v0);
    let max_edge = d01.max(d12).max(d23).max(d30);
    let min_edge = d01.min(d12).min(d23).min(d30);
    if min_edge < 1e-30 {
        return f64::INFINITY;
    }
    max_edge / min_edge
}
pub(super) fn edge_len(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
