// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! OpenFOAM case-directory format reader/writer.
//!
//! Supports parsing and writing scalar and vector fields in the OpenFOAM
//! dictionary format, as well as reading a simplified polyMesh directory
//! (points, faces, owner, neighbour).
//!
//! # Format overview
//!
//! An OpenFOAM field file is a dictionary with an `internalField` section
//! and a `boundaryField` section containing named patches.
//!
//! ```text
//! // * * * *
//! dimensions      [0 0 0 0 0 0 0];
//!
//! internalField   nonuniform List`scalar`
//! 3
//! (
//! 1.0
//! 2.0
//! 3.0
//! );
//!
//! boundaryField
//! {
//!     inlet
//!     {
//!         type            fixedValue;
//!         value           uniform 0;
//!     }
//! }
//! ```

use std::fmt::Write as FmtWrite;

use crate::Error as IoError;

// ── Boundary patch ────────────────────────────────────────────────────────────

/// An OpenFOAM boundary patch with a name and type string.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamBoundaryPatch {
    /// Patch name (e.g. `"inlet"`, `"outlet"`, `"wall"`).
    pub name: String,
    /// Patch type string (e.g. `"fixedValue"`, `"zeroGradient"`).
    pub patch_type: String,
    /// Uniform scalar value for this patch (if applicable).
    pub uniform_value: Option<f64>,
}

/// Create a named boundary patch.
///
/// This is a convenience constructor for [`FoamBoundaryPatch`].
pub fn foam_boundary_patch(
    name: impl Into<String>,
    patch_type: impl Into<String>,
    uniform_value: Option<f64>,
) -> FoamBoundaryPatch {
    FoamBoundaryPatch {
        name: name.into(),
        patch_type: patch_type.into(),
        uniform_value,
    }
}

// ── Field ─────────────────────────────────────────────────────────────────────

/// An OpenFOAM scalar or vector field.
///
/// The internal field holds one value per cell. Boundary patches provide
/// boundary conditions on named patch faces.
#[derive(Debug, Clone)]
pub struct FoamField {
    /// Field name (e.g. `"p"`, `"U"`).
    pub name: String,
    /// Internal field values (one entry per cell for scalar,
    /// three entries per cell for vector stored as `[vx0, vy0, vz0, …]`).
    pub internal_field: Vec<f64>,
    /// Number of components (1 for scalar, 3 for vector).
    pub n_components: usize,
    /// Boundary patches.
    pub patches: Vec<FoamBoundaryPatch>,
}

impl FoamField {
    /// Create a new scalar field with the given internal values.
    pub fn new_scalar(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            internal_field: values,
            n_components: 1,
            patches: Vec::new(),
        }
    }

    /// Create a new vector field. `values` must be a flat `[vx, vy, vz, …]` array.
    pub fn new_vector(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            internal_field: values,
            n_components: 3,
            patches: Vec::new(),
        }
    }

    /// Number of cells (logical entries) in the internal field.
    pub fn n_cells(&self) -> usize {
        self.internal_field
            .len()
            .checked_div(self.n_components)
            .unwrap_or(0)
    }
}

// ── Parser helpers ────────────────────────────────────────────────────────────

/// Extract the content of a `(…)` block from a whitespace-normalised string.
///
/// Returns the text between the first `(` and the matching `)`.
fn extract_paren_block(s: &str) -> Option<&str> {
    let start = s.find('(')?;
    let end = s.rfind(')')?;
    if end > start {
        Some(&s[start + 1..end])
    } else {
        None
    }
}

/// Parse a `nonuniform List`scalar` block.
fn parse_scalar_list(list_text: &str) -> Result<Vec<f64>, IoError> {
    let body = extract_paren_block(list_text)
        .ok_or_else(|| IoError::Parse("missing ( ) in scalar list".into()))?;
    body.split_whitespace()
        .filter(|s| !s.is_empty())
        .map(|tok| {
            tok.parse::<f64>()
                .map_err(|e| IoError::Parse(e.to_string()))
        })
        .collect()
}

/// Parse a `nonuniform List`vector` block — each entry is `(vx vy vz)`.
fn parse_vector_list(list_text: &str) -> Result<Vec<f64>, IoError> {
    let body = extract_paren_block(list_text)
        .ok_or_else(|| IoError::Parse("missing outer ( ) in vector list".into()))?;
    let mut values: Vec<f64> = Vec::new();
    let mut remaining = body.trim();
    while !remaining.is_empty() {
        if let Some(open) = remaining.find('(') {
            let close = remaining[open..].find(')').map(|i| open + i);
            if let Some(close) = close {
                let triple = &remaining[open + 1..close];
                for tok in triple.split_whitespace() {
                    values.push(
                        tok.parse::<f64>()
                            .map_err(|e| IoError::Parse(e.to_string()))?,
                    );
                }
                remaining = &remaining[close + 1..];
            } else {
                break;
            }
        } else {
            break;
        }
    }
    Ok(values)
}

// ── Public parsers ────────────────────────────────────────────────────────────

/// Parse an OpenFOAM scalar field file (simplified).
///
/// Recognises the `internalField nonuniform List`scalar` section.
///
/// # Errors
/// Returns `Error::Parse` on malformed input.
pub fn parse_foam_scalar(src: &str, name: impl Into<String>) -> Result<FoamField, IoError> {
    // Find the List<scalar> block
    let marker = "List<scalar>";
    if let Some(pos) = src.find(marker) {
        let after = &src[pos + marker.len()..];
        // skip the count line then read (…)
        let values = parse_scalar_list(after)?;
        let mut field = FoamField::new_scalar(name, values);
        // parse simple boundary patches
        field.patches = parse_boundary_patches(src);
        return Ok(field);
    }
    // Try uniform value
    let uniform_marker = "internalField   uniform";
    if let Some(pos) = src.find(uniform_marker) {
        let after = &src[pos + uniform_marker.len()..];
        let tok = after
            .split_whitespace()
            .next()
            .ok_or_else(|| IoError::Parse("missing uniform value".into()))?;
        let val: f64 = tok
            .trim_end_matches(';')
            .parse()
            .map_err(|e: std::num::ParseFloatError| IoError::Parse(e.to_string()))?;
        let n_cells = detect_n_cells(src);
        let mut field = FoamField::new_scalar(name, vec![val; n_cells]);
        field.patches = parse_boundary_patches(src);
        return Ok(field);
    }
    Err(IoError::Parse(
        "could not find internalField in foam scalar file".into(),
    ))
}

/// Parse an OpenFOAM vector field file (simplified).
///
/// Recognises the `internalField nonuniform List`vector` section.
///
/// # Errors
/// Returns `Error::Parse` on malformed input.
pub fn parse_foam_vector(src: &str, name: impl Into<String>) -> Result<FoamField, IoError> {
    let marker = "List<vector>";
    if let Some(pos) = src.find(marker) {
        let after = &src[pos + marker.len()..];
        let values = parse_vector_list(after)?;
        let mut field = FoamField::new_vector(name, values);
        field.patches = parse_boundary_patches(src);
        return Ok(field);
    }
    Err(IoError::Parse(
        "could not find List<vector> in foam vector file".into(),
    ))
}

/// Detect number of cells from `// nCells: N` comment if present, otherwise 0.
fn detect_n_cells(src: &str) -> usize {
    for line in src.lines() {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("// nCells:")
            && let Ok(n) = rest.trim().parse::<usize>()
        {
            return n;
        }
    }
    1
}

/// Parse boundary patches from an OpenFOAM field text.
///
/// Very simplified: finds `patchName { type X; value uniform V; }` blocks.
fn parse_boundary_patches(src: &str) -> Vec<FoamBoundaryPatch> {
    let mut patches = Vec::new();
    // find boundaryField {
    let bf_marker = "boundaryField";
    let Some(bf_pos) = src.find(bf_marker) else {
        return patches;
    };
    let after_bf = &src[bf_pos + bf_marker.len()..];
    // find first {
    let Some(open) = after_bf.find('{') else {
        return patches;
    };
    let inner = &after_bf[open + 1..];
    // walk through each patch block: "name { type X; }"
    let mut cursor = inner;
    loop {
        cursor = cursor.trim_start();
        if cursor.is_empty() || cursor.starts_with('}') {
            break;
        }
        // patch name is the first token
        let tok_end = cursor
            .find(|c: char| c.is_whitespace() || c == '{')
            .unwrap_or(cursor.len());
        let patch_name = cursor[..tok_end].trim().to_string();
        if patch_name.is_empty() || patch_name == "}" {
            break;
        }
        // advance past the name
        cursor = &cursor[tok_end..];
        // find opening brace
        let Some(blk_open) = cursor.find('{') else {
            break;
        };
        let blk_start = blk_open + 1;
        let Some(blk_close) = cursor[blk_start..].find('}') else {
            break;
        };
        let block = &cursor[blk_start..blk_start + blk_close];
        cursor = &cursor[blk_start + blk_close + 1..];
        // extract type
        let patch_type = block
            .lines()
            .find_map(|l| {
                let l = l.trim();
                if l.starts_with("type") {
                    Some(
                        l.trim_start_matches("type")
                            .trim()
                            .trim_end_matches(';')
                            .trim()
                            .to_string(),
                    )
                } else {
                    None
                }
            })
            .unwrap_or_else(|| "unknown".to_string());
        // extract uniform value
        let uniform_value = block.lines().find_map(|l| {
            let l = l.trim();
            if l.starts_with("value") && l.contains("uniform") {
                let parts: Vec<&str> = l.split_whitespace().collect();
                // "value  uniform 0;"
                if parts.len() >= 3 {
                    parts[2].trim_end_matches(';').parse::<f64>().ok()
                } else {
                    None
                }
            } else {
                None
            }
        });
        patches.push(FoamBoundaryPatch {
            name: patch_name,
            patch_type,
            uniform_value,
        });
    }
    patches
}

/// Write a [`FoamField`] in OpenFOAM dictionary format.
pub fn write_foam_field(field: &FoamField) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "// OpenFOAM field: {}", field.name);
    let _ = writeln!(out, "dimensions      [0 0 0 0 0 0 0];");
    let _ = writeln!(out);
    if field.n_components == 1 {
        let _ = writeln!(out, "internalField   nonuniform List<scalar>");
        let _ = writeln!(out, "{}", field.n_cells());
        let _ = writeln!(out, "(");
        for &v in &field.internal_field {
            let _ = writeln!(out, "{:.15e}", v);
        }
        let _ = writeln!(out, ");");
    } else {
        let _ = writeln!(out, "internalField   nonuniform List<vector>");
        let n = field.n_cells();
        let _ = writeln!(out, "{}", n);
        let _ = writeln!(out, "(");
        for i in 0..n {
            let base = i * 3;
            let (vx, vy, vz) = (
                field.internal_field[base],
                field.internal_field[base + 1],
                field.internal_field[base + 2],
            );
            let _ = writeln!(out, "({:.15e} {:.15e} {:.15e})", vx, vy, vz);
        }
        let _ = writeln!(out, ");");
    }
    let _ = writeln!(out);
    let _ = writeln!(out, "boundaryField");
    let _ = writeln!(out, "{{");
    for p in &field.patches {
        let _ = writeln!(out, "    {}", p.name);
        let _ = writeln!(out, "    {{");
        let _ = writeln!(out, "        type            {};", p.patch_type);
        if let Some(val) = p.uniform_value {
            let _ = writeln!(out, "        value           uniform {:.15e};", val);
        }
        let _ = writeln!(out, "    }}");
    }
    let _ = writeln!(out, "}}");
    out
}

// ── Mesh ──────────────────────────────────────────────────────────────────────

/// An OpenFOAM polyMesh in memory.
///
/// This is a simplified in-memory representation of the files found in
/// `constant/polyMesh/`.
#[derive(Debug, Clone)]
pub struct FoamMesh {
    /// Points: flat `[x0, y0, z0, x1, y1, z1, …]`.
    pub points: Vec<f64>,
    /// Faces: each face is a list of point indices.
    pub faces: Vec<Vec<usize>>,
    /// Owner cell index for each face.
    pub owner: Vec<usize>,
    /// Neighbour cell index for each internal face.
    pub neighbour: Vec<usize>,
}

impl FoamMesh {
    /// Number of mesh points.
    pub fn n_points(&self) -> usize {
        self.points.len() / 3
    }

    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of internal faces (those with a neighbour).
    pub fn n_internal_faces(&self) -> usize {
        self.neighbour.len()
    }
}

/// Parse a simplified OpenFOAM `points` file (ASCII).
///
/// Format: `n`\n(\n(x y z)\n…\n)`
fn parse_foam_points(src: &str) -> Result<Vec<f64>, IoError> {
    let body = extract_paren_block(src)
        .ok_or_else(|| IoError::Parse("missing ( ) in points file".into()))?;
    let mut pts: Vec<f64> = Vec::new();
    for chunk in body.split('(').skip(1) {
        let close = chunk
            .find(')')
            .ok_or_else(|| IoError::Parse("missing ) in point".into()))?;
        let triple = &chunk[..close];
        for tok in triple.split_whitespace() {
            pts.push(
                tok.parse::<f64>()
                    .map_err(|e| IoError::Parse(e.to_string()))?,
            );
        }
    }
    Ok(pts)
}

/// Parse a simplified OpenFOAM `owner` or `neighbour` file (ASCII).
fn parse_foam_index_list(src: &str) -> Result<Vec<usize>, IoError> {
    let body = extract_paren_block(src)
        .ok_or_else(|| IoError::Parse("missing ( ) in index list".into()))?;
    body.split_whitespace()
        .filter(|s| !s.is_empty())
        .map(|tok| {
            tok.parse::<usize>()
                .map_err(|e| IoError::Parse(e.to_string()))
        })
        .collect()
}

/// Parse a simplified OpenFOAM `faces` file (ASCII).
///
/// Each face entry: `n(p0 p1 … pn-1)`
fn parse_foam_faces(src: &str) -> Result<Vec<Vec<usize>>, IoError> {
    let body = extract_paren_block(src)
        .ok_or_else(|| IoError::Parse("missing outer ( ) in faces file".into()))?;
    let mut faces: Vec<Vec<usize>> = Vec::new();
    let mut cursor = body.trim();
    loop {
        cursor = cursor.trim_start();
        if cursor.is_empty() {
            break;
        }
        // The face entry starts with an optional count followed by '('
        let Some(open) = cursor.find('(') else { break };
        let close = cursor[open..].find(')').map(|i| open + i);
        let Some(close) = close else { break };
        let face_str = &cursor[open + 1..close];
        let pts: Vec<usize> = face_str
            .split_whitespace()
            .map(|t| {
                t.parse::<usize>()
                    .map_err(|e| IoError::Parse(e.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        faces.push(pts);
        cursor = &cursor[close + 1..];
    }
    Ok(faces)
}

/// Read an OpenFOAM polyMesh directory from four string slices.
///
/// All four string arguments correspond to the ASCII content of the
/// `points`, `faces`, `owner`, and `neighbour` files respectively.
///
/// # Errors
/// Returns `Error::Parse` on malformed input.
pub fn read_foam_mesh(
    points_src: &str,
    faces_src: &str,
    owner_src: &str,
    neighbour_src: &str,
) -> Result<FoamMesh, IoError> {
    let points = parse_foam_points(points_src)?;
    let faces = parse_foam_faces(faces_src)?;
    let owner = parse_foam_index_list(owner_src)?;
    let neighbour = parse_foam_index_list(neighbour_src)?;
    Ok(FoamMesh {
        points,
        faces,
        owner,
        neighbour,
    })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const SCALAR_FILE: &str = r#"
// OpenFOAM field
dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
3
(
1.0
2.5
3.0
);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0;
    }
    outlet
    {
        type            zeroGradient;
    }
}
"#;

    const VECTOR_FILE: &str = r#"
// OpenFOAM vector field
dimensions      [0 1 -1 0 0 0 0];

internalField   nonuniform List<vector>
2
(
(1.0 2.0 3.0)
(4.0 5.0 6.0)
);

boundaryField
{
    wall
    {
        type            noSlip;
    }
}
"#;

    const POINTS_FILE: &str = r#"
3
(
(0.0 0.0 0.0)
(1.0 0.0 0.0)
(0.5 1.0 0.0)
)
"#;

    const FACES_FILE: &str = r#"
1
(
3(0 1 2)
)
"#;

    const OWNER_FILE: &str = r#"
1
(
0
)
"#;

    const NEIGHBOUR_FILE: &str = r#"
0
(
)
"#;

    #[test]
    fn test_parse_scalar_n_cells() {
        let f = parse_foam_scalar(SCALAR_FILE, "p").unwrap();
        assert_eq!(f.n_cells(), 3);
    }

    #[test]
    fn test_parse_scalar_values() {
        let f = parse_foam_scalar(SCALAR_FILE, "p").unwrap();
        assert!((f.internal_field[0] - 1.0).abs() < 1e-10);
        assert!((f.internal_field[1] - 2.5).abs() < 1e-10);
        assert!((f.internal_field[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_scalar_name() {
        let f = parse_foam_scalar(SCALAR_FILE, "pressure").unwrap();
        assert_eq!(f.name, "pressure");
    }

    #[test]
    fn test_parse_scalar_n_components() {
        let f = parse_foam_scalar(SCALAR_FILE, "p").unwrap();
        assert_eq!(f.n_components, 1);
    }

    #[test]
    fn test_parse_scalar_patches() {
        let f = parse_foam_scalar(SCALAR_FILE, "p").unwrap();
        assert_eq!(f.patches.len(), 2);
        assert_eq!(f.patches[0].name, "inlet");
        assert_eq!(f.patches[0].patch_type, "fixedValue");
        assert_eq!(f.patches[0].uniform_value, Some(0.0));
    }

    #[test]
    fn test_parse_scalar_outlet_no_value() {
        let f = parse_foam_scalar(SCALAR_FILE, "p").unwrap();
        assert_eq!(f.patches[1].name, "outlet");
        assert!(f.patches[1].uniform_value.is_none());
    }

    #[test]
    fn test_parse_scalar_error_on_empty() {
        assert!(parse_foam_scalar("", "p").is_err());
    }

    #[test]
    fn test_parse_vector_n_cells() {
        let f = parse_foam_vector(VECTOR_FILE, "U").unwrap();
        assert_eq!(f.n_cells(), 2);
    }

    #[test]
    fn test_parse_vector_n_components() {
        let f = parse_foam_vector(VECTOR_FILE, "U").unwrap();
        assert_eq!(f.n_components, 3);
    }

    #[test]
    fn test_parse_vector_values() {
        let f = parse_foam_vector(VECTOR_FILE, "U").unwrap();
        assert!((f.internal_field[0] - 1.0).abs() < 1e-10);
        assert!((f.internal_field[1] - 2.0).abs() < 1e-10);
        assert!((f.internal_field[2] - 3.0).abs() < 1e-10);
        assert!((f.internal_field[3] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_vector_error_on_empty() {
        assert!(parse_foam_vector("", "U").is_err());
    }

    #[test]
    fn test_write_scalar_roundtrip() {
        let original = FoamField {
            name: "p".to_string(),
            internal_field: vec![1.5, 2.5, 3.5],
            n_components: 1,
            patches: vec![foam_boundary_patch("inlet", "fixedValue", Some(0.0))],
        };
        let s = write_foam_field(&original);
        let parsed = parse_foam_scalar(&s, "p").unwrap();
        assert_eq!(parsed.n_cells(), 3);
        assert!((parsed.internal_field[0] - 1.5).abs() < 1e-10);
        assert!((parsed.internal_field[2] - 3.5).abs() < 1e-10);
    }

    #[test]
    fn test_write_vector_roundtrip() {
        let original = FoamField::new_vector("U", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let s = write_foam_field(&original);
        let parsed = parse_foam_vector(&s, "U").unwrap();
        assert_eq!(parsed.n_cells(), 2);
        for (i, &expected) in [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0].iter().enumerate() {
            assert!(
                (parsed.internal_field[i] - expected).abs() < 1e-10,
                "mismatch at {i}"
            );
        }
    }

    #[test]
    fn test_write_field_contains_boundary_field() {
        let f = FoamField::new_scalar("p", vec![0.0]);
        let s = write_foam_field(&f);
        assert!(s.contains("boundaryField"));
    }

    #[test]
    fn test_foam_boundary_patch_constructor() {
        let p = foam_boundary_patch("wall", "noSlip", None);
        assert_eq!(p.name, "wall");
        assert_eq!(p.patch_type, "noSlip");
        assert!(p.uniform_value.is_none());
    }

    #[test]
    fn test_foam_boundary_patch_with_value() {
        let p = foam_boundary_patch("inlet", "fixedValue", Some(1.0));
        assert_eq!(p.uniform_value, Some(1.0));
    }

    #[test]
    fn test_read_foam_mesh_n_points() {
        let m = read_foam_mesh(POINTS_FILE, FACES_FILE, OWNER_FILE, NEIGHBOUR_FILE).unwrap();
        assert_eq!(m.n_points(), 3);
    }

    #[test]
    fn test_read_foam_mesh_point_values() {
        let m = read_foam_mesh(POINTS_FILE, FACES_FILE, OWNER_FILE, NEIGHBOUR_FILE).unwrap();
        assert!((m.points[0] - 0.0).abs() < 1e-12);
        assert!((m.points[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_read_foam_mesh_n_faces() {
        let m = read_foam_mesh(POINTS_FILE, FACES_FILE, OWNER_FILE, NEIGHBOUR_FILE).unwrap();
        assert_eq!(m.n_faces(), 1);
    }

    #[test]
    fn test_read_foam_mesh_face_connectivity() {
        let m = read_foam_mesh(POINTS_FILE, FACES_FILE, OWNER_FILE, NEIGHBOUR_FILE).unwrap();
        assert_eq!(m.faces[0], vec![0, 1, 2]);
    }

    #[test]
    fn test_read_foam_mesh_owner() {
        let m = read_foam_mesh(POINTS_FILE, FACES_FILE, OWNER_FILE, NEIGHBOUR_FILE).unwrap();
        assert_eq!(m.owner, vec![0]);
    }

    #[test]
    fn test_read_foam_mesh_no_neighbours() {
        let m = read_foam_mesh(POINTS_FILE, FACES_FILE, OWNER_FILE, NEIGHBOUR_FILE).unwrap();
        assert_eq!(m.n_internal_faces(), 0);
    }

    #[test]
    fn test_foam_field_new_scalar() {
        let f = FoamField::new_scalar("T", vec![300.0, 310.0]);
        assert_eq!(f.n_components, 1);
        assert_eq!(f.n_cells(), 2);
    }

    #[test]
    fn test_foam_field_new_vector() {
        let f = FoamField::new_vector("U", vec![1.0, 0.0, 0.0, 2.0, 0.0, 0.0]);
        assert_eq!(f.n_components, 3);
        assert_eq!(f.n_cells(), 2);
    }

    #[test]
    fn test_foam_field_debug() {
        let f = FoamField::new_scalar("p", vec![]);
        let s = format!("{f:?}");
        assert!(s.contains("FoamField"));
    }

    #[test]
    fn test_foam_mesh_debug() {
        let m = FoamMesh {
            points: vec![],
            faces: vec![],
            owner: vec![],
            neighbour: vec![],
        };
        let s = format!("{m:?}");
        assert!(s.contains("FoamMesh"));
    }

    #[test]
    fn test_write_field_has_internal_field_keyword() {
        let f = FoamField::new_scalar("rho", vec![1.0]);
        let s = write_foam_field(&f);
        assert!(s.contains("internalField"));
    }

    #[test]
    fn test_parse_vector_patches() {
        let f = parse_foam_vector(VECTOR_FILE, "U").unwrap();
        assert_eq!(f.patches.len(), 1);
        assert_eq!(f.patches[0].name, "wall");
    }

    #[test]
    fn test_foam_boundary_patch_clone() {
        let p = foam_boundary_patch("inlet", "fixedValue", Some(5.0));
        let p2 = p.clone();
        assert_eq!(p2.name, "inlet");
    }

    #[test]
    fn test_foam_mesh_clone() {
        let m = FoamMesh {
            points: vec![0.0, 0.0, 0.0],
            faces: vec![vec![0]],
            owner: vec![0],
            neighbour: vec![],
        };
        let m2 = m.clone();
        assert_eq!(m2.n_points(), 1);
    }
}
