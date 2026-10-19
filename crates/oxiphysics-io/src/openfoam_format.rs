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
use crate::openfoam::strip_foam_comments;

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

/// A parenthesised OpenFOAM list `[N] ( ... )` located inside a larger text.
struct ListBlock<'a> {
    /// The size written in front of the opening parenthesis, when present.
    declared: Option<usize>,
    /// The text between the opening parenthesis and its matching closing one.
    body: &'a str,
    /// The text after the matching closing parenthesis.
    rest: &'a str,
}

impl ListBlock<'_> {
    /// Fail when the declared size differs from the number of entries found.
    fn check_count(&self, found: usize, what: &str) -> Result<(), IoError> {
        match self.declared {
            Some(declared) if declared != found => Err(IoError::Parse(format!(
                "{what} list declares {declared} entries but holds {found}"
            ))),
            _ => Ok(()),
        }
    }
}

/// Find the first `( ... )` block of `s`.
///
/// The closing parenthesis is the one that matches the first opening
/// parenthesis (nested ones are counted), so anything that follows the list,
/// such as a `boundaryField` section full of parentheses, is never swallowed
/// into it.  `what` names the list in error messages.
fn read_list_block<'a>(s: &'a str, what: &str) -> Result<ListBlock<'a>, IoError> {
    let open = s
        .find('(')
        .ok_or_else(|| IoError::Parse(format!("missing ( ) in {what}")))?;
    let mut depth = 0usize;
    let mut close = None;
    for (offset, c) in s[open..].char_indices() {
        match c {
            '(' => depth += 1,
            ')' => {
                depth = depth.saturating_sub(1);
                if depth == 0 {
                    close = Some(open + offset);
                    break;
                }
            }
            _ => {}
        }
    }
    let close = close.ok_or_else(|| IoError::Parse(format!("unbalanced parentheses in {what}")))?;
    Ok(ListBlock {
        declared: declared_size(&s[..open])?,
        body: &s[open + 1..close],
        rest: &s[close + 1..],
    })
}

/// The list size written directly in front of an opening parenthesis, such as
/// the `3` of `List<scalar> 3 (` or of `3(`.  Digits that end an identifier or
/// a number are not a size.
fn declared_size(before: &str) -> Result<Option<usize>, IoError> {
    let before = before.trim_end();
    let digits_start = before
        .char_indices()
        .rev()
        .take_while(|(_, c)| c.is_ascii_digit())
        .last()
        .map(|(index, _)| index);
    let Some(start) = digits_start else {
        return Ok(None);
    };
    if before[..start]
        .chars()
        .next_back()
        .is_some_and(|c| c.is_alphanumeric() || matches!(c, '_' | '.' | '-' | '+' | '$'))
    {
        return Ok(None);
    }
    before[start..]
        .parse::<usize>()
        .map(Some)
        .map_err(|e| IoError::Parse(format!("bad list size '{}': {e}", &before[start..])))
}

/// Parse a scalar list `[N] ( v0 v1 ... )`; the declared size must match.
fn parse_scalar_list(list_text: &str) -> Result<Vec<f64>, IoError> {
    let block = read_list_block(list_text, "scalar list")?;
    let values = block
        .body
        .split_whitespace()
        .map(|tok| {
            tok.parse::<f64>()
                .map_err(|e| IoError::Parse(format!("bad scalar '{tok}': {e}")))
        })
        .collect::<Result<Vec<f64>, IoError>>()?;
    block.check_count(values.len(), "scalar")?;
    Ok(values)
}

/// Parse consecutive `(x y z)` groups into a flat `[x, y, z, ...]` vector and
/// the number of groups.
fn parse_vector_groups(body: &str) -> Result<(Vec<f64>, usize), IoError> {
    let mut values = Vec::new();
    let mut count = 0usize;
    let mut remaining = body.trim();
    while !remaining.is_empty() {
        let inner = remaining.strip_prefix('(').ok_or_else(|| {
            IoError::Parse(format!(
                "expected '(' at the start of a vector, found '{}'",
                remaining.chars().take(20).collect::<String>()
            ))
        })?;
        let close = inner
            .find(')')
            .ok_or_else(|| IoError::Parse("unterminated vector, missing ')'".into()))?;
        let components = inner[..close]
            .split_whitespace()
            .map(|tok| {
                tok.parse::<f64>()
                    .map_err(|e| IoError::Parse(format!("bad vector component '{tok}': {e}")))
            })
            .collect::<Result<Vec<f64>, IoError>>()?;
        if components.len() != 3 {
            return Err(IoError::Parse(format!(
                "a vector needs 3 components, found {}",
                components.len()
            )));
        }
        values.extend_from_slice(&components);
        count += 1;
        remaining = inner[close + 1..].trim_start();
    }
    Ok((values, count))
}

/// Parse a vector list `[N] ( (x y z) ... )` into a flat `[x, y, z, ...]`
/// vector; the declared size must match.
fn parse_vector_list(list_text: &str) -> Result<Vec<f64>, IoError> {
    let block = read_list_block(list_text, "vector list")?;
    let (values, count) = parse_vector_groups(block.body)?;
    block.check_count(count, "vector")?;
    Ok(values)
}

/// The text that follows the `internalField` keyword, or `None` when the
/// keyword is absent.  Mentions such as `$internalField` or a longer
/// identifier do not count.
fn internal_field_entry(src: &str) -> Option<&str> {
    const KEY: &str = "internalField";
    src.match_indices(KEY).find_map(|(pos, _)| {
        let before_ok = src[..pos]
            .chars()
            .next_back()
            .is_none_or(|c| !(c.is_alphanumeric() || matches!(c, '_' | '$' | ':')));
        let after = &src[pos + KEY.len()..];
        let after_ok = after.chars().next().is_none_or(char::is_whitespace);
        (before_ok && after_ok).then_some(after)
    })
}

/// Split the leading word (ended by whitespace, `(` or `;`) off `s`.
fn split_keyword(s: &str) -> (&str, &str) {
    let s = s.trim_start();
    let end = s
        .find(|c: char| c.is_whitespace() || c == '(' || c == ';')
        .unwrap_or(s.len());
    s.split_at(end)
}

/// The text of a `nonuniform` entry starting at its list: skips the optional
/// `List<element>` marker and refuses to run on into later entries.
fn nonuniform_list<'a>(rest: &'a str, element: &str) -> Result<&'a str, IoError> {
    let rest = rest.trim_start();
    let marker = format!("List<{element}>");
    let list_text = if rest.starts_with("List<") {
        rest.strip_prefix(marker.as_str())
            .ok_or_else(|| IoError::Parse(format!("expected a {marker} internalField")))?
    } else {
        rest
    };
    match (list_text.find('('), list_text.find(';')) {
        (None, _) => Err(IoError::Parse(
            "nonuniform internalField without a list".into(),
        )),
        (Some(open), Some(end)) if end < open => Err(IoError::Parse(
            "nonuniform internalField without a list".into(),
        )),
        _ => Ok(list_text),
    }
}

/// Largest cell count a uniform `internalField` is expanded to.
const MAX_UNIFORM_CELLS: usize = 1 << 26;

/// Repeat the components of a uniform value once per cell.
fn expand_uniform(components: &[f64], n_cells: usize) -> Result<Vec<f64>, IoError> {
    if n_cells > MAX_UNIFORM_CELLS {
        return Err(IoError::Parse(format!(
            "cannot expand a uniform field to {n_cells} cells (limit {MAX_UNIFORM_CELLS})"
        )));
    }
    let mut values = Vec::with_capacity(n_cells * components.len());
    for _ in 0..n_cells {
        values.extend_from_slice(components);
    }
    Ok(values)
}

/// Values of the `internalField` entry `entry` of a scalar field.
fn parse_scalar_internal_field(entry: &str, src: &str) -> Result<Vec<f64>, IoError> {
    let (keyword, rest) = split_keyword(entry);
    match keyword {
        "uniform" => {
            let token = rest
                .split(|c: char| c.is_whitespace() || c == ';')
                .find(|t| !t.is_empty())
                .ok_or_else(|| IoError::Parse("missing uniform value".into()))?;
            let value: f64 = token
                .parse()
                .map_err(|e: std::num::ParseFloatError| IoError::Parse(e.to_string()))?;
            expand_uniform(&[value], detect_n_cells(src))
        }
        "nonuniform" => parse_scalar_list(nonuniform_list(rest, "scalar")?),
        other => Err(IoError::Parse(format!(
            "unsupported internalField form '{other}'"
        ))),
    }
}

/// Values (flat `[x, y, z, ...]`) of the `internalField` entry `entry` of a
/// vector field.
fn parse_vector_internal_field(entry: &str, src: &str) -> Result<Vec<f64>, IoError> {
    let (keyword, rest) = split_keyword(entry);
    match keyword {
        "uniform" => {
            let rest = rest.trim_start();
            if !rest.starts_with('(') {
                return Err(IoError::Parse(
                    "a uniform vector internalField needs (x y z)".into(),
                ));
            }
            let block = read_list_block(rest, "uniform vector")?;
            let components = block
                .body
                .split_whitespace()
                .map(|tok| {
                    tok.parse::<f64>()
                        .map_err(|e| IoError::Parse(format!("bad vector component '{tok}': {e}")))
                })
                .collect::<Result<Vec<f64>, IoError>>()?;
            if components.len() != 3 {
                return Err(IoError::Parse(format!(
                    "a vector needs 3 components, found {}",
                    components.len()
                )));
            }
            expand_uniform(&components, detect_n_cells(src))
        }
        "nonuniform" => parse_vector_list(nonuniform_list(rest, "vector")?),
        other => Err(IoError::Parse(format!(
            "unsupported internalField form '{other}'"
        ))),
    }
}

// ── Public parsers ────────────────────────────────────────────────────────────

/// Parse an OpenFOAM scalar field file (simplified).
///
/// Recognises the `internalField` entry in its `uniform <value>` and
/// `nonuniform List<scalar> N ( ... )` forms (the size may be glued to the
/// parenthesis, `N( ... )`).  Only the `internalField` entry is read, so
/// parentheses and lists in the `boundaryField` section never reach the cell
/// values; the declared list size must match the number of values.
///
/// # Errors
/// Returns `Error::Parse` on malformed input.
pub fn parse_foam_scalar(src: &str, name: impl Into<String>) -> Result<FoamField, IoError> {
    let cleaned = strip_foam_comments(src);
    let values = match internal_field_entry(&cleaned) {
        Some(entry) => parse_scalar_internal_field(entry, src)?,
        None => {
            // Bare snippet without the keyword: read the first scalar list.
            let marker = "List<scalar>";
            let pos = cleaned.find(marker).ok_or_else(|| {
                IoError::Parse("could not find internalField in foam scalar file".into())
            })?;
            parse_scalar_list(&cleaned[pos + marker.len()..])?
        }
    };
    let mut field = FoamField::new_scalar(name, values);
    field.patches = parse_boundary_patches(src);
    Ok(field)
}

/// Parse an OpenFOAM vector field file (simplified).
///
/// Recognises the `internalField` entry in its `uniform (x y z)` and
/// `nonuniform List<vector> N ( (x y z) ... )` forms.  Only the
/// `internalField` entry is read, so the vectors of `value uniform (...)`
/// entries in the `boundaryField` section are never mistaken for cells; the
/// declared list size must match the number of vectors.
///
/// # Errors
/// Returns `Error::Parse` on malformed input.
pub fn parse_foam_vector(src: &str, name: impl Into<String>) -> Result<FoamField, IoError> {
    let cleaned = strip_foam_comments(src);
    let values = match internal_field_entry(&cleaned) {
        Some(entry) => parse_vector_internal_field(entry, src)?,
        None => {
            // Bare snippet without the keyword: read the first vector list.
            let marker = "List<vector>";
            let pos = cleaned.find(marker).ok_or_else(|| {
                IoError::Parse("could not find internalField in foam vector file".into())
            })?;
            parse_vector_list(&cleaned[pos + marker.len()..])?
        }
    };
    let mut field = FoamField::new_vector(name, values);
    field.patches = parse_boundary_patches(src);
    Ok(field)
}

/// Detect the number of cells from a `// nCells: N` comment.
///
/// A field file does not store its cell count, so without that comment a
/// `uniform` internal field is represented by a single cell (`1`).
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
    let cleaned = strip_foam_comments(src);
    let block = read_list_block(&cleaned, "points file")?;
    let (points, count) = parse_vector_groups(block.body)?;
    block.check_count(count, "points")?;
    Ok(points)
}

/// Parse a simplified OpenFOAM `owner` or `neighbour` file (ASCII).
fn parse_foam_index_list(src: &str) -> Result<Vec<usize>, IoError> {
    let cleaned = strip_foam_comments(src);
    let block = read_list_block(&cleaned, "index list")?;
    let indices = block
        .body
        .split_whitespace()
        .map(|tok| {
            tok.parse::<usize>()
                .map_err(|e| IoError::Parse(format!("bad cell index '{tok}': {e}")))
        })
        .collect::<Result<Vec<usize>, IoError>>()?;
    block.check_count(indices.len(), "index")?;
    Ok(indices)
}

/// Parse a simplified OpenFOAM `faces` file (ASCII).
///
/// Each face entry: `n(p0 p1 … pn-1)`
fn parse_foam_faces(src: &str) -> Result<Vec<Vec<usize>>, IoError> {
    let cleaned = strip_foam_comments(src);
    let outer = read_list_block(&cleaned, "faces file")?;
    let mut faces: Vec<Vec<usize>> = Vec::new();
    let mut cursor = outer.body.trim();
    while !cursor.is_empty() {
        let open = cursor
            .find('(')
            .ok_or_else(|| IoError::Parse("missing ( ) in face".into()))?;
        let prefix = cursor[..open].trim();
        if !prefix.bytes().all(|b| b.is_ascii_digit()) {
            return Err(IoError::Parse(format!(
                "unexpected text '{prefix}' in the faces list"
            )));
        }
        let face = read_list_block(cursor, "face")?;
        let vertices = face
            .body
            .split_whitespace()
            .map(|tok| {
                tok.parse::<usize>()
                    .map_err(|e| IoError::Parse(format!("bad vertex index '{tok}': {e}")))
            })
            .collect::<Result<Vec<usize>, IoError>>()?;
        face.check_count(vertices.len(), "face")?;
        faces.push(vertices);
        cursor = face.rest.trim_start();
    }
    outer.check_count(faces.len(), "faces")?;
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
mod list_tests;

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
