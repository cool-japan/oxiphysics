// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simplified OpenFOAM polyMesh format I/O.
//!
//! Supports reading polyMesh directories and writing basic OpenFOAM cases.
//! Uses plain `f64` arrays and `Vec` — no external linear-algebra crates.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// FoamHeader
// ---------------------------------------------------------------------------

/// FoamFile header block, common to every OpenFOAM dictionary file.
pub struct FoamHeader {
    /// OpenFOAM dictionary version, e.g. `"2.0"`.
    pub foam_file_version: String,
    /// Data format: `"ascii"` or `"binary"`.
    pub format: String,
    /// Dictionary class, e.g. `"polyBoundaryMesh"`.
    pub class_name: String,
    /// Object name, e.g. `"boundary"`.
    pub object_name: String,
}

impl FoamHeader {
    /// Parse a string that contains a `FoamFile { … }` block.
    pub fn parse(s: &str) -> Result<Self, String> {
        let start = s
            .find("FoamFile")
            .ok_or_else(|| "No FoamFile block found".to_string())?;
        let brace_open = s[start..]
            .find('{')
            .ok_or_else(|| "No opening brace for FoamFile".to_string())?
            + start;
        let brace_close = s[brace_open..]
            .find('}')
            .ok_or_else(|| "No closing brace for FoamFile".to_string())?
            + brace_open;

        let block = &s[brace_open + 1..brace_close];

        let mut kv: HashMap<&str, &str> = HashMap::new();
        for line in block.lines() {
            let line = line.trim().trim_end_matches(';');
            let mut parts = line.splitn(2, char::is_whitespace);
            if let (Some(k), Some(v)) = (parts.next(), parts.next()) {
                kv.insert(k.trim(), v.trim());
            }
        }

        Ok(FoamHeader {
            foam_file_version: kv.get("version").unwrap_or(&"2.0").to_string(),
            format: kv.get("format").unwrap_or(&"ascii").to_string(),
            class_name: kv.get("class").unwrap_or(&"").to_string(),
            object_name: kv.get("object").unwrap_or(&"").to_string(),
        })
    }
}

impl std::fmt::Display for FoamHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "FoamFile\n{{\n    version     {};\n    format      {};\n    class       {};\n    object      {};\n}}\n",
            self.foam_file_version, self.format, self.class_name, self.object_name,
        )
    }
}

impl std::str::FromStr for FoamHeader {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ---------------------------------------------------------------------------
// Basic mesh types
// ---------------------------------------------------------------------------

/// A point in 3-D space.
pub type FoamPoint = [f64; 3];

/// A polygonal face described by an ordered list of node indices.
pub struct FoamFace {
    /// Indices into the `points` array (variable count per face).
    pub node_indices: Vec<usize>,
}

/// One patch (boundary region) in an OpenFOAM polyMesh.
pub struct FoamBoundaryPatch {
    /// Patch name, e.g. `"inlet"`.
    pub name: String,
    /// Patch type, e.g. `"wall"`, `"inlet"`, `"outlet"`, `"symmetryPlane"`.
    pub patch_type: String,
    /// Number of boundary faces in this patch.
    pub n_faces: usize,
    /// Index of the first boundary face in the global face list.
    pub start_face: usize,
}

// ---------------------------------------------------------------------------
// FoamPolyMesh
// ---------------------------------------------------------------------------

/// An OpenFOAM polyMesh: points, faces, owner/neighbour connectivity, boundary.
pub struct FoamPolyMesh {
    /// All mesh points.
    pub points: Vec<FoamPoint>,
    /// All faces (internal + boundary).
    pub faces: Vec<FoamFace>,
    /// Owner cell index for every face.
    pub owner: Vec<usize>,
    /// Neighbour cell index for internal faces only.
    pub neighbour: Vec<usize>,
    /// Boundary patch descriptors.
    pub boundary: Vec<FoamBoundaryPatch>,
}

impl FoamPolyMesh {
    /// Number of cells = max owner index + 1.
    pub fn n_cells(&self) -> usize {
        self.owner.iter().copied().max().map(|m| m + 1).unwrap_or(0)
    }

    /// Number of internal faces = length of `neighbour`.
    pub fn n_internal_faces(&self) -> usize {
        self.neighbour.len()
    }

    /// Number of boundary faces = total faces − internal faces.
    pub fn n_boundary_faces(&self) -> usize {
        self.faces.len().saturating_sub(self.n_internal_faces())
    }

    /// Centroid of each face (average of its node positions).
    pub fn face_centers(&self) -> Vec<[f64; 3]> {
        self.faces
            .iter()
            .map(|f| {
                let n = f.node_indices.len() as f64;
                let mut c = [0.0_f64; 3];
                for &idx in &f.node_indices {
                    let p = self.points[idx];
                    c[0] += p[0];
                    c[1] += p[1];
                    c[2] += p[2];
                }
                [c[0] / n, c[1] / n, c[2] / n]
            })
            .collect()
    }

    /// Approximate area of each face (magnitude of cross-product for triangles/quads).
    pub fn face_areas(&self) -> Vec<f64> {
        self.faces
            .iter()
            .map(|f| face_area(&f.node_indices, &self.points))
            .collect()
    }

    /// Approximate cell centres as the average of their face centres.
    pub fn cell_centers(&self) -> Vec<[f64; 3]> {
        let fc = self.face_centers();
        let nc = self.n_cells();
        let mut sums = vec![[0.0_f64; 3]; nc];
        let mut counts = vec![0usize; nc];

        for (fi, &ow) in self.owner.iter().enumerate() {
            if ow < nc {
                let c = fc[fi];
                sums[ow][0] += c[0];
                sums[ow][1] += c[1];
                sums[ow][2] += c[2];
                counts[ow] += 1;
            }
        }
        for (fi, &nb) in self.neighbour.iter().enumerate() {
            if nb < nc {
                let c = fc[fi];
                sums[nb][0] += c[0];
                sums[nb][1] += c[1];
                sums[nb][2] += c[2];
                counts[nb] += 1;
            }
        }

        sums.iter()
            .zip(counts.iter())
            .map(|(s, &cnt)| {
                if cnt == 0 {
                    [0.0; 3]
                } else {
                    let n = cnt as f64;
                    [s[0] / n, s[1] / n, s[2] / n]
                }
            })
            .collect()
    }

    /// Approximate cell volumes (sum of (face area × distance face-centre to cell-centre) / 3).
    pub fn cell_volumes(&self) -> Vec<f64> {
        let fc = self.face_centers();
        let fa = self.face_areas();
        let cc = self.cell_centers();
        let nc = self.n_cells();
        let mut vols = vec![0.0_f64; nc];

        let contribute = |vols: &mut Vec<f64>, cell: usize, fi: usize| {
            let c = cc[cell];
            let f = fc[fi];
            let d = ((f[0] - c[0]).powi(2) + (f[1] - c[1]).powi(2) + (f[2] - c[2]).powi(2)).sqrt();
            vols[cell] += fa[fi] * d / 3.0;
        };

        for (fi, &ow) in self.owner.iter().enumerate() {
            if ow < nc {
                contribute(&mut vols, ow, fi);
            }
        }
        for (fi, &nb) in self.neighbour.iter().enumerate() {
            if nb < nc {
                contribute(&mut vols, nb, fi);
            }
        }
        vols
    }
}

// ---------------------------------------------------------------------------
// Internal geometry helpers
// ---------------------------------------------------------------------------

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn face_area(indices: &[usize], points: &[FoamPoint]) -> f64 {
    if indices.len() < 3 {
        return 0.0;
    }
    // Fan-triangulation from the first vertex.
    let p0 = points[indices[0]];
    let mut area = 0.0_f64;
    for i in 1..indices.len() - 1 {
        let p1 = points[indices[i]];
        let p2 = points[indices[i + 1]];
        let ab = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let ac = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        area += norm(cross(ab, ac)) * 0.5;
    }
    area
}

// ---------------------------------------------------------------------------
// Parsing helpers
// ---------------------------------------------------------------------------

/// Extract the data block that follows the `FoamFile` header.
/// Returns the content between the first `(` and its matching `)`.
fn extract_data_block(content: &str) -> Option<&str> {
    // Skip FoamFile block if present.
    let search_from = if let Some(pos) = content.find("FoamFile") {
        // find the closing brace
        let brace_open = content[pos..].find('{')? + pos;
        content[brace_open..].find('}')? + brace_open + 1
    } else {
        0
    };
    let rest = &content[search_from..];
    let start = rest.find('(')?;
    let end = rest.rfind(')')?;
    if end > start {
        Some(&rest[start..=end])
    } else {
        None
    }
}

/// Parse an OpenFOAM list of `usize` values.
///
/// Accepts the standard format:
/// ```text
/// N
/// (
///   v1
///   v2
///   …
/// )
/// ```
pub fn parse_foam_list_usize(s: &str) -> Result<Vec<usize>, String> {
    let block =
        extract_data_block(s).ok_or_else(|| "Could not find list block '(' … ')'".to_string())?;
    let inner = &block[1..block.len() - 1]; // strip outer parens
    let mut result = Vec::new();
    for token in inner.split_whitespace() {
        let v: usize = token
            .parse()
            .map_err(|_| format!("Cannot parse usize: {:?}", token))?;
        result.push(v);
    }
    Ok(result)
}

/// Parse an OpenFOAM list of `f64` values.
pub fn parse_foam_list_f64(s: &str) -> Result<Vec<f64>, String> {
    let block =
        extract_data_block(s).ok_or_else(|| "Could not find list block '(' … ')'".to_string())?;
    let inner = &block[1..block.len() - 1];
    let mut result = Vec::new();
    for token in inner.split_whitespace() {
        let v: f64 = token
            .parse()
            .map_err(|_| format!("Cannot parse f64: {:?}", token))?;
        result.push(v);
    }
    Ok(result)
}

/// Parse an OpenFOAM `points` file.
///
/// Format inside the data block:
/// ```text
/// (x y z)
/// (x y z)
/// …
/// ```
pub fn parse_foam_points(content: &str) -> Result<Vec<FoamPoint>, String> {
    let block =
        extract_data_block(content).ok_or_else(|| "Could not find points block".to_string())?;

    let mut points = Vec::new();
    // Each point looks like   (x y z)
    let inner = &block[1..block.len() - 1];
    let chars = inner.char_indices().peekable();
    for (i, ch) in chars {
        if ch == '(' {
            // find matching ')'
            let sub = &inner[i..];
            let end = sub
                .find(')')
                .ok_or_else(|| "Unmatched '(' in points".to_string())?;
            let triplet = &sub[1..end];
            let vals: Vec<f64> = triplet
                .split_whitespace()
                .map(|t| t.parse::<f64>().map_err(|_| format!("Bad f64: {:?}", t)))
                .collect::<Result<_, _>>()?;
            if vals.len() < 3 {
                return Err(format!("Point needs 3 coords, got {:?}", vals));
            }
            points.push([vals[0], vals[1], vals[2]]);
        }
    }
    Ok(points)
}

/// Parse an OpenFOAM `faces` file.
///
/// Each face is encoded as `N(i0 i1 … iN-1)`.
pub fn parse_foam_faces(content: &str) -> Result<Vec<FoamFace>, String> {
    let block =
        extract_data_block(content).ok_or_else(|| "Could not find faces block".to_string())?;
    let inner = &block[1..block.len() - 1];

    let mut faces = Vec::new();
    let mut i = 0;
    let bytes = inner.as_bytes();
    while i < bytes.len() {
        // Skip whitespace
        if bytes[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        // Collect digits (face node count)
        if bytes[i].is_ascii_digit() {
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let _count: usize = inner[start..i]
                .parse()
                .map_err(|_| "Bad face count".to_string())?;
            // Expect '('
            while i < bytes.len() && bytes[i] != b'(' {
                i += 1;
            }
            if i >= bytes.len() {
                return Err("Missing '(' after face count".to_string());
            }
            i += 1; // skip '('
            // Collect until ')'
            let j = i;
            while i < bytes.len() && bytes[i] != b')' {
                i += 1;
            }
            let indices: Vec<usize> = inner[j..i]
                .split_whitespace()
                .map(|t| {
                    t.parse::<usize>()
                        .map_err(|_| format!("Bad index: {:?}", t))
                })
                .collect::<Result<_, _>>()?;
            faces.push(FoamFace {
                node_indices: indices,
            });
            i += 1; // skip ')'
        } else {
            i += 1;
        }
    }
    Ok(faces)
}

/// Parse an OpenFOAM `boundary` file.
pub fn parse_foam_boundary(content: &str) -> Result<Vec<FoamBoundaryPatch>, String> {
    // Skip FoamFile header.
    let search_from = if let Some(pos) = content.find("FoamFile") {
        let brace_open = content[pos..]
            .find('{')
            .ok_or_else(|| "No '{' in FoamFile".to_string())?
            + pos;
        content[brace_open..]
            .find('}')
            .ok_or_else(|| "No '}' in FoamFile".to_string())?
            + brace_open
            + 1
    } else {
        0
    };
    let rest = &content[search_from..];

    // Skip leading count and outer parentheses/braces.
    // boundary format uses { } not ( )
    let outer_open = rest
        .find('(')
        .unwrap_or_else(|| rest.find('{').unwrap_or(0));
    let body = &rest[outer_open..];

    let mut patches = Vec::new();
    // Pattern: name\n{\n  type X;\n  nFaces N;\n  startFace S;\n}
    let mut chars = body.char_indices().peekable();
    while let Some((i, ch)) = chars.next() {
        if ch == '{' || ch == '(' || ch == ')' || ch == '}' {
            continue;
        }
        if ch.is_alphabetic() || ch == '_' {
            // Collect patch name
            let sub = &body[i..];
            let name_end = sub.find(|c: char| c.is_whitespace()).unwrap_or(sub.len());
            let name = sub[..name_end].to_string();
            if name == "FoamFile" {
                continue;
            }
            // Find opening brace for this patch
            let brace_start = match sub.find('{') {
                Some(p) => p,
                None => continue,
            };
            let brace_body_start = brace_start + 1;
            let brace_end = match sub[brace_body_start..].find('}') {
                Some(p) => p + brace_body_start,
                None => continue,
            };
            let patch_body = &sub[brace_body_start..brace_end];

            let mut patch_type = String::new();
            let mut n_faces = 0usize;
            let mut start_face = 0usize;

            for line in patch_body.lines() {
                let line = line.trim().trim_end_matches(';');
                let mut parts = line.splitn(2, char::is_whitespace);
                match (parts.next(), parts.next()) {
                    (Some("type"), Some(v)) => patch_type = v.trim().to_string(),
                    (Some("nFaces"), Some(v)) => n_faces = v.trim().parse().unwrap_or(0),
                    (Some("startFace"), Some(v)) => start_face = v.trim().parse().unwrap_or(0),
                    _ => {}
                }
            }

            if !name.is_empty() && !patch_type.is_empty() {
                patches.push(FoamBoundaryPatch {
                    name,
                    patch_type,
                    n_faces,
                    start_face,
                });
            }

            // Advance the outer iterator past the consumed region.
            let consumed = i + brace_end + 1;
            // Skip chars until we've advanced past 'consumed'.
            while chars.peek().map(|&(j, _)| j < consumed).unwrap_or(false) {
                chars.next();
            }
        }
    }
    Ok(patches)
}

// ---------------------------------------------------------------------------
// Writing functions
// ---------------------------------------------------------------------------

/// Generate a minimal FoamFile header block.
pub fn write_foam_header(class: &str, object: &str) -> String {
    FoamHeader {
        foam_file_version: "2.0".to_string(),
        format: "ascii".to_string(),
        class_name: class.to_string(),
        object_name: object.to_string(),
    }
    .to_string()
}

/// Serialise a points list to OpenFOAM ASCII format.
pub fn write_foam_points(points: &[FoamPoint]) -> String {
    let mut s = write_foam_header("vectorField", "points");
    s.push('\n');
    s.push_str(&format!("{}\n(\n", points.len()));
    for p in points {
        s.push_str(&format!("    ({} {} {})\n", p[0], p[1], p[2]));
    }
    s.push_str(")\n");
    s
}

/// Serialise a faces list to OpenFOAM ASCII format.
pub fn write_foam_faces(faces: &[FoamFace]) -> String {
    let mut s = write_foam_header("faceList", "faces");
    s.push('\n');
    s.push_str(&format!("{}\n(\n", faces.len()));
    for f in faces {
        let indices: Vec<String> = f.node_indices.iter().map(|i| i.to_string()).collect();
        s.push_str(&format!(
            "    {}({})\n",
            f.node_indices.len(),
            indices.join(" ")
        ));
    }
    s.push_str(")\n");
    s
}

/// Serialise owner and neighbour lists to OpenFOAM ASCII format.
/// Returns `(owner_string, neighbour_string)`.
pub fn write_foam_owner_neighbour(owner: &[usize], neighbour: &[usize]) -> (String, String) {
    let mut o = write_foam_header("labelList", "owner");
    o.push('\n');
    o.push_str(&format!("{}\n(\n", owner.len()));
    for &v in owner {
        o.push_str(&format!("    {}\n", v));
    }
    o.push_str(")\n");

    let mut n = write_foam_header("labelList", "neighbour");
    n.push('\n');
    n.push_str(&format!("{}\n(\n", neighbour.len()));
    for &v in neighbour {
        n.push_str(&format!("    {}\n", v));
    }
    n.push_str(")\n");

    (o, n)
}

/// Write a `volScalarField` file (e.g. pressure `p` or temperature `T`).
pub fn write_scalar_field(name: &str, values: &[f64], boundary: &[FoamBoundaryPatch]) -> String {
    let mut s = write_foam_header("volScalarField", name);
    s.push('\n');
    s.push_str("dimensions      [0 0 0 0 0 0 0];\n\n");
    s.push_str(&format!(
        "internalField   nonuniform List<scalar>\n{}\n(\n",
        values.len()
    ));
    for &v in values {
        s.push_str(&format!("    {}\n", v));
    }
    s.push_str(");\n\nboundaryField\n{\n");
    for patch in boundary {
        s.push_str(&format!(
            "    {}\n    {{\n        type            zeroGradient;\n    }}\n",
            patch.name
        ));
    }
    s.push_str("}\n");
    s
}

// ---------------------------------------------------------------------------
// FoamScalarField
// ---------------------------------------------------------------------------

/// A named scalar field on the mesh (internal + boundary).
pub struct FoamScalarField {
    /// Field name, e.g. `"p"` or `"T"`.
    pub name: String,
    /// Internal-cell values (one per cell).
    pub internal_field: Vec<f64>,
    /// Boundary values (one `Vec`f64` per patch; may be empty for zeroGradient).
    pub boundary_values: Vec<Vec<f64>>,
}

impl FoamScalarField {
    /// Create a uniform field with all cells set to `value`.
    pub fn uniform(name: &str, n_cells: usize, value: f64) -> Self {
        FoamScalarField {
            name: name.to_string(),
            internal_field: vec![value; n_cells],
            boundary_values: Vec::new(),
        }
    }

    /// Serialise to an OpenFOAM `volScalarField` string.
    pub fn to_foam_string(&self) -> String {
        let mut s = write_foam_header("volScalarField", &self.name);
        s.push('\n');
        s.push_str("dimensions      [0 0 0 0 0 0 0];\n\n");
        if self.internal_field.is_empty() {
            s.push_str("internalField   uniform 0;\n");
        } else {
            s.push_str(&format!(
                "internalField   nonuniform List<scalar>\n{}\n(\n",
                self.internal_field.len()
            ));
            for &v in &self.internal_field {
                s.push_str(&format!("    {}\n", v));
            }
            s.push_str(");\n");
        }
        s.push_str("\nboundaryField\n{\n");
        for (i, bv) in self.boundary_values.iter().enumerate() {
            let patch_name = format!("patch{}", i);
            if bv.is_empty() {
                s.push_str(&format!(
                    "    {}\n    {{\n        type            zeroGradient;\n    }}\n",
                    patch_name
                ));
            } else {
                s.push_str(&format!(
                    "    {}\n    {{\n        type            fixedValue;\n        value           nonuniform List<scalar>\n{}\n(\n",
                    patch_name, bv.len()
                ));
                for &v in bv {
                    s.push_str(&format!("            {}\n", v));
                }
                s.push_str("        );\n    }\n");
            }
        }
        s.push_str("}\n");
        s
    }
}

// ---------------------------------------------------------------------------
// SimpleCaseWriter
// ---------------------------------------------------------------------------

/// Writes a minimal OpenFOAM case directory structure in memory.
pub struct SimpleCaseWriter {
    /// Path to the case directory (used for `file_list` output).
    pub case_dir: String,
    /// The polyMesh to write.
    pub mesh: FoamPolyMesh,
}

impl SimpleCaseWriter {
    /// Create a new writer for `case_dir` with the given mesh.
    pub fn new(case_dir: &str, mesh: FoamPolyMesh) -> Self {
        SimpleCaseWriter {
            case_dir: case_dir.to_string(),
            mesh,
        }
    }

    /// Return all mesh file contents concatenated (useful for testing without disk I/O).
    pub fn write_mesh_to_string(&self) -> String {
        let mut out = String::new();
        out.push_str("=== constant/polyMesh/points ===\n");
        out.push_str(&write_foam_points(&self.mesh.points));
        out.push_str("\n=== constant/polyMesh/faces ===\n");
        out.push_str(&write_foam_faces(&self.mesh.faces));
        let (owner_str, neighbour_str) =
            write_foam_owner_neighbour(&self.mesh.owner, &self.mesh.neighbour);
        out.push_str("\n=== constant/polyMesh/owner ===\n");
        out.push_str(&owner_str);
        out.push_str("\n=== constant/polyMesh/neighbour ===\n");
        out.push_str(&neighbour_str);
        out
    }

    /// List the files that would be written to disk.
    pub fn file_list(&self) -> Vec<String> {
        let base = &self.case_dir;
        vec![
            format!("{}/constant/polyMesh/points", base),
            format!("{}/constant/polyMesh/faces", base),
            format!("{}/constant/polyMesh/owner", base),
            format!("{}/constant/polyMesh/neighbour", base),
            format!("{}/constant/polyMesh/boundary", base),
            format!("{}/system/controlDict", base),
            format!("{}/system/fvSchemes", base),
            format!("{}/system/fvSolution", base),
        ]
    }
}

// ---------------------------------------------------------------------------
// OpenFOAM field reading
// ---------------------------------------------------------------------------

/// Parse an OpenFOAM `volVectorField` internal field (list of 3-vectors).
///
/// Accepts blocks like:
/// ```text
/// internalField   nonuniform List`vector`
/// N
/// (
///   (x y z)
///   …
/// )
/// ```
pub fn parse_foam_vector_field(content: &str) -> Result<Vec<[f64; 3]>, String> {
    parse_foam_points(content)
}

/// Parse an OpenFOAM `volScalarField` internal field from its text representation.
///
/// Handles both `uniform VALUE` and `nonuniform List`scalar` forms.
pub fn parse_foam_scalar_field_internal(content: &str) -> Result<Vec<f64>, String> {
    // Check for "uniform VALUE".
    for line in content.lines() {
        let line = line.trim();
        if let Some(after_field) = line.strip_prefix("internalField") {
            let rest = &after_field.trim_start();
            if let Some(rest2) = rest.strip_prefix("uniform") {
                let v: f64 = rest2
                    .trim()
                    .trim_end_matches(';')
                    .parse()
                    .map_err(|_| format!("Cannot parse uniform scalar: {:?}", rest2))?;
                return Ok(vec![v]);
            }
        }
    }
    // Fall back to list parsing.
    parse_foam_list_f64(content)
}

// ---------------------------------------------------------------------------
// OpenFOAM mesh export helpers
// ---------------------------------------------------------------------------

/// Generate the `constant/polyMesh/boundary` file content for a mesh.
pub fn write_foam_boundary(patches: &[FoamBoundaryPatch]) -> String {
    let mut s = write_foam_header("polyBoundaryMesh", "boundary");
    s.push('\n');
    s.push_str(&format!("{}\n(\n", patches.len()));
    for p in patches {
        s.push_str(&format!(
            "    {}\n    {{\n        type            {};\n        nFaces          {};\n        startFace       {};\n    }}\n",
            p.name, p.patch_type, p.n_faces, p.start_face
        ));
    }
    s.push_str(")\n");
    s
}

/// Generate a minimal `system/controlDict` content for an OpenFOAM case.
pub fn write_foam_control_dict(
    start_time: f64,
    end_time: f64,
    delta_t: f64,
    write_interval: f64,
) -> String {
    let mut s = write_foam_header("dictionary", "controlDict");
    s.push('\n');
    s.push_str("application     icoFoam;\n\n");
    s.push_str(&format!(
        "startFrom       latestTime;\nstartTime       {};\n",
        start_time
    ));
    s.push_str(&format!(
        "stopAt          endTime;\nendTime         {};\n",
        end_time
    ));
    s.push_str(&format!("deltaT          {};\n", delta_t));
    s.push_str(&format!(
        "writeControl    timeStep;\nwriteInterval   {};\n",
        write_interval
    ));
    s.push_str("purgeWrite      0;\nwriteFormat     ascii;\nwritePrecision  6;\n");
    s
}

/// Generate a minimal `system/fvSchemes` dictionary.
pub fn write_foam_fv_schemes() -> String {
    let mut s = write_foam_header("dictionary", "fvSchemes");
    s.push('\n');
    s.push_str("ddtSchemes\n{\n    default         Euler;\n}\n\n");
    s.push_str("gradSchemes\n{\n    default         Gauss linear;\n}\n\n");
    s.push_str("divSchemes\n{\n    default         none;\n    div(phi,U)      Gauss linearUpwind grad(U);\n}\n\n");
    s.push_str("laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n");
    s.push_str("interpolationSchemes\n{\n    default         linear;\n}\n\n");
    s.push_str("snGradSchemes\n{\n    default         corrected;\n}\n");
    s
}

/// Generate a minimal `system/fvSolution` dictionary.
pub fn write_foam_fv_solution() -> String {
    let mut s = write_foam_header("dictionary", "fvSolution");
    s.push('\n');
    s.push_str("solvers\n{\n");
    s.push_str("    p\n    {\n        solver          PCG;\n        preconditioner  DIC;\n        tolerance       1e-06;\n        relTol          0.05;\n    }\n\n");
    s.push_str("    U\n    {\n        solver          smoothSolver;\n        smoother        symGaussSeidel;\n        tolerance       1e-05;\n        relTol          0;\n    }\n}\n\n");
    s.push_str("PISO\n{\n    nCorrectors     2;\n    nNonOrthogonalCorrectors 0;\n}\n");
    s
}

// ---------------------------------------------------------------------------
// Boundary condition writing
// ---------------------------------------------------------------------------

/// Standard OpenFOAM boundary condition types.
#[derive(Debug, Clone, PartialEq)]
pub enum FoamBcType {
    /// Zero-gradient (Neumann) boundary condition.
    ZeroGradient,
    /// Fixed scalar value (Dirichlet) boundary condition.
    FixedValue(f64),
    /// Fixed vector value boundary condition.
    FixedVectorValue([f64; 3]),
    /// Inlet-outlet condition with a given reference value.
    InletOutlet(f64),
    /// No-slip wall condition for velocity.
    NoSlip,
    /// Slip wall condition.
    Slip,
    /// Empty (2-D / unused patch) condition.
    Empty,
}

impl FoamBcType {
    /// Serialise to the OpenFOAM dictionary sub-block (without patch name).
    pub fn to_foam_block(&self) -> String {
        match self {
            FoamBcType::ZeroGradient => "        type            zeroGradient;\n".to_string(),
            FoamBcType::FixedValue(v) => format!(
                "        type            fixedValue;\n        value           uniform {};\n",
                v
            ),
            FoamBcType::FixedVectorValue(v) => format!(
                "        type            fixedValue;\n        value           uniform ({} {} {});\n",
                v[0], v[1], v[2]
            ),
            FoamBcType::InletOutlet(v) => format!(
                "        type            inletOutlet;\n        inletValue      uniform {};\n        value           uniform {};\n",
                v, v
            ),
            FoamBcType::NoSlip => "        type            noSlip;\n".to_string(),
            FoamBcType::Slip => "        type            slip;\n".to_string(),
            FoamBcType::Empty => "        type            empty;\n".to_string(),
        }
    }
}

/// Write an OpenFOAM scalar field file with per-patch boundary conditions.
pub fn write_scalar_field_with_bc(
    name: &str,
    dimensions: &str,
    internal_value: &str,
    patches: &[(&str, FoamBcType)],
) -> String {
    let mut s = write_foam_header("volScalarField", name);
    s.push('\n');
    s.push_str(&format!("dimensions      {};\n\n", dimensions));
    s.push_str(&format!(
        "internalField   {};\n\nboundaryField\n{{\n",
        internal_value
    ));
    for (patch_name, bc) in patches {
        s.push_str(&format!(
            "    {}\n    {{\n{}",
            patch_name,
            bc.to_foam_block()
        ));
        s.push_str("    }\n");
    }
    s.push_str("}\n");
    s
}

/// Write an OpenFOAM vector field file.
pub fn write_vector_field(
    name: &str,
    dimensions: &str,
    values: &[[f64; 3]],
    patches: &[(&str, FoamBcType)],
) -> String {
    let mut s = write_foam_header("volVectorField", name);
    s.push('\n');
    s.push_str(&format!("dimensions      {};\n\n", dimensions));
    if values.is_empty() {
        s.push_str("internalField   uniform (0 0 0);\n");
    } else {
        s.push_str(&format!(
            "internalField   nonuniform List<vector>\n{}\n(\n",
            values.len()
        ));
        for v in values {
            s.push_str(&format!("    ({} {} {})\n", v[0], v[1], v[2]));
        }
        s.push_str(");\n");
    }
    s.push_str("\nboundaryField\n{\n");
    for (patch_name, bc) in patches {
        s.push_str(&format!(
            "    {}\n    {{\n{}",
            patch_name,
            bc.to_foam_block()
        ));
        s.push_str("    }\n");
    }
    s.push_str("}\n");
    s
}

// ---------------------------------------------------------------------------
// Parallel foam decomposition
// ---------------------------------------------------------------------------

/// Decompose a scalar field into `n_procs` subdomains.
///
/// Uses a simple contiguous partitioning: processor `p` owns cells
/// `[p * chunk_size, (p+1) * chunk_size)`.  Returns one `Vec`f64` per process.
pub fn decompose_scalar_field(values: &[f64], n_procs: usize) -> Vec<Vec<f64>> {
    if n_procs == 0 || values.is_empty() {
        return vec![values.to_vec()];
    }
    let n = values.len();
    let chunk = n.div_ceil(n_procs);
    (0..n_procs)
        .map(|p| {
            let start = (p * chunk).min(n);
            let end = ((p + 1) * chunk).min(n);
            values[start..end].to_vec()
        })
        .collect()
}

/// Reconstruct a scalar field from decomposed sub-domain slices.
pub fn reconstruct_scalar_field(parts: &[Vec<f64>]) -> Vec<f64> {
    parts.iter().flat_map(|v| v.iter().copied()).collect()
}

/// Decompose cell ownership for parallel processing.
///
/// Returns `n_procs` Vec`usize` each containing the global cell indices
/// that belong to processor `p`.
pub fn decompose_cell_indices(n_cells: usize, n_procs: usize) -> Vec<Vec<usize>> {
    if n_procs == 0 {
        return vec![(0..n_cells).collect()];
    }
    let chunk = n_cells.div_ceil(n_procs);
    (0..n_procs)
        .map(|p| {
            let start = (p * chunk).min(n_cells);
            let end = ((p + 1) * chunk).min(n_cells);
            (start..end).collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Foam case structure
// ---------------------------------------------------------------------------

/// Full in-memory representation of an OpenFOAM case.
pub struct FoamCase {
    /// Case directory path.
    pub case_dir: String,
    /// The polyMesh.
    pub mesh: FoamPolyMesh,
    /// Named scalar fields at time 0.
    pub scalar_fields: Vec<FoamScalarField>,
    /// Simulation start time (s).
    pub start_time: f64,
    /// Simulation end time (s).
    pub end_time: f64,
    /// Time step size (s).
    pub delta_t: f64,
    /// Write interval (in time steps or seconds, depending on writeControl).
    pub write_interval: f64,
}

impl FoamCase {
    /// Create a new case.
    pub fn new(case_dir: &str, mesh: FoamPolyMesh) -> Self {
        Self {
            case_dir: case_dir.to_string(),
            mesh,
            scalar_fields: Vec::new(),
            start_time: 0.0,
            end_time: 1.0,
            delta_t: 0.001,
            write_interval: 0.1,
        }
    }

    /// Add a scalar field to the initial conditions.
    pub fn add_scalar_field(&mut self, field: FoamScalarField) {
        self.scalar_fields.push(field);
    }

    /// Return a list of all files this case would write to disk.
    pub fn file_list(&self) -> Vec<String> {
        let base = &self.case_dir;
        let mut files = vec![
            format!("{}/constant/polyMesh/points", base),
            format!("{}/constant/polyMesh/faces", base),
            format!("{}/constant/polyMesh/owner", base),
            format!("{}/constant/polyMesh/neighbour", base),
            format!("{}/constant/polyMesh/boundary", base),
            format!("{}/system/controlDict", base),
            format!("{}/system/fvSchemes", base),
            format!("{}/system/fvSolution", base),
        ];
        for sf in &self.scalar_fields {
            files.push(format!("{}/0/{}", base, sf.name));
        }
        files
    }

    /// Generate the controlDict content.
    pub fn control_dict_content(&self) -> String {
        write_foam_control_dict(
            self.start_time,
            self.end_time,
            self.delta_t,
            self.write_interval,
        )
    }

    /// Generate the boundary file content.
    pub fn boundary_content(&self) -> String {
        write_foam_boundary(&self.mesh.boundary)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_mesh() -> FoamPolyMesh {
        // Two tetrahedra sharing one internal face.
        // 5 points, 7 faces (1 internal + 6 boundary), 2 cells.
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
            [1.5, 0.5, 1.0],
        ];
        let faces = vec![
            FoamFace {
                node_indices: vec![0, 1, 3],
            }, // internal face
            FoamFace {
                node_indices: vec![0, 2, 3],
            },
            FoamFace {
                node_indices: vec![1, 2, 3],
            },
            FoamFace {
                node_indices: vec![0, 1, 2],
            },
            FoamFace {
                node_indices: vec![1, 4, 3],
            },
            FoamFace {
                node_indices: vec![1, 2, 4],
            },
            FoamFace {
                node_indices: vec![3, 4, 2],
            },
        ];
        let owner = vec![0, 0, 0, 0, 1, 1, 1];
        let neighbour = vec![1]; // face 0 is shared
        let boundary = vec![FoamBoundaryPatch {
            name: "wall".to_string(),
            patch_type: "wall".to_string(),
            n_faces: 6,
            start_face: 1,
        }];
        FoamPolyMesh {
            points,
            faces,
            owner,
            neighbour,
            boundary,
        }
    }

    #[test]
    fn test_foam_header_to_string_contains_foamfile() {
        let h = FoamHeader {
            foam_file_version: "2.0".to_string(),
            format: "ascii".to_string(),
            class_name: "polyBoundaryMesh".to_string(),
            object_name: "boundary".to_string(),
        };
        let s = h.to_string();
        assert!(s.contains("FoamFile"));
        assert!(s.contains("polyBoundaryMesh"));
    }

    #[test]
    fn test_parse_foam_list_usize() {
        let input = "3\n(\n1\n2\n3\n)";
        let result = parse_foam_list_usize(input).expect("parse failed");
        assert_eq!(result, vec![1, 2, 3]);
    }

    #[test]
    fn test_parse_foam_points() {
        let input = "1\n(\n    (3.0 2.0 1.0)\n)";
        let pts = parse_foam_points(input).expect("parse failed");
        assert_eq!(pts.len(), 1);
        assert!((pts[0][0] - 3.0).abs() < 1e-12);
        assert!((pts[0][1] - 2.0).abs() < 1e-12);
        assert!((pts[0][2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_n_cells() {
        let mesh = simple_mesh();
        assert_eq!(mesh.n_cells(), 2);
    }

    #[test]
    fn test_write_foam_header_contains_class() {
        let s = write_foam_header("volScalarField", "p");
        assert!(s.contains("FoamFile"));
        assert!(s.contains("volScalarField"));
        assert!(s.contains("p"));
    }

    #[test]
    fn test_foam_scalar_field_uniform() {
        let f = FoamScalarField::uniform("T", 4, 300.0);
        assert_eq!(f.internal_field.len(), 4);
        assert!(f.internal_field.iter().all(|&v| (v - 300.0).abs() < 1e-12));
    }

    #[test]
    fn test_face_centers_centroid() {
        let mesh = simple_mesh();
        let fc = mesh.face_centers();
        // Face 0: nodes 0,1,3 → ([0,0,0],[1,0,0],[0.5,0.5,1.0]) → centre=(0.5, 1/6, 1/3)
        // Actually: (0+1+0.5)/3, (0+0+0.5)/3, (0+0+1)/3
        assert!((fc[0][0] - (0.0 + 1.0 + 0.5) / 3.0).abs() < 1e-9);
        assert!((fc[0][1] - (0.0 + 0.0 + 0.5) / 3.0).abs() < 1e-9);
        assert!((fc[0][2] - (0.0 + 0.0 + 1.0) / 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_simple_case_writer_file_list_nonempty() {
        let tmpdir = std::env::temp_dir();
        let writer = SimpleCaseWriter::new(
            tmpdir.join("test_case").to_str().unwrap_or(""),
            simple_mesh(),
        );
        let files = writer.file_list();
        assert!(!files.is_empty());
        assert!(files.iter().any(|f| f.contains("points")));
    }

    // -----------------------------------------------------------------------
    // Boundary writing tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_write_foam_boundary_contains_patch_name() {
        let mesh = simple_mesh();
        let s = write_foam_boundary(&mesh.boundary);
        assert!(
            s.contains("wall"),
            "Boundary output should contain patch name 'wall'"
        );
        assert!(
            s.contains("FoamFile"),
            "Boundary output should contain FoamFile header"
        );
    }

    #[test]
    fn test_write_foam_control_dict() {
        let s = write_foam_control_dict(0.0, 1.0, 0.001, 0.1);
        assert!(s.contains("endTime"), "controlDict should contain endTime");
        assert!(s.contains("deltaT"), "controlDict should contain deltaT");
    }

    #[test]
    fn test_write_foam_fv_schemes() {
        let s = write_foam_fv_schemes();
        assert!(
            s.contains("ddtSchemes"),
            "fvSchemes should contain ddtSchemes"
        );
        assert!(
            s.contains("gradSchemes"),
            "fvSchemes should contain gradSchemes"
        );
    }

    #[test]
    fn test_write_foam_fv_solution() {
        let s = write_foam_fv_solution();
        assert!(s.contains("solvers"), "fvSolution should contain solvers");
        assert!(s.contains("PISO"), "fvSolution should contain PISO block");
    }

    // -----------------------------------------------------------------------
    // Boundary condition tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_foam_bc_zero_gradient() {
        let bc = FoamBcType::ZeroGradient;
        let s = bc.to_foam_block();
        assert!(s.contains("zeroGradient"));
    }

    #[test]
    fn test_foam_bc_fixed_value() {
        let bc = FoamBcType::FixedValue(42.0);
        let s = bc.to_foam_block();
        assert!(s.contains("fixedValue"));
        assert!(s.contains("42"));
    }

    #[test]
    fn test_foam_bc_no_slip() {
        let bc = FoamBcType::NoSlip;
        let s = bc.to_foam_block();
        assert!(s.contains("noSlip"));
    }

    #[test]
    fn test_write_scalar_field_with_bc() {
        let patches = vec![
            ("inlet", FoamBcType::FixedValue(1.0)),
            ("outlet", FoamBcType::ZeroGradient),
        ];
        let s = write_scalar_field_with_bc("p", "[0 2 -2 0 0 0 0]", "uniform 0", &patches);
        assert!(s.contains("inlet"));
        assert!(s.contains("outlet"));
        assert!(s.contains("fixedValue"));
        assert!(s.contains("zeroGradient"));
    }

    #[test]
    fn test_write_vector_field() {
        let vals = vec![[1.0, 0.0, 0.0], [0.5, 0.5, 0.0]];
        let patches: Vec<(&str, FoamBcType)> = vec![("walls", FoamBcType::NoSlip)];
        let s = write_vector_field("U", "[0 1 -1 0 0 0 0]", &vals, &patches);
        assert!(s.contains("volVectorField"));
        assert!(s.contains("noSlip"));
    }

    // -----------------------------------------------------------------------
    // Parallel decomposition tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_decompose_scalar_field_total_count() {
        let values: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let parts = decompose_scalar_field(&values, 4);
        let total: usize = parts.iter().map(|p| p.len()).sum();
        assert_eq!(total, 100, "Decomposed parts should cover all cells");
    }

    #[test]
    fn test_reconstruct_scalar_field() {
        let original: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let parts = decompose_scalar_field(&original, 5);
        let reconstructed = reconstruct_scalar_field(&parts);
        assert_eq!(reconstructed, original);
    }

    #[test]
    fn test_decompose_cell_indices() {
        let n_cells = 100;
        let n_procs = 4;
        let parts = decompose_cell_indices(n_cells, n_procs);
        assert_eq!(parts.len(), n_procs);
        let total: usize = parts.iter().map(|p| p.len()).sum();
        assert_eq!(total, n_cells);
        // All indices should be unique and cover 0..n_cells.
        let mut all: Vec<usize> = parts.into_iter().flatten().collect();
        all.sort_unstable();
        let expected: Vec<usize> = (0..n_cells).collect();
        assert_eq!(all, expected);
    }

    #[test]
    fn test_decompose_single_proc() {
        let values = vec![1.0, 2.0, 3.0];
        let parts = decompose_scalar_field(&values, 1);
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0], values);
    }

    // -----------------------------------------------------------------------
    // FoamCase tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_foam_case_file_list() {
        let mesh = simple_mesh();
        let tmpdir = std::env::temp_dir();
        let case = FoamCase::new(tmpdir.join("mycase").to_str().unwrap_or(""), mesh);
        let files = case.file_list();
        assert!(files.iter().any(|f| f.contains("controlDict")));
        assert!(files.iter().any(|f| f.contains("fvSchemes")));
    }

    #[test]
    fn test_foam_case_with_scalar_field() {
        let mesh = simple_mesh();
        let tmpdir = std::env::temp_dir();
        let mut case = FoamCase::new(tmpdir.join("mycase2").to_str().unwrap_or(""), mesh);
        case.add_scalar_field(FoamScalarField::uniform("p", 10, 0.0));
        let files = case.file_list();
        assert!(files.iter().any(|f| f.ends_with("/0/p")));
    }

    #[test]
    fn test_foam_case_control_dict_content() {
        let mesh = simple_mesh();
        let tmpdir = std::env::temp_dir();
        let mut case = FoamCase::new(tmpdir.join("mycase3").to_str().unwrap_or(""), mesh);
        case.end_time = 2.5;
        let s = case.control_dict_content();
        assert!(s.contains("2.5"), "controlDict should contain the end time");
    }

    #[test]
    fn test_foam_case_boundary_content() {
        let mesh = simple_mesh();
        let tmpdir = std::env::temp_dir();
        let case = FoamCase::new(tmpdir.join("mycase4").to_str().unwrap_or(""), mesh);
        let s = case.boundary_content();
        assert!(
            s.contains("wall"),
            "boundary content should include patch name"
        );
    }

    // -----------------------------------------------------------------------
    // Field parsing tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_foam_vector_field() {
        let input = "1\n(\n    (1.0 2.0 3.0)\n)";
        let vecs = parse_foam_vector_field(input).expect("parse failed");
        assert_eq!(vecs.len(), 1);
        assert!((vecs[0][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_parse_foam_scalar_field_internal_uniform() {
        let content = "internalField   uniform 101325;\n";
        let vals = parse_foam_scalar_field_internal(content).expect("parse failed");
        assert_eq!(vals.len(), 1);
        assert!((vals[0] - 101325.0).abs() < 1e-8);
    }

    #[test]
    fn test_write_foam_faces_round_trip() {
        let faces = vec![
            FoamFace {
                node_indices: vec![0, 1, 2],
            },
            FoamFace {
                node_indices: vec![1, 3, 2],
            },
        ];
        let s = write_foam_faces(&faces);
        // Parse back.
        let parsed = parse_foam_faces(&s).expect("round-trip parse failed");
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].node_indices, vec![0, 1, 2]);
        assert_eq!(parsed[1].node_indices, vec![1, 3, 2]);
    }
}
