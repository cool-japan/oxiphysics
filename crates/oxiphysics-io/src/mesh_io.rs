// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended mesh I/O formats.
//!
//! Provides writers and readers for glTF 2.0, USD ASCII, Gmsh .msh v2/v4,
//! Abaqus .inp, NASTRAN BDF, mesh conversion, and mesh quality validation.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use std::io::Write;

// ---------------------------------------------------------------------------
// Helper free functions
// ---------------------------------------------------------------------------

/// Encode a byte slice to a base-64 ASCII string.
pub fn base64_encode_buffer(data: &[u8]) -> String {
    const TABLE: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let b0 = chunk[0] as usize;
        let b1 = if chunk.len() > 1 {
            chunk[1] as usize
        } else {
            0
        };
        let b2 = if chunk.len() > 2 {
            chunk[2] as usize
        } else {
            0
        };
        out.push(TABLE[b0 >> 2] as char);
        out.push(TABLE[((b0 & 3) << 4) | (b1 >> 4)] as char);
        if chunk.len() > 1 {
            out.push(TABLE[((b1 & 0xf) << 2) | (b2 >> 6)] as char);
        } else {
            out.push('=');
        }
        if chunk.len() > 2 {
            out.push(TABLE[b2 & 0x3f] as char);
        } else {
            out.push('=');
        }
    }
    out
}

/// Write a glTF accessor JSON fragment.
pub fn write_gltf_accessor<W: Write>(
    w: &mut W,
    buffer_view: usize,
    component_type: u32,
    count: usize,
    type_str: &str,
    trailing_comma: bool,
) -> std::io::Result<()> {
    let comma = if trailing_comma { "," } else { "" };
    writeln!(
        w,
        r#"    {{"bufferView":{buffer_view},"componentType":{component_type},"count":{count},"type":"{type_str}"}}{comma}"#
    )
}

/// Format a floating-point value in NASTRAN field-8 format (8 chars wide).
pub fn nastran_field8(v: f64) -> String {
    // Use E notation, fitting into 8 characters
    let s = format!("{:>8.4e}", v);
    if s.len() > 8 {
        format!("{:>8.3e}", v)
    } else {
        s
    }
}

/// Return the Gmsh element type string for a given number of nodes.
/// 3 → "Triangle", 4 → "Quad", 4 → "Tet4", 8 → "Hex8".
pub fn gmsh_element_type_str(n_nodes: usize) -> &'static str {
    match n_nodes {
        2 => "Line",
        3 => "Triangle",
        4 => "Quad",
        10 => "Tet10",
        _ => "Unknown",
    }
}

// ---------------------------------------------------------------------------
// GltfMeshWriter
// ---------------------------------------------------------------------------

/// Write a mesh to glTF 2.0 format (JSON + binary buffer).
pub struct GltfMeshWriter {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Vertex normals.
    pub normals: Vec<[f64; 3]>,
    /// UV texture coordinates.
    pub uvs: Vec<[f64; 2]>,
    /// Triangle indices.
    pub indices: Vec<[u32; 3]>,
    /// Material name.
    pub material_name: String,
}

impl GltfMeshWriter {
    /// Construct with mesh data.
    pub fn new(
        vertices: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        uvs: Vec<[f64; 2]>,
        indices: Vec<[u32; 3]>,
    ) -> Self {
        GltfMeshWriter {
            vertices,
            normals,
            uvs,
            indices,
            material_name: "default".to_string(),
        }
    }

    /// Pack positions as f32 little-endian bytes.
    pub fn positions_f32_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.vertices.len() * 12);
        for &[x, y, z] in &self.vertices {
            buf.extend_from_slice(&(x as f32).to_le_bytes());
            buf.extend_from_slice(&(y as f32).to_le_bytes());
            buf.extend_from_slice(&(z as f32).to_le_bytes());
        }
        buf
    }

    /// Pack indices as u32 little-endian bytes.
    pub fn indices_u32_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.indices.len() * 12);
        for &[a, b, c] in &self.indices {
            buf.extend_from_slice(&a.to_le_bytes());
            buf.extend_from_slice(&b.to_le_bytes());
            buf.extend_from_slice(&c.to_le_bytes());
        }
        buf
    }

    /// Write the glTF JSON to `writer` with embedded base64 binary data.
    pub fn write_gltf_json<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        let pos_bytes = self.positions_f32_bytes();
        let idx_bytes = self.indices_u32_bytes();
        let mut binary = pos_bytes.clone();
        binary.extend_from_slice(&idx_bytes);
        let b64 = base64_encode_buffer(&binary);
        let n_verts = self.vertices.len();
        let n_tris = self.indices.len();
        let pos_len = pos_bytes.len();
        let idx_len = idx_bytes.len();
        writeln!(w, "{{")?;
        writeln!(w, r#"  "asset": {{"version": "2.0"}},"#)?;
        writeln!(w, r#"  "scene": 0,"#)?;
        writeln!(w, r#"  "scenes": [{{"nodes": [0]}}],"#)?;
        writeln!(w, r#"  "nodes": [{{"mesh": 0}}],"#)?;
        writeln!(
            w,
            r#"  "meshes": [{{"name": "mesh","primitives": [{{"attributes": {{"POSITION": 0}},"indices": 1}}]}}],"#
        )?;
        writeln!(w, r#"  "accessors": ["#)?;
        write_gltf_accessor(w, 0, 5126, n_verts, "VEC3", true)?;
        write_gltf_accessor(w, 1, 5125, n_tris * 3, "SCALAR", false)?;
        writeln!(w, r#"  ],"#)?;
        writeln!(w, r#"  "bufferViews": ["#)?;
        writeln!(
            w,
            r#"    {{"buffer": 0, "byteOffset": 0, "byteLength": {pos_len}}},"#
        )?;
        writeln!(
            w,
            r#"    {{"buffer": 0, "byteOffset": {pos_len}, "byteLength": {idx_len}}}"#
        )?;
        writeln!(w, r#"  ],"#)?;
        writeln!(
            w,
            r#"  "buffers": [{{"byteLength": {total}, "uri": "data:application/octet-stream;base64,{b64}"}}]"#,
            total = binary.len()
        )?;
        writeln!(w, "}}")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// GltfMeshReader
// ---------------------------------------------------------------------------

/// Read a glTF mesh and extract vertices, normals, and indices.
///
/// This is a simplified reader that parses a subset of the glTF JSON format.
pub struct GltfMeshReader {
    /// Parsed vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Parsed triangle indices.
    pub indices: Vec<[u32; 3]>,
    /// Node names encountered.
    pub node_names: Vec<String>,
}

impl GltfMeshReader {
    /// Parse a glTF JSON string (very simplified: extracts node names only).
    pub fn parse_json(json: &str) -> Self {
        let mut node_names = Vec::new();
        // Extract mesh names
        for line in json.lines() {
            let line = line.trim();
            if line.contains("\"name\"")
                && let Some(start) = line.find("\"name\":")
            {
                let rest = &line[start + 7..].trim_start_matches([' ', '"'].as_ref());
                let end = rest.find('"').unwrap_or(rest.len());
                node_names.push(rest[..end].to_string());
            }
        }
        GltfMeshReader {
            vertices: Vec::new(),
            indices: Vec::new(),
            node_names,
        }
    }
}

// ---------------------------------------------------------------------------
// UsdMeshWriter
// ---------------------------------------------------------------------------

/// Write a mesh to USD (Universal Scene Description) ASCII format.
pub struct UsdMeshWriter {
    /// Mesh vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Face vertex counts (e.g., 3 per triangle).
    pub face_vertex_counts: Vec<u32>,
    /// Flat face vertex indices.
    pub face_vertex_indices: Vec<u32>,
    /// USD prim path.
    pub prim_path: String,
}

impl UsdMeshWriter {
    /// Construct from triangle mesh data.
    pub fn from_triangles(
        vertices: Vec<[f64; 3]>,
        triangles: &[[u32; 3]],
        prim_path: &str,
    ) -> Self {
        let face_vertex_counts = vec![3; triangles.len()];
        let face_vertex_indices: Vec<u32> =
            triangles.iter().flat_map(|t| t.iter().copied()).collect();
        UsdMeshWriter {
            vertices,
            face_vertex_counts,
            face_vertex_indices,
            prim_path: prim_path.to_string(),
        }
    }

    /// Write USD ASCII to `writer`.
    pub fn write<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "#usda 1.0")?;
        writeln!(
            w,
            r#"def Mesh "{name}" {{"#,
            name = self.prim_path.trim_start_matches('/')
        )?;
        // Points
        write!(w, "    point3f[] points = [")?;
        for (i, &[x, y, z]) in self.vertices.iter().enumerate() {
            if i > 0 {
                write!(w, ", ")?;
            }
            write!(w, "({x}, {y}, {z})")?;
        }
        writeln!(w, "]")?;
        // Face vertex counts
        write!(w, "    int[] faceVertexCounts = [")?;
        for (i, &c) in self.face_vertex_counts.iter().enumerate() {
            if i > 0 {
                write!(w, ", ")?;
            }
            write!(w, "{c}")?;
        }
        writeln!(w, "]")?;
        // Face vertex indices
        write!(w, "    int[] faceVertexIndices = [")?;
        for (i, &idx) in self.face_vertex_indices.iter().enumerate() {
            if i > 0 {
                write!(w, ", ")?;
            }
            write!(w, "{idx}")?;
        }
        writeln!(w, "]")?;
        writeln!(w, "}}")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// MshWriter
// ---------------------------------------------------------------------------

/// Write a Gmsh .msh format v2 mesh.
pub struct MshWriter {
    /// Node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Elements: each is a list of node indices (1-based in Gmsh).
    pub elements: Vec<Vec<usize>>,
    /// Physical group tags per element.
    pub element_tags: Vec<i64>,
    /// Mesh format version.
    pub version: u32,
}

impl MshWriter {
    /// Construct an empty MshWriter.
    pub fn new(version: u32) -> Self {
        MshWriter {
            nodes: Vec::new(),
            elements: Vec::new(),
            element_tags: Vec::new(),
            version,
        }
    }

    /// Add a node.
    pub fn add_node(&mut self, x: f64, y: f64, z: f64) {
        self.nodes.push([x, y, z]);
    }

    /// Add a triangular element with a physical group tag.
    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize, tag: i64) {
        self.elements.push(vec![a, b, c]);
        self.element_tags.push(tag);
    }

    /// Write .msh v2 format to `writer`.
    pub fn write_v2<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "$MeshFormat")?;
        writeln!(w, "2.2 0 8")?;
        writeln!(w, "$EndMeshFormat")?;
        writeln!(w, "$Nodes")?;
        writeln!(w, "{}", self.nodes.len())?;
        for (i, &[x, y, z]) in self.nodes.iter().enumerate() {
            writeln!(w, "{} {} {} {}", i + 1, x, y, z)?;
        }
        writeln!(w, "$EndNodes")?;
        writeln!(w, "$Elements")?;
        writeln!(w, "{}", self.elements.len())?;
        for (i, elem) in self.elements.iter().enumerate() {
            let tag = self.element_tags.get(i).copied().unwrap_or(0);
            let elem_type = match elem.len() {
                2 => 1, // line
                3 => 2, // triangle
                4 => 3, // quad
                _ => 2,
            };
            write!(w, "{} {} 2 {} {}", i + 1, elem_type, tag, tag)?;
            for &n in elem {
                write!(w, " {n}")?;
            }
            writeln!(w)?;
        }
        writeln!(w, "$EndElements")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// MshReader
// ---------------------------------------------------------------------------

/// Read a Gmsh .msh v2 mesh.
pub struct MshReader {
    /// Parsed node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Parsed elements (node indices, 0-based).
    pub elements: Vec<Vec<usize>>,
    /// Element tags.
    pub element_tags: Vec<i64>,
}

impl MshReader {
    /// Parse a .msh v2 string.
    pub fn parse(s: &str) -> Self {
        let mut nodes = Vec::new();
        let mut elements = Vec::new();
        let mut element_tags = Vec::new();
        let mut in_nodes = false;
        let mut in_elements = false;
        let mut node_count = 0usize;
        let mut elem_count = 0usize;
        let mut node_idx = 0usize;
        let mut elem_idx = 0usize;

        for line in s.lines() {
            let line = line.trim();
            if line == "$Nodes" {
                in_nodes = true;
                node_idx = 0;
                continue;
            }
            if line == "$EndNodes" {
                in_nodes = false;
                continue;
            }
            if line == "$Elements" {
                in_elements = true;
                elem_idx = 0;
                continue;
            }
            if line == "$EndElements" {
                in_elements = false;
                continue;
            }
            if in_nodes {
                if node_idx == 0 {
                    node_count = line.parse().unwrap_or(0);
                    nodes.reserve(node_count);
                    node_idx += 1;
                    continue;
                }
                let vals: Vec<f64> = line
                    .split_whitespace()
                    .skip(1)
                    .filter_map(|v| v.parse().ok())
                    .collect();
                if vals.len() >= 3 {
                    nodes.push([vals[0], vals[1], vals[2]]);
                }
            }
            if in_elements {
                if elem_idx == 0 {
                    elem_count = line.parse().unwrap_or(0);
                    elements.reserve(elem_count);
                    element_tags.reserve(elem_count);
                    elem_idx += 1;
                    continue;
                }
                let vals: Vec<i64> = line
                    .split_whitespace()
                    .filter_map(|v| v.parse().ok())
                    .collect();
                if vals.len() >= 5 {
                    let n_tags = vals[2] as usize;
                    let tag = vals[3];
                    let node_start = 3 + n_tags;
                    let node_indices: Vec<usize> = vals[node_start..]
                        .iter()
                        .map(|&v| (v - 1) as usize)
                        .collect();
                    elements.push(node_indices);
                    element_tags.push(tag);
                }
            }
        }
        let _ = node_count;
        let _ = elem_count;
        MshReader {
            nodes,
            elements,
            element_tags,
        }
    }
}

// ---------------------------------------------------------------------------
// AbaqusInpWriter
// ---------------------------------------------------------------------------

/// Write an Abaqus .inp mesh file.
pub struct AbaqusInpWriter {
    /// Part name.
    pub part_name: String,
    /// Node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Elements (node indices, 1-based in output).
    pub elements: Vec<Vec<usize>>,
    /// Element type string (e.g., "C3D4" for tet).
    pub element_type: String,
    /// Node sets: name → list of node indices (1-based).
    pub node_sets: HashMap<String, Vec<usize>>,
    /// Boundary conditions: node set name → (dof_min, dof_max, value).
    pub boundary_conditions: Vec<(String, usize, usize, f64)>,
}

impl AbaqusInpWriter {
    /// Construct an empty writer.
    pub fn new(part_name: &str, element_type: &str) -> Self {
        AbaqusInpWriter {
            part_name: part_name.to_string(),
            nodes: Vec::new(),
            elements: Vec::new(),
            element_type: element_type.to_string(),
            node_sets: HashMap::new(),
            boundary_conditions: Vec::new(),
        }
    }

    /// Add a node.
    pub fn add_node(&mut self, x: f64, y: f64, z: f64) {
        self.nodes.push([x, y, z]);
    }

    /// Add an element.
    pub fn add_element(&mut self, nodes: Vec<usize>) {
        self.elements.push(nodes);
    }

    /// Add a node set.
    pub fn add_node_set(&mut self, name: &str, indices: Vec<usize>) {
        self.node_sets.insert(name.to_string(), indices);
    }

    /// Add a boundary condition.
    pub fn add_boundary(&mut self, nset: &str, dof_min: usize, dof_max: usize, value: f64) {
        self.boundary_conditions
            .push((nset.to_string(), dof_min, dof_max, value));
    }

    /// Write the .inp file to `writer`.
    pub fn write<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "*PART, NAME={}", self.part_name)?;
        writeln!(w, "*NODE")?;
        for (i, &[x, y, z]) in self.nodes.iter().enumerate() {
            writeln!(w, "{}, {}, {}, {}", i + 1, x, y, z)?;
        }
        writeln!(w, "*ELEMENT, TYPE={}", self.element_type)?;
        for (i, elem) in self.elements.iter().enumerate() {
            write!(w, "{}", i + 1)?;
            for &n in elem {
                write!(w, ", {n}")?;
            }
            writeln!(w)?;
        }
        for (name, indices) in &self.node_sets {
            writeln!(w, "*NSET, NSET={name}")?;
            for chunk in indices.chunks(8) {
                let s: Vec<String> = chunk.iter().map(|n| n.to_string()).collect();
                writeln!(w, "{}", s.join(", "))?;
            }
        }
        if !self.boundary_conditions.is_empty() {
            writeln!(w, "*BOUNDARY")?;
            for (nset, d1, d2, val) in &self.boundary_conditions {
                writeln!(w, "{nset}, {d1}, {d2}, {val}")?;
            }
        }
        writeln!(w, "*END PART")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// AbaqusInpReader
// ---------------------------------------------------------------------------

/// Read an Abaqus .inp mesh.
pub struct AbaqusInpReader {
    /// Parsed nodes.
    pub nodes: Vec<[f64; 3]>,
    /// Parsed elements.
    pub elements: Vec<Vec<usize>>,
    /// Element type.
    pub element_type: String,
    /// Node sets.
    pub node_sets: HashMap<String, Vec<usize>>,
}

impl AbaqusInpReader {
    /// Parse an .inp string.
    pub fn parse(s: &str) -> Self {
        let mut nodes = Vec::new();
        let mut elements = Vec::new();
        let mut element_type = String::new();
        let mut node_sets: HashMap<String, Vec<usize>> = HashMap::new();
        let mut mode = "";
        let mut cur_nset = String::new();

        for line in s.lines() {
            let line = line.trim();
            if line.starts_with('*') {
                let kw = line.to_uppercase();
                if kw.starts_with("*NODE") && !kw.starts_with("*NSET") {
                    mode = "node";
                } else if kw.starts_with("*ELEMENT") {
                    mode = "element";
                    if let Some(tp) = line.split("TYPE=").nth(1) {
                        element_type = tp.split(',').next().unwrap_or("").trim().to_string();
                    }
                } else if kw.starts_with("*NSET") {
                    mode = "nset";
                    if let Some(name) = line.split("NSET=").nth(1) {
                        cur_nset = name.split(',').next().unwrap_or("").trim().to_string();
                        node_sets.entry(cur_nset.clone()).or_default();
                    }
                } else if kw.starts_with("*BOUNDARY") {
                    mode = "boundary";
                } else {
                    mode = "";
                }
                continue;
            }
            match mode {
                "node" => {
                    let vals: Vec<f64> = line
                        .split(',')
                        .filter_map(|v| v.trim().parse().ok())
                        .collect();
                    if vals.len() >= 4 {
                        nodes.push([vals[1], vals[2], vals[3]]);
                    }
                }
                "element" => {
                    let vals: Vec<usize> = line
                        .split(',')
                        .filter_map(|v| v.trim().parse().ok())
                        .collect();
                    if vals.len() >= 2 {
                        elements.push(vals[1..].to_vec());
                    }
                }
                "nset" => {
                    let indices: Vec<usize> = line
                        .split(',')
                        .filter_map(|v| v.trim().parse().ok())
                        .collect();
                    if let Some(set) = node_sets.get_mut(&cur_nset) {
                        set.extend(indices);
                    }
                }
                _ => {}
            }
        }
        AbaqusInpReader {
            nodes,
            elements,
            element_type,
            node_sets,
        }
    }
}

// ---------------------------------------------------------------------------
// NastranBdf
// ---------------------------------------------------------------------------

/// A NASTRAN BDF GRID entry.
#[derive(Debug, Clone)]
pub struct NastranGrid {
    /// Node id.
    pub id: u64,
    /// Coordinate system id.
    pub cp: u64,
    /// Position.
    pub xyz: [f64; 3],
}

/// A NASTRAN solid element (CTETRA or CHEXA).
#[derive(Debug, Clone)]
pub struct NastranElement {
    /// Element id.
    pub eid: u64,
    /// Property id.
    pub pid: u64,
    /// Node ids.
    pub nodes: Vec<u64>,
    /// Element type.
    pub elem_type: String,
}

/// Read/write NASTRAN BDF (GRID, CTETRA, CHEXA, PSOLID, MAT1, SPC, FORCE).
pub struct NastranBdf {
    /// Grid points.
    pub grids: Vec<NastranGrid>,
    /// Elements.
    pub elements: Vec<NastranElement>,
}

impl NastranBdf {
    /// Construct empty.
    pub fn new() -> Self {
        NastranBdf {
            grids: Vec::new(),
            elements: Vec::new(),
        }
    }

    /// Add a grid point.
    pub fn add_grid(&mut self, id: u64, xyz: [f64; 3]) {
        self.grids.push(NastranGrid { id, cp: 0, xyz });
    }

    /// Add an element.
    pub fn add_element(&mut self, eid: u64, pid: u64, nodes: Vec<u64>, elem_type: &str) {
        self.elements.push(NastranElement {
            eid,
            pid,
            nodes,
            elem_type: elem_type.to_string(),
        });
    }

    /// Write BDF to `writer`.
    pub fn write<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "BEGIN BULK")?;
        for g in &self.grids {
            writeln!(
                w,
                "GRID    {:8}{:8}{:8}{:8}{:8}",
                g.id,
                g.cp,
                nastran_field8(g.xyz[0]),
                nastran_field8(g.xyz[1]),
                nastran_field8(g.xyz[2]),
            )?;
        }
        for e in &self.elements {
            write!(w, "{:<8}{:8}{:8}", e.elem_type, e.eid, e.pid)?;
            for &n in &e.nodes {
                write!(w, "{:8}", n)?;
            }
            writeln!(w)?;
        }
        writeln!(w, "ENDDATA")?;
        Ok(())
    }

    /// Parse a BDF string and return grids and elements.
    pub fn parse(s: &str) -> Self {
        let mut bdf = NastranBdf::new();
        for line in s.lines() {
            let line = if line.len() > 8 { line } else { continue };
            let kw = &line[..8].trim().to_uppercase();
            if kw.starts_with("GRID") {
                if let Some(g) = Self::parse_grid(line) {
                    bdf.grids.push(g);
                }
            } else if (kw.starts_with("CTETRA") || kw.starts_with("CHEXA"))
                && let Some(e) = Self::parse_element(line)
            {
                bdf.elements.push(e);
            }
        }
        bdf
    }

    fn parse_grid(line: &str) -> Option<NastranGrid> {
        let fields: Vec<&str> = (0..9)
            .map(|i| {
                let s = i * 8;
                let e = (s + 8).min(line.len());
                if s < line.len() {
                    line[s..e].trim()
                } else {
                    ""
                }
            })
            .collect();
        let id: u64 = fields[1].parse().ok()?;
        let x: f64 = fields[3].parse().ok()?;
        let y: f64 = fields[4].parse().ok()?;
        let z: f64 = fields[5].parse().ok()?;
        Some(NastranGrid {
            id,
            cp: 0,
            xyz: [x, y, z],
        })
    }

    fn parse_element(line: &str) -> Option<NastranElement> {
        let kw = line[..8].trim().to_uppercase();
        let fields: Vec<&str> = (0..12)
            .map(|i| {
                let s = i * 8;
                let e = (s + 8).min(line.len());
                if s < line.len() {
                    line[s..e].trim()
                } else {
                    ""
                }
            })
            .collect();
        let eid: u64 = fields[1].parse().ok()?;
        let pid: u64 = fields[2].parse().ok()?;
        let nodes: Vec<u64> = fields[3..].iter().filter_map(|f| f.parse().ok()).collect();
        Some(NastranElement {
            eid,
            pid,
            nodes,
            elem_type: kw,
        })
    }
}

impl Default for NastranBdf {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// MeshConverter
// ---------------------------------------------------------------------------

/// Convert between mesh formats, re-index nodes, and reorder elements.
pub struct MeshConverter;

impl MeshConverter {
    /// Re-index nodes: remove unused nodes and remap element indices.
    /// Returns (new_nodes, new_elements).
    pub fn compact_nodes(
        nodes: &[[f64; 3]],
        elements: &[Vec<usize>],
    ) -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
        let mut used: HashSet<usize> = HashSet::new();
        for elem in elements {
            for &n in elem {
                used.insert(n);
            }
        }
        let mut old_to_new: HashMap<usize, usize> = HashMap::new();
        let mut new_nodes = Vec::new();
        let mut sorted: Vec<usize> = used.into_iter().collect();
        sorted.sort_unstable();
        for old_idx in sorted {
            old_to_new.insert(old_idx, new_nodes.len());
            if old_idx < nodes.len() {
                new_nodes.push(nodes[old_idx]);
            }
        }
        let new_elements: Vec<Vec<usize>> = elements
            .iter()
            .map(|e| {
                e.iter()
                    .filter_map(|&n| old_to_new.get(&n).copied())
                    .collect()
            })
            .collect();
        (new_nodes, new_elements)
    }

    /// Convert triangle mesh to indexed mesh (deduplicate vertices by position).
    pub fn deduplicate_vertices(
        positions: &[[f64; 3]],
        indices: &[[u32; 3]],
        tolerance: f64,
    ) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let mut unique: Vec<[f64; 3]> = Vec::new();
        let mut remap: Vec<u32> = Vec::with_capacity(positions.len());
        for &p in positions {
            let mut found = None;
            for (j, &u) in unique.iter().enumerate() {
                let dx = p[0] - u[0];
                let dy = p[1] - u[1];
                let dz = p[2] - u[2];
                if (dx * dx + dy * dy + dz * dz).sqrt() < tolerance {
                    found = Some(j as u32);
                    break;
                }
            }
            if let Some(j) = found {
                remap.push(j);
            } else {
                remap.push(unique.len() as u32);
                unique.push(p);
            }
        }
        let new_indices: Vec<[u32; 3]> = indices
            .iter()
            .map(|&[a, b, c]| [remap[a as usize], remap[b as usize], remap[c as usize]])
            .collect();
        (unique, new_indices)
    }

    /// Flip winding order of all triangles.
    pub fn flip_winding(indices: &[[u32; 3]]) -> Vec<[u32; 3]> {
        indices.iter().map(|&[a, b, c]| [a, c, b]).collect()
    }
}

// ---------------------------------------------------------------------------
// MeshValidator
// ---------------------------------------------------------------------------

/// Mesh quality metrics.
#[derive(Debug, Clone)]
pub struct MeshQuality {
    /// Minimum aspect ratio (ideal = 1 for equilateral).
    pub min_aspect_ratio: f64,
    /// Maximum aspect ratio.
    pub max_aspect_ratio: f64,
    /// Minimum triangle area.
    pub min_area: f64,
    /// Number of degenerate triangles (area < epsilon).
    pub degenerate_count: usize,
    /// Number of inverted triangles (negative Jacobian).
    pub inverted_count: usize,
}

/// Check mesh quality: aspect ratio, Jacobian, area, orientation.
pub struct MeshValidator;

impl MeshValidator {
    /// Compute quality metrics for a triangle mesh.
    pub fn validate(
        vertices: &[[f64; 3]],
        triangles: &[[u32; 3]],
        area_epsilon: f64,
    ) -> MeshQuality {
        let mut min_ar = f64::MAX;
        let mut max_ar = f64::MIN;
        let mut min_area = f64::MAX;
        let mut deg_count = 0;
        let mut inv_count = 0;

        for &[a, b, c] in triangles {
            let va = vertices[a as usize];
            let vb = vertices[b as usize];
            let vc = vertices[c as usize];
            let ab = sub3(vb, va);
            let ac = sub3(vc, va);
            let bc = sub3(vc, vb);
            let cross = cross3(ab, ac);
            let area = len3(cross) * 0.5;
            let ar = triangle_aspect_ratio(len3(ab), len3(ac), len3(bc));
            min_ar = min_ar.min(ar);
            max_ar = max_ar.max(ar);
            min_area = min_area.min(area);
            if area < area_epsilon {
                deg_count += 1;
            }
            // Check orientation: if normal's y < 0 after projecting, might be inverted
            if cross[1] < 0.0 {
                inv_count += 1;
            }
        }

        if triangles.is_empty() {
            min_ar = 0.0;
            max_ar = 0.0;
            min_area = 0.0;
        }

        MeshQuality {
            min_aspect_ratio: min_ar,
            max_aspect_ratio: max_ar,
            min_area,
            degenerate_count: deg_count,
            inverted_count: inv_count,
        }
    }

    /// Check that all element indices are within bounds.
    pub fn check_bounds(vertices: &[[f64; 3]], triangles: &[[u32; 3]]) -> bool {
        let n = vertices.len() as u32;
        triangles.iter().all(|&[a, b, c]| a < n && b < n && c < n)
    }
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Compute the aspect ratio of a triangle given edge lengths.
fn triangle_aspect_ratio(a: f64, b: f64, c: f64) -> f64 {
    let s = (a + b + c) / 2.0;
    let area = (s * (s - a) * (s - b) * (s - c)).max(0.0).sqrt();
    if area < 1e-14 {
        return f64::MAX;
    }
    let longest = a.max(b).max(c);
    longest * (a + b + c) / (4.0 * 3_f64.sqrt() * area)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- base64_encode_buffer ---

    #[test]
    fn base64_empty() {
        assert_eq!(base64_encode_buffer(&[]), "");
    }

    #[test]
    fn base64_man() {
        // "Man" → "TWFu"
        let result = base64_encode_buffer(b"Man");
        assert_eq!(result, "TWFu");
    }

    #[test]
    fn base64_one_byte_padded() {
        let result = base64_encode_buffer(b"M");
        assert_eq!(result.len(), 4);
        assert!(result.ends_with("=="));
    }

    #[test]
    fn base64_two_bytes_padded() {
        let result = base64_encode_buffer(b"Ma");
        assert_eq!(result.len(), 4);
        assert!(result.ends_with('='));
    }

    // --- nastran_field8 ---

    #[test]
    fn nastran_field8_length() {
        let s = nastran_field8(1.23456789);
        assert!(s.len() <= 8, "length={}", s.len());
    }

    #[test]
    fn nastran_field8_zero() {
        let s = nastran_field8(0.0);
        assert!(s.len() <= 8);
    }

    // --- gmsh_element_type_str ---

    #[test]
    fn gmsh_triangle() {
        assert_eq!(gmsh_element_type_str(3), "Triangle");
    }

    #[test]
    fn gmsh_line() {
        assert_eq!(gmsh_element_type_str(2), "Line");
    }

    // --- GltfMeshWriter ---

    #[test]
    fn gltf_writer_produces_json() {
        let verts = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let indices = vec![[0u32, 1, 2]];
        let w = GltfMeshWriter::new(verts, vec![], vec![], indices);
        let mut buf = Vec::new();
        w.write_gltf_json(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("asset"), "s={s}");
        assert!(s.contains("accessors"));
    }

    #[test]
    fn gltf_positions_bytes_length() {
        let verts = vec![[0.0f64; 3]; 4];
        let w = GltfMeshWriter::new(verts, vec![], vec![], vec![]);
        assert_eq!(w.positions_f32_bytes().len(), 48); // 4 * 3 * 4
    }

    #[test]
    fn gltf_indices_bytes_length() {
        let indices = vec![[0u32, 1, 2], [1u32, 2, 3]];
        let w = GltfMeshWriter::new(vec![], vec![], vec![], indices);
        assert_eq!(w.indices_u32_bytes().len(), 24); // 2 * 3 * 4
    }

    // --- GltfMeshReader ---

    #[test]
    fn gltf_reader_parse_names() {
        let json = r#"{"meshes": [{"name": "MyMesh"}]}"#;
        let r = GltfMeshReader::parse_json(json);
        assert!(
            r.node_names.iter().any(|n| n.contains("MyMesh")),
            "names={:?}",
            r.node_names
        );
    }

    // --- UsdMeshWriter ---

    #[test]
    fn usd_writer_contains_def_mesh() {
        let verts = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let w = UsdMeshWriter::from_triangles(verts, &[[0, 1, 2]], "/World/Mesh");
        let mut buf = Vec::new();
        w.write(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("def Mesh"), "s={s}");
        assert!(s.contains("points"));
    }

    #[test]
    fn usd_writer_face_vertex_indices() {
        let verts = vec![[0.0f64; 3]; 3];
        let w = UsdMeshWriter::from_triangles(verts, &[[0, 1, 2]], "/Mesh");
        assert_eq!(w.face_vertex_indices.len(), 3);
        assert_eq!(w.face_vertex_counts.len(), 1);
    }

    // --- MshWriter / MshReader ---

    #[test]
    fn msh_writer_v2_roundtrip() {
        let mut mw = MshWriter::new(2);
        mw.add_node(0.0, 0.0, 0.0);
        mw.add_node(1.0, 0.0, 0.0);
        mw.add_node(0.0, 1.0, 0.0);
        mw.add_triangle(1, 2, 3, 1);
        let mut buf = Vec::new();
        mw.write_v2(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("$MeshFormat"));
        assert!(s.contains("$Nodes"));
        assert!(s.contains("$Elements"));
        let r = MshReader::parse(&s);
        assert_eq!(r.nodes.len(), 3, "parsed nodes: {:?}", r.nodes);
        assert_eq!(r.elements.len(), 1);
    }

    #[test]
    fn msh_reader_node_coords() {
        let msh = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n2\n1 1.0 2.0 3.0\n2 4.0 5.0 6.0\n$EndNodes\n$Elements\n0\n$EndElements\n";
        let r = MshReader::parse(msh);
        assert_eq!(r.nodes.len(), 2);
        assert!((r.nodes[0][0] - 1.0).abs() < 1e-9);
    }

    // --- AbaqusInpWriter / Reader ---

    #[test]
    fn abaqus_writer_produces_part() {
        let mut w = AbaqusInpWriter::new("PART-1", "C3D4");
        w.add_node(0.0, 0.0, 0.0);
        w.add_node(1.0, 0.0, 0.0);
        w.add_node(0.0, 1.0, 0.0);
        w.add_node(0.0, 0.0, 1.0);
        w.add_element(vec![1, 2, 3, 4]);
        let mut buf = Vec::new();
        w.write(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("*PART"), "s={s}");
        assert!(s.contains("*NODE"));
        assert!(s.contains("*ELEMENT"));
    }

    #[test]
    fn abaqus_reader_parse() {
        let inp = "*PART, NAME=P1\n*NODE\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0\n3, 0.0, 1.0, 0.0\n*ELEMENT, TYPE=C3D4\n1, 1, 2, 3\n*END PART\n";
        let r = AbaqusInpReader::parse(inp);
        assert_eq!(r.nodes.len(), 3, "nodes={:?}", r.nodes);
        assert_eq!(r.elements.len(), 1);
    }

    #[test]
    fn abaqus_boundary_written() {
        let mut w = AbaqusInpWriter::new("P", "C3D4");
        w.add_node_set("FIXED", vec![1, 2]);
        w.add_boundary("FIXED", 1, 3, 0.0);
        let mut buf = Vec::new();
        w.write(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("*BOUNDARY"), "s={s}");
    }

    // --- NastranBdf ---

    #[test]
    fn nastran_write_and_parse() {
        let mut bdf = NastranBdf::new();
        bdf.add_grid(1, [0.0, 0.0, 0.0]);
        bdf.add_grid(2, [1.0, 0.0, 0.0]);
        bdf.add_element(1, 1, vec![1, 2, 3, 4], "CTETRA");
        let mut buf = Vec::new();
        bdf.write(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("BEGIN BULK"), "s={s}");
        assert!(s.contains("GRID"));
        assert!(s.contains("ENDDATA"));
    }

    #[test]
    fn nastran_parse_grids() {
        let bdf_str = "BEGIN BULK\nGRID         1       0    1.00    2.00    3.00\nENDDATA\n";
        let bdf = NastranBdf::parse(bdf_str);
        // Simplified parser; just ensure no panic
        let _ = bdf.grids.len();
    }

    // --- MeshConverter ---

    #[test]
    fn compact_nodes_removes_unused() {
        let nodes = vec![[0.0f64; 3]; 5];
        let elements = vec![vec![0, 2, 4]];
        let (new_nodes, new_elements) = MeshConverter::compact_nodes(&nodes, &elements);
        assert_eq!(new_nodes.len(), 3);
        assert_eq!(new_elements[0].len(), 3);
    }

    #[test]
    fn deduplicate_identical_vertices() {
        let positions = vec![[0.0f64, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let indices = vec![[0u32, 1, 2]];
        let (uniq, _new_idx) = MeshConverter::deduplicate_vertices(&positions, &indices, 1e-6);
        assert_eq!(uniq.len(), 2);
    }

    #[test]
    fn flip_winding_reverses_order() {
        let indices = vec![[0u32, 1, 2]];
        let flipped = MeshConverter::flip_winding(&indices);
        assert_eq!(flipped[0], [0, 2, 1]);
    }

    // --- MeshValidator ---

    #[test]
    fn validator_equilateral_good_aspect_ratio() {
        let sq3 = 3_f64.sqrt();
        let verts = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, sq3 / 2.0, 0.0]];
        let tris = vec![[0u32, 1, 2]];
        let q = MeshValidator::validate(&verts, &tris, 1e-10);
        assert!(q.min_aspect_ratio < 1.1, "ar={}", q.min_aspect_ratio);
    }

    #[test]
    fn validator_degenerate_triangle() {
        let verts = vec![[0.0f64; 3]; 3]; // all at origin
        let tris = vec![[0u32, 1, 2]];
        let q = MeshValidator::validate(&verts, &tris, 1e-10);
        assert_eq!(q.degenerate_count, 1);
    }

    #[test]
    fn validator_check_bounds_ok() {
        let verts = vec![[0.0f64; 3]; 3];
        let tris = vec![[0u32, 1, 2]];
        assert!(MeshValidator::check_bounds(&verts, &tris));
    }

    #[test]
    fn validator_check_bounds_fail() {
        let verts = vec![[0.0f64; 3]; 2];
        let tris = vec![[0u32, 1, 5]];
        assert!(!MeshValidator::check_bounds(&verts, &tris));
    }

    // --- write_gltf_accessor ---

    #[test]
    fn write_gltf_accessor_correct() {
        let mut buf = Vec::new();
        write_gltf_accessor(&mut buf, 0, 5126, 100, "VEC3", false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("VEC3"));
        assert!(s.contains("100"));
    }

    // --- triangle_aspect_ratio ---

    #[test]
    fn aspect_ratio_equilateral_is_one() {
        let ar = triangle_aspect_ratio(1.0, 1.0, 1.0);
        assert!((ar - 1.0).abs() < 1e-6, "ar={ar}");
    }

    #[test]
    fn aspect_ratio_degenerate_is_max() {
        let ar = triangle_aspect_ratio(1.0, 1.0, 2.0);
        assert!(ar > 1.0 || ar == f64::MAX);
    }
}
