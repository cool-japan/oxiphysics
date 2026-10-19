// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh export and import utilities.
//!
//! Provides:
//! - `MeshData`: vertex positions, face indices, normals, UVs, material IDs.
//! - `ObjExporter`: write Wavefront OBJ + MTL files.
//! - `StlExporter`: binary and ASCII STL export with normals.
//! - `PlyExporter`: PLY header + binary/ASCII body.
//! - `GltfExporter`: glTF 2.0 JSON + binary buffer.
//! - `MeshImporter`: parse OBJ, binary/ASCII STL, and PLY.
//! - `MeshSimplification`: quadric error metric edge-collapse simplification.

use std::fmt::Write as FmtWrite;

mod qem;

// ============================================================================
// MeshData
// ============================================================================

/// Container for mesh geometry and per-vertex attributes.
#[derive(Debug, Clone, Default)]
pub struct MeshData {
    /// Vertex positions as `[f64; 3]` triples.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle face indices: each entry is `[v0, v1, v2]`.
    pub faces: Vec<[usize; 3]>,
    /// Per-vertex normals (may be empty).
    pub normals: Vec<[f64; 3]>,
    /// Per-vertex UV coordinates (may be empty).
    pub uvs: Vec<[f64; 2]>,
    /// Per-face material ID (may be empty).
    pub material_ids: Vec<usize>,
}

impl MeshData {
    /// Create an empty `MeshData`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        idx
    }

    /// Add a triangular face.
    pub fn add_face(&mut self, v0: usize, v1: usize, v2: usize) {
        self.faces.push([v0, v1, v2]);
    }

    /// Add a triangular face with material ID.
    pub fn add_face_with_material(&mut self, v0: usize, v1: usize, v2: usize, mat: usize) {
        self.faces.push([v0, v1, v2]);
        self.material_ids.push(mat);
    }

    /// Compute per-face normals and store as per-vertex normals (flat shading).
    pub fn compute_flat_normals(&mut self) {
        self.normals = vec![[0.0; 3]; self.vertices.len()];
        for face in &self.faces {
            let v0 = self.vertices[face[0]];
            let v1 = self.vertices[face[1]];
            let v2 = self.vertices[face[2]];
            let n = triangle_normal_f64(v0, v1, v2);
            for &vi in face {
                self.normals[vi] = n;
            }
        }
    }

    /// Compute smooth per-vertex normals by averaging face normals.
    pub fn compute_smooth_normals(&mut self) {
        let mut acc = vec![[0.0f64; 3]; self.vertices.len()];
        let mut count = vec![0u32; self.vertices.len()];
        for face in &self.faces {
            let v0 = self.vertices[face[0]];
            let v1 = self.vertices[face[1]];
            let v2 = self.vertices[face[2]];
            let n = triangle_normal_f64(v0, v1, v2);
            for &vi in face {
                acc[vi][0] += n[0];
                acc[vi][1] += n[1];
                acc[vi][2] += n[2];
                count[vi] += 1;
            }
        }
        self.normals = acc
            .iter()
            .zip(count.iter())
            .map(|(a, &c)| {
                if c == 0 {
                    [0.0; 3]
                } else {
                    let inv = 1.0 / c as f64;
                    normalize3_f64([a[0] * inv, a[1] * inv, a[2] * inv])
                }
            })
            .collect();
    }

    /// Return the axis-aligned bounding box of all vertices as `(min, max)`.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for &v in &self.vertices {
            for i in 0..3 {
                mn[i] = mn[i].min(v[i]);
                mx[i] = mx[i].max(v[i]);
            }
        }
        (mn, mx)
    }

    /// Number of triangles in the mesh.
    pub fn num_triangles(&self) -> usize {
        self.faces.len()
    }
}

// ============================================================================
// Geometry helpers
// ============================================================================

/// Compute the unit normal of a triangle given three vertices (f64).
fn triangle_normal_f64(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    let a = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let b = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    normalize3_f64(n)
}

/// Normalise a 3-D f64 vector; returns zero vector if length < eps.
fn normalize3_f64(v: [f64; 3]) -> [f64; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if l < 1e-12 {
        [0.0; 3]
    } else {
        [v[0] / l, v[1] / l, v[2] / l]
    }
}

// ============================================================================
// Material
// ============================================================================

/// Simple Phong/PBR material description used by exporters.
#[derive(Debug, Clone)]
pub struct Material {
    /// Material name.
    pub name: String,
    /// Diffuse / base colour (RGB in \[0,1\]).
    pub diffuse: [f64; 3],
    /// Specular colour (RGB in \[0,1\]).
    pub specular: [f64; 3],
    /// Shininess / roughness exponent.
    pub shininess: f64,
    /// Opacity (1.0 = fully opaque).
    pub opacity: f64,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            name: "default".to_string(),
            diffuse: [0.8, 0.8, 0.8],
            specular: [0.2, 0.2, 0.2],
            shininess: 32.0,
            opacity: 1.0,
        }
    }
}

impl Material {
    /// Create a new material with the given name and default parameters.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }
}

// ============================================================================
// ObjExporter
// ============================================================================

/// Wavefront OBJ + MTL exporter.
#[derive(Debug, Default)]
pub struct ObjExporter {
    /// Materials to embed in the `.mtl` file.
    pub materials: Vec<Material>,
}

impl ObjExporter {
    /// Create a new `ObjExporter`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a material.
    pub fn add_material(&mut self, mat: Material) {
        self.materials.push(mat);
    }

    /// Serialise `mesh` to OBJ format and return the string.
    ///
    /// If `mtl_name` is provided a `mtllib` directive is emitted.
    pub fn export_obj(&self, mesh: &MeshData, mtl_name: Option<&str>) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "# OxiPhysics OBJ export");
        if let Some(name) = mtl_name {
            let _ = writeln!(out, "mtllib {name}");
        }

        // Vertices
        for &v in &mesh.vertices {
            let _ = writeln!(out, "v {} {} {}", v[0], v[1], v[2]);
        }

        // Normals
        for &n in &mesh.normals {
            let _ = writeln!(out, "vn {} {} {}", n[0], n[1], n[2]);
        }

        // UVs
        for &uv in &mesh.uvs {
            let _ = writeln!(out, "vt {} {}", uv[0], uv[1]);
        }

        // Group faces by material
        let has_mats = !mesh.material_ids.is_empty() && !self.materials.is_empty();
        let mut current_mat: Option<usize> = None;

        for (fi, &face) in mesh.faces.iter().enumerate() {
            if has_mats {
                let mat_id = *mesh.material_ids.get(fi).unwrap_or(&0);
                if current_mat != Some(mat_id) {
                    if let Some(mat) = self.materials.get(mat_id) {
                        let _ = writeln!(out, "usemtl {}", mat.name);
                    }
                    current_mat = Some(mat_id);
                }
            }
            // OBJ uses 1-based indices
            let v0 = face[0] + 1;
            let v1 = face[1] + 1;
            let v2 = face[2] + 1;
            if !mesh.normals.is_empty() && !mesh.uvs.is_empty() {
                let _ = writeln!(out, "f {v0}/{v0}/{v0} {v1}/{v1}/{v1} {v2}/{v2}/{v2}");
            } else if !mesh.normals.is_empty() {
                let _ = writeln!(out, "f {v0}//{v0} {v1}//{v1} {v2}//{v2}");
            } else if !mesh.uvs.is_empty() {
                let _ = writeln!(out, "f {v0}/{v0} {v1}/{v1} {v2}/{v2}");
            } else {
                let _ = writeln!(out, "f {v0} {v1} {v2}");
            }
        }
        out
    }

    /// Serialise the material list to MTL format.
    pub fn export_mtl(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "# OxiPhysics MTL export");
        for mat in &self.materials {
            let _ = writeln!(out, "newmtl {}", mat.name);
            writeln!(
                out,
                "Kd {} {} {}",
                mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]
            )
            .expect("operation should succeed");
            writeln!(
                out,
                "Ks {} {} {}",
                mat.specular[0], mat.specular[1], mat.specular[2]
            )
            .expect("operation should succeed");
            let _ = writeln!(out, "Ns {}", mat.shininess);
            let _ = writeln!(out, "d {}", mat.opacity);
        }
        out
    }
}

// ============================================================================
// StlExporter
// ============================================================================

/// Binary and ASCII STL exporter.
#[derive(Debug, Default)]
pub struct StlExporter;

impl StlExporter {
    /// Create a new `StlExporter`.
    pub fn new() -> Self {
        Self
    }

    /// Export `mesh` to binary STL bytes.
    ///
    /// Layout: 80-byte header, 4-byte count, 50 bytes per triangle.
    pub fn export_binary(&self, mesh: &MeshData, solid_name: &str) -> Vec<u8> {
        let n = mesh.faces.len();
        let mut buf = Vec::with_capacity(84 + n * 50);

        // 80-byte header
        let mut header = [0u8; 80];
        let name_bytes = solid_name.as_bytes();
        let copy_len = name_bytes.len().min(80);
        header[..copy_len].copy_from_slice(&name_bytes[..copy_len]);
        buf.extend_from_slice(&header);

        // Triangle count
        buf.extend_from_slice(&(n as u32).to_le_bytes());

        for &face in &mesh.faces {
            let v0 = mesh.vertices[face[0]];
            let v1 = mesh.vertices[face[1]];
            let v2 = mesh.vertices[face[2]];

            // Per-triangle normal from geometry (or from normals array)
            let n_vec = if !mesh.normals.is_empty() {
                mesh.normals[face[0]]
            } else {
                triangle_normal_f64(v0, v1, v2)
            };

            for &c in &n_vec {
                buf.extend_from_slice(&(c as f32).to_le_bytes());
            }
            for &v in &[v0, v1, v2] {
                for &c in &v {
                    buf.extend_from_slice(&(c as f32).to_le_bytes());
                }
            }
            buf.extend_from_slice(&0u16.to_le_bytes()); // attribute
        }
        buf
    }

    /// Export `mesh` to ASCII STL string.
    pub fn export_ascii(&self, mesh: &MeshData, solid_name: &str) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "solid {solid_name}");
        for &face in &mesh.faces {
            let v0 = mesh.vertices[face[0]];
            let v1 = mesh.vertices[face[1]];
            let v2 = mesh.vertices[face[2]];
            let n_vec = triangle_normal_f64(v0, v1, v2);
            let _ = writeln!(out, "  facet normal {} {} {}", n_vec[0], n_vec[1], n_vec[2]);
            let _ = writeln!(out, "    outer loop");
            let _ = writeln!(out, "      vertex {} {} {}", v0[0], v0[1], v0[2]);
            let _ = writeln!(out, "      vertex {} {} {}", v1[0], v1[1], v1[2]);
            let _ = writeln!(out, "      vertex {} {} {}", v2[0], v2[1], v2[2]);
            let _ = writeln!(out, "    endloop");
            let _ = writeln!(out, "  endfacet");
        }
        let _ = writeln!(out, "endsolid {solid_name}");
        out
    }

    /// Compute a simple CRC32 checksum of `data`.
    pub fn crc32(data: &[u8]) -> u32 {
        let mut crc: u32 = 0xFFFF_FFFF;
        for &byte in data {
            crc ^= byte as u32;
            for _ in 0..8 {
                if crc & 1 != 0 {
                    crc = (crc >> 1) ^ 0xEDB8_8320;
                } else {
                    crc >>= 1;
                }
            }
        }
        crc ^ 0xFFFF_FFFF
    }
}

// ============================================================================
// PlyExporter
// ============================================================================

/// PLY (Polygon File Format) exporter supporting binary and ASCII modes.
#[derive(Debug, Default)]
pub struct PlyExporter {
    /// Write in binary little-endian format when `true`; ASCII when `false`.
    pub binary: bool,
}

impl PlyExporter {
    /// Create a new `PlyExporter`.
    pub fn new(binary: bool) -> Self {
        Self { binary }
    }

    /// Build the PLY header string.
    fn build_header(&self, mesh: &MeshData) -> String {
        let mut h = String::new();
        let _ = writeln!(h, "ply");
        if self.binary {
            let _ = writeln!(h, "format binary_little_endian 1.0");
        } else {
            let _ = writeln!(h, "format ascii 1.0");
        }
        let _ = writeln!(h, "comment OxiPhysics PLY export");
        let _ = writeln!(h, "element vertex {}", mesh.vertices.len());
        let _ = writeln!(h, "property float x");
        let _ = writeln!(h, "property float y");
        let _ = writeln!(h, "property float z");
        if !mesh.normals.is_empty() {
            let _ = writeln!(h, "property float nx");
            let _ = writeln!(h, "property float ny");
            let _ = writeln!(h, "property float nz");
        }
        if !mesh.uvs.is_empty() {
            let _ = writeln!(h, "property float s");
            let _ = writeln!(h, "property float t");
        }
        let _ = writeln!(h, "element face {}", mesh.faces.len());
        let _ = writeln!(h, "property list uchar int vertex_indices");
        let _ = writeln!(h, "end_header");
        h
    }

    /// Export `mesh` to ASCII PLY bytes.
    pub fn export_ascii_bytes(&self, mesh: &MeshData) -> Vec<u8> {
        let header = self.build_header(mesh);
        let mut out = header;
        for (i, &v) in mesh.vertices.iter().enumerate() {
            let mut line = format!("{} {} {}", v[0] as f32, v[1] as f32, v[2] as f32);
            if !mesh.normals.is_empty() {
                let n = mesh.normals[i];
                line.push_str(&format!(" {} {} {}", n[0] as f32, n[1] as f32, n[2] as f32));
            }
            if !mesh.uvs.is_empty() {
                let uv = mesh.uvs[i];
                line.push_str(&format!(" {} {}", uv[0] as f32, uv[1] as f32));
            }
            let _ = writeln!(out, "{line}");
        }
        for &face in &mesh.faces {
            let _ = writeln!(out, "3 {} {} {}", face[0], face[1], face[2]);
        }
        out.into_bytes()
    }

    /// Export `mesh` to binary PLY bytes (little-endian).
    pub fn export_binary_bytes(&self, mesh: &MeshData) -> Vec<u8> {
        let header = self.build_header(mesh);
        let mut buf = header.into_bytes();

        for (i, &v) in mesh.vertices.iter().enumerate() {
            for &c in &v {
                buf.extend_from_slice(&(c as f32).to_le_bytes());
            }
            if !mesh.normals.is_empty() {
                let n = mesh.normals[i];
                for &c in &n {
                    buf.extend_from_slice(&(c as f32).to_le_bytes());
                }
            }
            if !mesh.uvs.is_empty() {
                let uv = mesh.uvs[i];
                for &c in &uv {
                    buf.extend_from_slice(&(c as f32).to_le_bytes());
                }
            }
        }
        for &face in &mesh.faces {
            buf.push(3u8); // list count
            for &vi in &face {
                buf.extend_from_slice(&(vi as i32).to_le_bytes());
            }
        }
        buf
    }

    /// Export using the mode selected at construction time.
    pub fn export(&self, mesh: &MeshData) -> Vec<u8> {
        if self.binary {
            self.export_binary_bytes(mesh)
        } else {
            self.export_ascii_bytes(mesh)
        }
    }
}

// ============================================================================
// GltfExporter
// ============================================================================

/// A minimal glTF 2.0 exporter producing JSON + binary buffer.
#[derive(Debug, Default)]
pub struct GltfExporter {
    /// Scene name embedded in the JSON.
    pub scene_name: String,
}

impl GltfExporter {
    /// Create a new `GltfExporter` with the given scene name.
    pub fn new(scene_name: impl Into<String>) -> Self {
        Self {
            scene_name: scene_name.into(),
        }
    }

    /// Export `mesh` to a `(json, binary_buffer)` pair.
    ///
    /// The binary buffer holds interleaved position data; the JSON references
    /// it via a buffer view and accessor.
    pub fn export(&self, mesh: &MeshData) -> (String, Vec<u8>) {
        // Build binary buffer: float32 positions
        let mut bin: Vec<u8> = Vec::new();
        for &v in &mesh.vertices {
            for &c in &v {
                bin.extend_from_slice(&(c as f32).to_le_bytes());
            }
        }
        // Face indices as uint32
        let indices_offset = bin.len();
        for &face in &mesh.faces {
            for &vi in &face {
                bin.extend_from_slice(&(vi as u32).to_le_bytes());
            }
        }

        let vertex_count = mesh.vertices.len();
        let index_count = mesh.faces.len() * 3;
        let positions_byte_len = vertex_count * 3 * 4;
        let indices_byte_len = index_count * 4;
        let total_byte_len = bin.len();

        // Compute bounding box for accessor min/max
        let (bb_min, bb_max) = mesh.bounding_box();

        let json = format!(
            r#"{{
  "asset": {{"version": "2.0", "generator": "OxiPhysics"}},
  "scene": 0,
  "scenes": [{{"name": "{}", "nodes": [0]}}],
  "nodes": [{{"mesh": 0}}],
  "meshes": [{{
    "name": "mesh0",
    "primitives": [{{
      "attributes": {{"POSITION": 0}},
      "indices": 1,
      "mode": 4
    }}]
  }}],
  "accessors": [
    {{
      "bufferView": 0,
      "byteOffset": 0,
      "componentType": 5126,
      "count": {vertex_count},
      "type": "VEC3",
      "min": [{}, {}, {}],
      "max": [{}, {}, {}]
    }},
    {{
      "bufferView": 1,
      "byteOffset": 0,
      "componentType": 5125,
      "count": {index_count},
      "type": "SCALAR"
    }}
  ],
  "bufferViews": [
    {{"buffer": 0, "byteOffset": 0, "byteLength": {positions_byte_len}, "target": 34962}},
    {{"buffer": 0, "byteOffset": {indices_offset}, "byteLength": {indices_byte_len}, "target": 34963}}
  ],
  "buffers": [{{"byteLength": {total_byte_len}}}]
}}"#,
            self.scene_name,
            bb_min[0] as f32,
            bb_min[1] as f32,
            bb_min[2] as f32,
            bb_max[0] as f32,
            bb_max[1] as f32,
            bb_max[2] as f32,
        );

        (json, bin)
    }
}

// ============================================================================
// MeshImporter
// ============================================================================

/// Mesh importer supporting OBJ, binary STL, ASCII STL, and PLY.
#[derive(Debug, Default)]
pub struct MeshImporter;

impl MeshImporter {
    /// Create a new `MeshImporter`.
    pub fn new() -> Self {
        Self
    }

    /// Parse a Wavefront OBJ string into a `MeshData`.
    ///
    /// Supports `v`, `vn`, `vt`, and `f` directives.
    pub fn parse_obj(&self, src: &str) -> MeshData {
        let mut mesh = MeshData::new();
        let mut obj_normals: Vec<[f64; 3]> = Vec::new();
        let mut obj_uvs: Vec<[f64; 2]> = Vec::new();

        for line in src.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("v ") {
                let vals: Vec<f64> = rest
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 3 {
                    mesh.add_vertex([vals[0], vals[1], vals[2]]);
                }
            } else if let Some(rest) = line.strip_prefix("vn ") {
                let vals: Vec<f64> = rest
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 3 {
                    obj_normals.push([vals[0], vals[1], vals[2]]);
                }
            } else if let Some(rest) = line.strip_prefix("vt ") {
                let vals: Vec<f64> = rest
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 2 {
                    obj_uvs.push([vals[0], vals[1]]);
                }
            } else if let Some(rest) = line.strip_prefix("f ") {
                let parts: Vec<&str> = rest.split_whitespace().collect();
                if parts.len() >= 3 {
                    // Parse each token as vi/vt/vn (1-based)
                    let parse_index = |tok: &str| -> usize {
                        let s = tok.split('/').next().unwrap_or("1");
                        s.parse::<usize>().unwrap_or(1).saturating_sub(1)
                    };
                    let i0 = parse_index(parts[0]);
                    let i1 = parse_index(parts[1]);
                    let i2 = parse_index(parts[2]);
                    mesh.add_face(i0, i1, i2);
                    // Fan triangulation for quads / ngons
                    for k in 3..parts.len() {
                        let ik = parse_index(parts[k]);
                        mesh.add_face(i0, parse_index(parts[k - 1]), ik);
                    }
                }
            }
        }

        // Transfer normals if counts match
        if !obj_normals.is_empty() && obj_normals.len() == mesh.vertices.len() {
            mesh.normals = obj_normals;
        }
        if !obj_uvs.is_empty() && obj_uvs.len() == mesh.vertices.len() {
            mesh.uvs = obj_uvs;
        }

        mesh
    }

    /// Parse a binary STL byte slice into a `MeshData`.
    ///
    /// Expects the standard 80-byte header + 4-byte count + 50 bytes per triangle.
    pub fn parse_stl_binary(&self, data: &[u8]) -> Option<MeshData> {
        if data.len() < 84 {
            return None;
        }
        let count = u32::from_le_bytes(data[80..84].try_into().ok()?) as usize;
        let expected_len = 84 + count * 50;
        if data.len() < expected_len {
            return None;
        }

        let mut mesh = MeshData::new();
        let mut offset = 84usize;
        for _ in 0..count {
            let read_f32 = |buf: &[u8], off: usize| -> f32 {
                f32::from_le_bytes(
                    buf[off..off + 4]
                        .try_into()
                        .expect("slice length must match"),
                )
            };
            let nx = read_f32(data, offset) as f64;
            let ny = read_f32(data, offset + 4) as f64;
            let nz = read_f32(data, offset + 8) as f64;
            let v0 = [
                read_f32(data, offset + 12) as f64,
                read_f32(data, offset + 16) as f64,
                read_f32(data, offset + 20) as f64,
            ];
            let v1 = [
                read_f32(data, offset + 24) as f64,
                read_f32(data, offset + 28) as f64,
                read_f32(data, offset + 32) as f64,
            ];
            let v2 = [
                read_f32(data, offset + 36) as f64,
                read_f32(data, offset + 40) as f64,
                read_f32(data, offset + 44) as f64,
            ];
            let base = mesh.vertices.len();
            mesh.add_vertex(v0);
            mesh.add_vertex(v1);
            mesh.add_vertex(v2);
            mesh.add_face(base, base + 1, base + 2);
            mesh.normals.push([nx, ny, nz]);
            mesh.normals.push([nx, ny, nz]);
            mesh.normals.push([nx, ny, nz]);
            offset += 50;
        }
        Some(mesh)
    }

    /// Parse an ASCII STL string into a `MeshData`.
    pub fn parse_stl_ascii(&self, src: &str) -> MeshData {
        let mut mesh = MeshData::new();
        let mut current_normal = [0.0f64; 3];
        let mut verts_in_loop: Vec<[f64; 3]> = Vec::new();

        for line in src.lines() {
            let line = line.trim();
            if line.starts_with("facet normal") {
                let vals: Vec<f64> = line
                    .split_whitespace()
                    .skip(2)
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 3 {
                    current_normal = [vals[0], vals[1], vals[2]];
                }
            } else if line.starts_with("vertex") {
                let vals: Vec<f64> = line
                    .split_whitespace()
                    .skip(1)
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 3 {
                    verts_in_loop.push([vals[0], vals[1], vals[2]]);
                }
            } else if line.starts_with("endfacet") {
                if verts_in_loop.len() == 3 {
                    let base = mesh.vertices.len();
                    for &v in &verts_in_loop {
                        mesh.add_vertex(v);
                        mesh.normals.push(current_normal);
                    }
                    mesh.add_face(base, base + 1, base + 2);
                }
                verts_in_loop.clear();
            }
        }
        mesh
    }

    /// Parse an ASCII PLY string into a `MeshData`.
    pub fn parse_ply_ascii(&self, src: &str) -> Option<MeshData> {
        let mut lines = src.lines();
        let mut vertex_count = 0usize;
        let mut face_count = 0usize;
        let mut has_normal = false;
        let mut has_uv = false;
        // Parse header
        for line in lines.by_ref() {
            let line = line.trim();
            if line == "end_header" {
                break;
            }
            if line.starts_with("element vertex") {
                vertex_count = line
                    .split_whitespace()
                    .last()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
            } else if line.starts_with("element face") {
                face_count = line
                    .split_whitespace()
                    .last()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
            } else if line.contains("property float nx") || line.contains("property float nz") {
                has_normal = true;
            } else if line.contains("property float s") || line.contains("property float t") {
                has_uv = true;
            }
        }

        let mut mesh = MeshData::new();

        // Read vertices
        for _ in 0..vertex_count {
            let line = lines.next()?.trim().to_string();
            let vals: Vec<f64> = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if vals.len() < 3 {
                return None;
            }
            mesh.add_vertex([vals[0], vals[1], vals[2]]);
            if has_normal && vals.len() >= 6 {
                mesh.normals.push([vals[3], vals[4], vals[5]]);
            }
            let uv_offset = if has_normal { 6 } else { 3 };
            if has_uv && vals.len() >= uv_offset + 2 {
                mesh.uvs.push([vals[uv_offset], vals[uv_offset + 1]]);
            }
        }

        // Read faces
        for _ in 0..face_count {
            let line = lines.next()?.trim().to_string();
            let vals: Vec<usize> = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if vals.is_empty() {
                continue;
            }
            let n = vals[0];
            if n >= 3 && vals.len() > n {
                mesh.add_face(vals[1], vals[2], vals[3]);
                for k in 3..n {
                    mesh.add_face(vals[1], vals[k], vals[k + 1]);
                }
            }
        }

        Some(mesh)
    }
}

// ============================================================================
// MeshSimplification
// ============================================================================

/// Mesh simplification via quadric error metrics (QEM) and edge collapse.
///
/// Implements the algorithm of Garland and Heckbert, "Surface Simplification Using
/// Quadric Error Metrics" (SIGGRAPH 1997).
#[derive(Debug, Default)]
pub struct MeshSimplification;

impl MeshSimplification {
    /// Create a new `MeshSimplification`.
    pub fn new() -> Self {
        Self
    }

    /// Simplify `mesh` towards `target_triangles` triangles by quadric-error edge
    /// collapses.
    ///
    /// Every vertex carries the quadric of the planes of its triangles (weighted by
    /// triangle area) plus one heavily weighted plane per boundary edge, so open
    /// boundaries keep their outline.  The cheapest edge is collapsed first.  The merged
    /// vertex goes to the point that minimizes the sum of the two quadrics (or, when that
    /// point is not well defined, to the best point of the edge), it takes over the sum of
    /// the two quadrics, and only the edges around it are costed again.
    ///
    /// A collapse is skipped when it would change the topology of the surface (pinch it,
    /// join two boundary pieces, cut off or delete a piece) or fold a neighbouring
    /// triangle over (its normal would turn by 90 degrees or more).  The number of
    /// connected pieces and the Euler characteristic of the mesh therefore never change.
    /// The simplification stops as soon as the mesh has at most `target_triangles`
    /// triangles *or* no allowed collapse remains; a target below what the surface
    /// permits is not reached (a closed piece stays at four triangles, a disc at one).  A
    /// collapse removes one triangle along a boundary edge and two along an interior
    /// edge, so the result can also end one triangle below the target.
    ///
    /// Connectivity comes from shared vertex indices, so vertices that are duplicated at
    /// the same position (as in an STL triangle soup) must be welded first; triangles
    /// that share no vertex are separate pieces that cannot be reduced.
    ///
    /// Triangles that use a vertex index outside the mesh, the same vertex twice or a
    /// vertex with a non-finite coordinate cannot be simplified and are dropped.  The
    /// result holds exactly the vertices its triangles use, in order of first use, and,
    /// when the input has one material ID per triangle, the material IDs of the
    /// surviving triangles; normals and UV coordinates are not carried over (see
    /// [`MeshData::compute_smooth_normals`]).
    ///
    /// A mesh that already has at most `target_triangles` triangles is returned
    /// unchanged.
    pub fn simplify(&self, mesh: &MeshData, target_triangles: usize) -> MeshData {
        if mesh.faces.len() <= target_triangles {
            return mesh.clone();
        }
        qem::simplify(mesh, target_triangles)
    }
}

/// Build a mesh from `faces` that holds only the vertices they reference.
///
/// Output vertices are numbered in order of first use, so a vertex shared by several
/// faces is stored exactly once. Faces that reference a vertex outside `vertices` are
/// dropped.
fn compact_mesh(vertices: &[[f64; 3]], faces: &[[usize; 3]]) -> MeshData {
    let mut mesh = MeshData::new();
    let mut remap = vec![usize::MAX; vertices.len()];
    for face in faces {
        if face.iter().any(|&vi| vi >= vertices.len()) {
            continue;
        }
        let mut new_face = [0usize; 3];
        for (slot, &vi) in new_face.iter_mut().zip(face.iter()) {
            if remap[vi] == usize::MAX {
                remap[vi] = mesh.add_vertex(vertices[vi]);
            }
            *slot = remap[vi];
        }
        mesh.add_face(new_face[0], new_face[1], new_face[2]);
    }
    mesh
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: simple tetrahedron MeshData
    fn tetra() -> MeshData {
        let mut m = MeshData::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_vertex([0.0, 0.0, 1.0]);
        m.add_face(0, 1, 2);
        m.add_face(0, 1, 3);
        m.add_face(0, 2, 3);
        m.add_face(1, 2, 3);
        m
    }

    // Helper: unit square (two triangles)
    fn quad() -> MeshData {
        let mut m = MeshData::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([1.0, 1.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_face(0, 1, 2);
        m.add_face(0, 2, 3);
        m
    }

    // -- MeshData tests -------------------------------------------------------

    #[test]
    fn test_mesh_add_vertex() {
        let mut m = MeshData::new();
        let idx = m.add_vertex([1.0, 2.0, 3.0]);
        assert_eq!(idx, 0);
        assert_eq!(m.vertices.len(), 1);
    }

    #[test]
    fn test_mesh_add_face() {
        let m = tetra();
        assert_eq!(m.num_triangles(), 4);
    }

    #[test]
    fn test_mesh_bounding_box() {
        let m = tetra();
        let (mn, mx) = m.bounding_box();
        assert!((mn[0] - 0.0).abs() < 1e-10);
        assert!((mx[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mesh_compute_flat_normals() {
        let mut m = quad();
        m.compute_flat_normals();
        assert_eq!(m.normals.len(), m.vertices.len());
        // All normals should point in Z direction
        for n in &m.normals {
            assert!((n[2].abs() - 1.0).abs() < 1e-5);
        }
    }

    #[test]
    fn test_mesh_compute_smooth_normals() {
        let mut m = quad();
        m.compute_smooth_normals();
        assert_eq!(m.normals.len(), m.vertices.len());
    }

    #[test]
    fn test_mesh_add_face_with_material() {
        let mut m = MeshData::new();
        m.add_vertex([0.0; 3]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_face_with_material(0, 1, 2, 3);
        assert_eq!(m.material_ids[0], 3);
    }

    // -- ObjExporter tests ----------------------------------------------------

    #[test]
    fn test_obj_export_basic() {
        let exporter = ObjExporter::new();
        let m = quad();
        let obj = exporter.export_obj(&m, None);
        assert!(obj.contains("v 0 0 0"));
        assert!(obj.contains("f 1 2 3"));
    }

    #[test]
    fn test_obj_export_with_mtllib() {
        let exporter = ObjExporter::new();
        let m = quad();
        let obj = exporter.export_obj(&m, Some("scene.mtl"));
        assert!(obj.contains("mtllib scene.mtl"));
    }

    #[test]
    fn test_obj_export_normals() {
        let mut exporter = ObjExporter::new();
        exporter.add_material(Material::new("mat0"));
        let mut m = quad();
        m.compute_flat_normals();
        let obj = exporter.export_obj(&m, None);
        assert!(obj.contains("vn "));
    }

    #[test]
    fn test_mtl_export() {
        let mut exporter = ObjExporter::new();
        exporter.add_material(Material::new("red"));
        let mtl = exporter.export_mtl();
        assert!(mtl.contains("newmtl red"));
        assert!(mtl.contains("Kd "));
    }

    #[test]
    fn test_obj_usemtl_directive() {
        let mut exporter = ObjExporter::new();
        exporter.add_material(Material::new("mat0"));
        let mut m = quad();
        m.material_ids = vec![0, 0];
        let obj = exporter.export_obj(&m, None);
        assert!(obj.contains("usemtl mat0"));
    }

    // -- StlExporter tests ----------------------------------------------------

    #[test]
    fn test_stl_binary_header() {
        let exp = StlExporter::new();
        let m = quad();
        let bytes = exp.export_binary(&m, "test");
        assert!(bytes.len() >= 84);
        // First 80 bytes are header
        assert_eq!(&bytes[..4], b"test");
    }

    #[test]
    fn test_stl_binary_triangle_count() {
        let exp = StlExporter::new();
        let m = quad();
        let bytes = exp.export_binary(&m, "q");
        let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap());
        assert_eq!(count, 2);
    }

    #[test]
    fn test_stl_binary_size() {
        let exp = StlExporter::new();
        let m = quad();
        let bytes = exp.export_binary(&m, "q");
        assert_eq!(bytes.len(), 84 + 2 * 50);
    }

    #[test]
    fn test_stl_ascii_contains_solid() {
        let exp = StlExporter::new();
        let m = quad();
        let s = exp.export_ascii(&m, "mymesh");
        assert!(s.starts_with("solid mymesh"));
        assert!(s.contains("endsolid mymesh"));
    }

    #[test]
    fn test_stl_ascii_facet_count() {
        let exp = StlExporter::new();
        let m = quad();
        let s = exp.export_ascii(&m, "q");
        let count = s.matches("facet normal").count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_stl_crc32_known() {
        // CRC32 of empty slice should be 0x00000000
        let crc = StlExporter::crc32(&[]);
        assert_eq!(crc, 0x0000_0000);
    }

    // -- PlyExporter tests ----------------------------------------------------

    #[test]
    fn test_ply_ascii_header() {
        let exp = PlyExporter::new(false);
        let m = quad();
        let bytes = exp.export_ascii_bytes(&m);
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("ply"));
        assert!(s.contains("format ascii"));
        assert!(s.contains("element vertex 4"));
        assert!(s.contains("element face 2"));
    }

    #[test]
    fn test_ply_binary_header() {
        let exp = PlyExporter::new(true);
        let m = quad();
        let bytes = exp.export_binary_bytes(&m);
        let header_end = bytes
            .windows(11)
            .position(|w| w == b"end_header\n")
            .unwrap()
            + 11;
        let header = String::from_utf8_lossy(&bytes[..header_end]);
        assert!(header.contains("format binary_little_endian"));
    }

    #[test]
    fn test_ply_export_dispatch() {
        let exp_ascii = PlyExporter::new(false);
        let exp_bin = PlyExporter::new(true);
        let m = quad();
        let a = exp_ascii.export(&m);
        let b = exp_bin.export(&m);
        assert!(!a.is_empty());
        assert!(!b.is_empty());
    }

    #[test]
    fn test_ply_ascii_with_normals() {
        let mut m = quad();
        m.compute_flat_normals();
        let exp = PlyExporter::new(false);
        let bytes = exp.export_ascii_bytes(&m);
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("property float nx"));
    }

    // -- GltfExporter tests ---------------------------------------------------

    #[test]
    fn test_gltf_json_structure() {
        let exp = GltfExporter::new("test_scene");
        let m = quad();
        let (json, _bin) = exp.export(&m);
        assert!(json.contains("\"version\": \"2.0\""));
        assert!(json.contains("test_scene"));
        assert!(json.contains("POSITION"));
    }

    #[test]
    fn test_gltf_binary_buffer_size() {
        let exp = GltfExporter::new("s");
        let m = quad(); // 4 verts, 2 faces
        let (_json, bin) = exp.export(&m);
        // 4 * 3 * 4 bytes for positions + 2 * 3 * 4 bytes for indices
        let expected = 4 * 3 * 4 + 2 * 3 * 4;
        assert_eq!(bin.len(), expected);
    }

    #[test]
    fn test_gltf_accessor_count() {
        let exp = GltfExporter::new("s");
        let m = quad();
        let (json, _) = exp.export(&m);
        assert!(json.contains("\"count\": 4")); // vertex count
        assert!(json.contains("\"count\": 6")); // index count
    }

    // -- MeshImporter OBJ tests -----------------------------------------------

    #[test]
    fn test_import_obj_basic() {
        let src = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let imp = MeshImporter::new();
        let m = imp.parse_obj(src);
        assert_eq!(m.vertices.len(), 3);
        assert_eq!(m.faces.len(), 1);
    }

    #[test]
    fn test_import_obj_normals() {
        let src = "v 0 0 0\nv 1 0 0\nv 0 1 0\nvn 0 0 1\nvn 0 0 1\nvn 0 0 1\nf 1//1 2//2 3//3\n";
        let imp = MeshImporter::new();
        let m = imp.parse_obj(src);
        assert_eq!(m.vertices.len(), 3);
        // Normals are only stored when count matches vertices
        assert_eq!(m.normals.len(), 3);
    }

    #[test]
    fn test_import_obj_quad_fan() {
        let src = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\n";
        let imp = MeshImporter::new();
        let m = imp.parse_obj(src);
        assert_eq!(m.faces.len(), 2); // quad -> 2 tris
    }

    #[test]
    fn test_import_stl_binary_roundtrip() {
        let exp = StlExporter::new();
        let m = quad();
        let bytes = exp.export_binary(&m, "test");
        let imp = MeshImporter::new();
        let m2 = imp.parse_stl_binary(&bytes).unwrap();
        assert_eq!(m2.faces.len(), 2);
    }

    #[test]
    fn test_import_stl_ascii_roundtrip() {
        let exp = StlExporter::new();
        let m = quad();
        let s = exp.export_ascii(&m, "q");
        let imp = MeshImporter::new();
        let m2 = imp.parse_stl_ascii(&s);
        assert_eq!(m2.faces.len(), 2);
    }

    #[test]
    fn test_import_stl_binary_too_short() {
        let imp = MeshImporter::new();
        let result = imp.parse_stl_binary(&[0u8; 10]);
        assert!(result.is_none());
    }

    #[test]
    fn test_import_ply_ascii_roundtrip() {
        let exp = PlyExporter::new(false);
        let m = quad();
        let bytes = exp.export_ascii_bytes(&m);
        let s = String::from_utf8(bytes).unwrap();
        let imp = MeshImporter::new();
        let m2 = imp.parse_ply_ascii(&s).unwrap();
        assert_eq!(m2.vertices.len(), 4);
        assert_eq!(m2.faces.len(), 2);
    }

    #[test]
    fn test_import_ply_ascii_with_normals() {
        let exp = PlyExporter::new(false);
        let mut m = quad();
        m.compute_flat_normals();
        let bytes = exp.export_ascii_bytes(&m);
        let s = String::from_utf8(bytes).unwrap();
        let imp = MeshImporter::new();
        let m2 = imp.parse_ply_ascii(&s).unwrap();
        assert_eq!(m2.normals.len(), 4);
    }

    // -- MeshSimplification tests ---------------------------------------------

    #[test]
    fn test_simplify_already_simple() {
        let s = MeshSimplification::new();
        let m = quad();
        let out = s.simplify(&m, 10);
        assert_eq!(out.faces.len(), m.faces.len());
    }

    #[test]
    fn test_simplify_reduces_triangles() {
        // The diagonal of a quad joins two boundary vertices and is never collapsed; a
        // boundary edge is, which removes one of the two triangles.
        let out = MeshSimplification::new().simplify(&quad(), 1);
        assert_eq!(out.faces.len(), 1);
        assert_eq!(out.vertices.len(), 3);
        assert_consistent(&out);
    }

    #[test]
    fn test_simplify_never_reduces_a_closed_surface_below_a_tetrahedron() {
        // Collapsing any edge of a tetrahedron would leave two coincident triangles.
        let out = MeshSimplification::new().simplify(&tetra(), 2);
        assert_eq!(out.faces.len(), 4);
        assert_eq!(out.vertices.len(), 4);
        assert_consistent(&out);
    }

    #[test]
    fn test_simplify_empty_mesh() {
        let s = MeshSimplification::new();
        let m = MeshData::new();
        let out = s.simplify(&m, 100);
        assert_eq!(out.faces.len(), 0);
    }

    #[test]
    fn test_simplify_single_face() {
        let s = MeshSimplification::new();
        let mut m = MeshData::new();
        m.add_vertex([0.0; 3]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_face(0, 1, 2);
        let out = s.simplify(&m, 1);
        assert_eq!(out.faces.len(), 1);
    }

    /// Square grid of `n x n` cells in the z = 0 plane, two triangles per cell.
    fn grid_mesh(n: usize) -> MeshData {
        let mut m = MeshData::new();
        for j in 0..=n {
            for i in 0..=n {
                m.add_vertex([i as f64, j as f64, 0.0]);
            }
        }
        let stride = n + 1;
        for j in 0..n {
            for i in 0..n {
                let a = j * stride + i;
                let b = a + 1;
                let c = a + stride;
                let d = c + 1;
                m.add_face(a, b, d);
                m.add_face(a, d, c);
            }
        }
        m
    }

    /// Assert the structural invariants every simplified mesh must satisfy.
    fn assert_consistent(out: &MeshData) {
        let mut referenced = std::collections::BTreeSet::new();
        for face in &out.faces {
            assert!(
                face.iter().all(|&vi| vi < out.vertices.len()),
                "face {face:?} indexes past {} vertices",
                out.vertices.len()
            );
            assert!(
                face[0] != face[1] && face[1] != face[2] && face[0] != face[2],
                "degenerate face {face:?}"
            );
            referenced.extend(face.iter().copied());
        }
        assert_eq!(
            out.vertices.len(),
            referenced.len(),
            "every output vertex must be referenced by a face"
        );
        let mut positions = std::collections::BTreeSet::new();
        for v in &out.vertices {
            assert!(
                positions.insert([v[0].to_bits(), v[1].to_bits(), v[2].to_bits()]),
                "duplicate output vertex position {v:?}"
            );
        }
    }

    #[test]
    fn test_compact_mesh_keeps_shared_vertex_positions() {
        // Faces [0,1,2] and [2,3,4] share vertex 2, the most recently inserted one.
        let vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let out = compact_mesh(&vertices, &[[0, 1, 2], [2, 3, 4]]);
        assert_eq!(out.vertices, vertices.to_vec());
        assert_eq!(out.faces, vec![[0, 1, 2], [2, 3, 4]]);
    }

    #[test]
    fn test_compact_mesh_renumbers_in_first_use_order() {
        let vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
        ];
        let out = compact_mesh(&vertices, &[[5, 1, 3], [3, 1, 0]]);
        assert_eq!(
            out.vertices,
            vec![
                [5.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]
            ]
        );
        assert_eq!(out.faces, vec![[0, 1, 2], [2, 1, 3]]);
    }

    #[test]
    fn test_compact_mesh_drops_faces_with_unknown_vertices() {
        let vertices = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let out = compact_mesh(&vertices, &[[0, 1, 2], [0, 1, 9]]);
        assert_eq!(out.faces, vec![[0, 1, 2]]);
        assert_eq!(out.vertices.len(), 3);
    }

    #[test]
    fn test_simplify_grid_output_is_consistent() {
        let s = MeshSimplification::new();
        let m = grid_mesh(4);
        assert_eq!(m.faces.len(), 32);
        assert_consistent(&m);
        let out = s.simplify(&m, 16);
        assert!(!out.faces.is_empty());
        assert!(out.faces.len() <= 16, "{} faces", out.faces.len());
        assert_consistent(&out);
        // A planar mesh stays planar and inside the original bounding box.
        for v in &out.vertices {
            assert!(v[2].abs() < 1e-12);
            assert!((-1e-12..=4.0 + 1e-12).contains(&v[0]), "{v:?}");
            assert!((-1e-12..=4.0 + 1e-12).contains(&v[1]), "{v:?}");
        }
    }

    #[test]
    fn test_simplify_triangle_strip_is_consistent() {
        // Every face of a strip reuses the vertex inserted by the previous face.
        let mut m = MeshData::new();
        for i in 0..8 {
            m.add_vertex([i as f64, (i % 2) as f64, 0.0]);
        }
        for i in 0..6usize {
            if i.is_multiple_of(2) {
                m.add_face(i, i + 1, i + 2);
            } else {
                m.add_face(i + 1, i, i + 2);
            }
        }
        let out = MeshSimplification::new().simplify(&m, 4);
        assert!(out.faces.len() <= 4);
        assert_consistent(&out);
    }

    #[test]
    fn test_simplify_is_deterministic() {
        let s = MeshSimplification::new();
        let m = grid_mesh(5);
        let a = s.simplify(&m, 20);
        let b = s.simplify(&m, 20);
        assert_eq!(a.faces, b.faces);
        assert_eq!(a.vertices, b.vertices);
    }

    #[test]
    fn test_simplify_tolerates_non_finite_vertices() {
        let mut m = grid_mesh(3);
        m.vertices[5] = [f64::NAN, 0.0, 0.0];
        let out = MeshSimplification::new().simplify(&m, 8);
        assert!(out.faces.len() <= 8);
        for face in &out.faces {
            assert!(face.iter().all(|&vi| vi < out.vertices.len()));
        }
        assert!(
            out.vertices.iter().flatten().all(|c| c.is_finite()),
            "the triangles around the NaN vertex are dropped"
        );
    }

    // -- Material tests -------------------------------------------------------

    #[test]
    fn test_material_default() {
        let mat = Material::default();
        assert_eq!(mat.name, "default");
        assert!((mat.opacity - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_material_new() {
        let mat = Material::new("chrome");
        assert_eq!(mat.name, "chrome");
    }

    // -- Geometry helper tests -----------------------------------------------

    #[test]
    fn test_triangle_normal_z() {
        let n = triangle_normal_f64([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((n[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3_f64_unit() {
        let n = normalize3_f64([3.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3_f64_zero() {
        let n = normalize3_f64([0.0; 3]);
        assert_eq!(n, [0.0; 3]);
    }
}
