//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::collections::HashMap;

/// A single triangle in an STL mesh, with a normal vector and three vertices.
pub struct StlTriangle {
    /// Surface normal vector (unit vector).
    pub normal: [f32; 3],
    /// First vertex.
    pub v0: [f32; 3],
    /// Second vertex.
    pub v1: [f32; 3],
    /// Third vertex.
    pub v2: [f32; 3],
}
/// An STL mesh consisting of named triangles.
pub struct StlMesh {
    /// The triangles that make up this mesh.
    pub triangles: Vec<StlTriangle>,
    /// The name of the solid (used in ASCII output).
    pub name: String,
}
impl StlMesh {
    /// Create a new empty STL mesh with the given name.
    pub fn new(name: &str) -> Self {
        StlMesh {
            triangles: Vec::new(),
            name: name.to_string(),
        }
    }
    /// Add a triangle defined by three vertices; the normal is computed automatically.
    pub fn add_triangle(&mut self, v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) {
        let normal = compute_normal(v0, v1, v2);
        self.triangles.push(StlTriangle { normal, v0, v1, v2 });
    }
    /// Serialize this mesh to binary STL bytes.
    ///
    /// Layout: 80-byte header + 4-byte little-endian triangle count +
    /// (12 floats + 2-byte attribute) per triangle = 50 bytes each.
    pub fn to_binary_bytes(&self) -> Vec<u8> {
        let n = self.triangles.len();
        let mut buf = Vec::with_capacity(84 + n * 50);
        let mut header = [0u8; 80];
        let name_bytes = self.name.as_bytes();
        let copy_len = name_bytes.len().min(80);
        header[..copy_len].copy_from_slice(&name_bytes[..copy_len]);
        buf.extend_from_slice(&header);
        buf.extend_from_slice(&(n as u32).to_le_bytes());
        for tri in &self.triangles {
            for &f in &tri.normal {
                buf.extend_from_slice(&f.to_le_bytes());
            }
            for &f in &tri.v0 {
                buf.extend_from_slice(&f.to_le_bytes());
            }
            for &f in &tri.v1 {
                buf.extend_from_slice(&f.to_le_bytes());
            }
            for &f in &tri.v2 {
                buf.extend_from_slice(&f.to_le_bytes());
            }
            buf.extend_from_slice(&0u16.to_le_bytes());
        }
        buf
    }
    /// Parse a binary STL byte slice into an `StlMesh`.
    pub fn from_binary_bytes(data: &[u8]) -> Result<Self, String> {
        if data.len() < 84 {
            return Err(format!(
                "Binary STL too short: {} bytes (need at least 84)",
                data.len()
            ));
        }
        let header_bytes = &data[..80];
        let name = String::from_utf8_lossy(header_bytes)
            .trim_end_matches('\0')
            .trim()
            .to_string();
        let n = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
        let expected_len = 84 + n * 50;
        if data.len() < expected_len {
            return Err(format!(
                "Binary STL too short for {} triangles: got {} bytes, need {}",
                n,
                data.len(),
                expected_len
            ));
        }
        let mut triangles = Vec::with_capacity(n);
        for i in 0..n {
            let base = 84 + i * 50;
            let normal = read_f32x3(data, base)?;
            let v0 = read_f32x3(data, base + 12)?;
            let v1 = read_f32x3(data, base + 24)?;
            let v2 = read_f32x3(data, base + 36)?;
            triangles.push(StlTriangle { normal, v0, v1, v2 });
        }
        Ok(StlMesh { triangles, name })
    }
    /// Serialize this mesh to ASCII STL text.
    pub fn to_ascii(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("solid {}\n", self.name));
        for tri in &self.triangles {
            s.push_str(&format!(
                "  facet normal {:.6e} {:.6e} {:.6e}\n",
                tri.normal[0], tri.normal[1], tri.normal[2]
            ));
            s.push_str("    outer loop\n");
            s.push_str(&format!(
                "      vertex {:.6e} {:.6e} {:.6e}\n",
                tri.v0[0], tri.v0[1], tri.v0[2]
            ));
            s.push_str(&format!(
                "      vertex {:.6e} {:.6e} {:.6e}\n",
                tri.v1[0], tri.v1[1], tri.v1[2]
            ));
            s.push_str(&format!(
                "      vertex {:.6e} {:.6e} {:.6e}\n",
                tri.v2[0], tri.v2[1], tri.v2[2]
            ));
            s.push_str("    endloop\n");
            s.push_str("  endfacet\n");
        }
        s.push_str(&format!("endsolid {}\n", self.name));
        s
    }
    /// Parse an ASCII STL string into an `StlMesh`.
    pub fn from_ascii(s: &str) -> Result<Self, String> {
        let mut lines = s.lines().peekable();
        let first = lines
            .next()
            .ok_or_else(|| "Empty STL string".to_string())?
            .trim();
        let name = if let Some(rest) = first.strip_prefix("solid") {
            rest.trim().to_string()
        } else {
            return Err(format!("Expected 'solid', got: {first}"));
        };
        let mut triangles = Vec::new();
        while let Some(l) = lines.next() {
            let line = l.trim().to_string();
            if line.starts_with("endsolid") {
                break;
            }
            if line.is_empty() {
                continue;
            }
            if !line.starts_with("facet normal") {
                return Err(format!("Expected 'facet normal', got: {line}"));
            }
            let normal = parse_vec3_from_line(&line, "facet normal")?;
            let outer = lines
                .next()
                .ok_or_else(|| "Unexpected EOF after facet normal".to_string())?
                .trim()
                .to_string();
            if !outer.starts_with("outer loop") {
                return Err(format!("Expected 'outer loop', got: {outer}"));
            }
            let v0 = parse_vertex_line(
                lines
                    .next()
                    .ok_or_else(|| "Unexpected EOF reading vertex 0".to_string())?
                    .trim(),
            )?;
            let v1 = parse_vertex_line(
                lines
                    .next()
                    .ok_or_else(|| "Unexpected EOF reading vertex 1".to_string())?
                    .trim(),
            )?;
            let v2 = parse_vertex_line(
                lines
                    .next()
                    .ok_or_else(|| "Unexpected EOF reading vertex 2".to_string())?
                    .trim(),
            )?;
            let _endloop = lines
                .next()
                .ok_or_else(|| "Unexpected EOF reading endloop".to_string())?;
            let _endfacet = lines
                .next()
                .ok_or_else(|| "Unexpected EOF reading endfacet".to_string())?;
            triangles.push(StlTriangle { normal, v0, v1, v2 });
        }
        Ok(StlMesh { triangles, name })
    }
    /// Compute the total surface area of this mesh.
    pub fn surface_area(&self) -> f32 {
        self.triangles.iter().map(triangle_area).sum()
    }
    /// Compute the axis-aligned bounding box as `(min, max)`.
    ///
    /// Returns `([0.0;3], [0.0;3])` for an empty mesh.
    pub fn bounding_box(&self) -> ([f32; 3], [f32; 3]) {
        if self.triangles.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = [f32::MAX; 3];
        let mut hi = [f32::MIN; 3];
        for tri in &self.triangles {
            for v in [&tri.v0, &tri.v1, &tri.v2] {
                for k in 0..3 {
                    if v[k] < lo[k] {
                        lo[k] = v[k];
                    }
                    if v[k] > hi[k] {
                        hi[k] = v[k];
                    }
                }
            }
        }
        (lo, hi)
    }
    /// Return `true` if every edge is shared by exactly two triangles (watertight).
    pub fn is_watertight(&self) -> bool {
        let mut edge_counts: HashMap<EdgeKey, usize> = HashMap::new();
        for tri in &self.triangles {
            for &(a, b) in &[(tri.v0, tri.v1), (tri.v1, tri.v2), (tri.v2, tri.v0)] {
                let key = EdgeKey::new(a, b);
                *edge_counts.entry(key).or_insert(0) += 1;
            }
        }
        edge_counts.values().all(|&c| c == 2)
    }
}
/// Result of an STL mesh validation pass.
#[derive(Debug, Clone, Default)]
pub struct StlValidationReport {
    /// Number of boundary (non-manifold) edges.
    pub boundary_edge_count: usize,
    /// Whether the mesh appears watertight (no boundary edges).
    pub is_watertight: bool,
    /// Whether the mesh appears manifold (each edge shared by exactly 2 triangles).
    pub is_manifold: bool,
    /// Number of degenerate triangles (zero-area).
    pub degenerate_count: usize,
}
/// A simple indexed triangle mesh produced from the STL-to-mesh pipeline.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// Unique vertex positions.
    pub positions: Vec<[f32; 3]>,
    /// Per-vertex normals (averaged from all adjacent triangles).
    pub normals: Vec<[f32; 3]>,
    /// Triangle indices (3 indices per triangle).
    pub indices: Vec<[usize; 3]>,
}
/// Color encoded in the 2-byte attribute field of a binary STL triangle
/// using the VisCAM/SolidView convention:
/// bit 15 = valid flag, bits 14–10 = R (5-bit), 9–5 = G, 4–0 = B.
#[derive(Debug, Clone, Copy, Default)]
pub struct StlColor {
    /// Red channel (0–31).
    pub r: u8,
    /// Green channel (0–31).
    pub g: u8,
    /// Blue channel (0–31).
    pub b: u8,
}
impl StlColor {
    /// Create a color from 5-bit R/G/B channels.
    pub fn new(r: u8, g: u8, b: u8) -> Self {
        Self {
            r: r & 0x1F,
            g: g & 0x1F,
            b: b & 0x1F,
        }
    }
    /// Encode into a 2-byte attribute word (VisCAM format, bit 15 set).
    pub fn encode(&self) -> u16 {
        0x8000 | ((self.r as u16) << 10) | ((self.g as u16) << 5) | (self.b as u16)
    }
    /// Decode from a 2-byte attribute word.  Returns `None` if the valid bit
    /// (bit 15) is not set.
    pub fn decode(attr: u16) -> Option<Self> {
        if attr & 0x8000 == 0 {
            return None;
        }
        Some(Self {
            r: ((attr >> 10) & 0x1F) as u8,
            g: ((attr >> 5) & 0x1F) as u8,
            b: (attr & 0x1F) as u8,
        })
    }
    /// Convert to normalized (0.0–1.0) RGB triple.
    pub fn to_normalized(&self) -> [f32; 3] {
        [
            self.r as f32 / 31.0,
            self.g as f32 / 31.0,
            self.b as f32 / 31.0,
        ]
    }
}
/// Quality metric report for an STL mesh.
#[derive(Debug, Clone)]
pub struct StlQualityMetrics {
    /// Aspect ratio of each triangle (longest edge / shortest altitude).
    pub max_aspect_ratio: f32,
    /// Minimum triangle area.
    pub min_area: f32,
    /// Maximum triangle area.
    pub max_area: f32,
    /// Fraction of triangles that are degenerate (zero area).
    pub degenerate_fraction: f32,
    /// Whether any NaN/Inf coordinates exist.
    pub has_bad_geometry: bool,
}
/// Validation result for an STL mesh.
#[derive(Debug, Clone)]
pub struct StlValidation {
    /// Number of degenerate (zero-area) triangles.
    pub degenerate_count: usize,
    /// Number of non-manifold edges (shared by != 2 triangles).
    pub non_manifold_edges: usize,
    /// Whether the mesh is watertight.
    pub is_watertight: bool,
    /// Number of triangles with NaN or Inf vertices.
    pub nan_inf_count: usize,
    /// Whether all normals are unit-length.
    pub normals_valid: bool,
}
/// A canonicalised (sorted) directed edge key using bit-exact f32 representation.
#[derive(PartialEq, Eq, Hash)]
pub(super) struct EdgeKey {
    pub(super) a: [u32; 3],
    pub(super) b: [u32; 3],
}
impl EdgeKey {
    pub(super) fn new(p: [f32; 3], q: [f32; 3]) -> Self {
        let pa = p.map(f32::to_bits);
        let qa = q.map(f32::to_bits);
        if pa <= qa {
            EdgeKey { a: pa, b: qa }
        } else {
            EdgeKey { a: qa, b: pa }
        }
    }
}
/// Statistics for an STL mesh.
#[derive(Debug, Clone)]
pub struct StlStatistics {
    /// Number of triangles.
    pub triangle_count: usize,
    /// Total surface area.
    pub surface_area: f32,
    /// Bounding box min corner.
    pub bb_min: [f32; 3],
    /// Bounding box max corner.
    pub bb_max: [f32; 3],
    /// Bounding box size.
    pub bb_size: [f32; 3],
    /// Average triangle area.
    pub avg_triangle_area: f32,
    /// Minimum triangle area.
    pub min_triangle_area: f32,
    /// Maximum triangle area.
    pub max_triangle_area: f32,
    /// Average edge length.
    pub avg_edge_length: f32,
    /// Number of unique vertices (approximate, by bit-exact matching).
    pub approx_unique_vertices: usize,
}
/// Vertex welding result: unique vertices and remapped face indices.
#[derive(Debug, Clone)]
pub struct WeldedMesh {
    /// Unique vertex positions.
    pub vertices: Vec<[f32; 3]>,
    /// Triangle indices into `vertices`.
    pub triangles: Vec<[usize; 3]>,
}
