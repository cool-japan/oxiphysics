//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_strip, clippy::should_implement_trait)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

use crate::{Error, Result};
use oxiphysics_core::math::Vec3;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// A simple scene graph for OBJ scenes.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ObjScene {
    /// All nodes in the scene.
    pub nodes: Vec<ObjSceneNode>,
    /// All meshes in the scene.
    pub meshes: Vec<ObjMesh>,
    /// All materials in the scene.
    pub materials: Vec<ObjMaterial>,
}
#[allow(dead_code)]
impl ObjScene {
    /// Create an empty scene.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a mesh to the scene and return its index.
    pub fn add_mesh(&mut self, mesh: ObjMesh) -> usize {
        let idx = self.meshes.len();
        self.meshes.push(mesh);
        idx
    }
    /// Add a material to the scene and return its index.
    pub fn add_material(&mut self, mat: ObjMaterial) -> usize {
        let idx = self.materials.len();
        self.materials.push(mat);
        idx
    }
    /// Add a node to the scene and return its index.
    pub fn add_node(
        &mut self,
        name: &str,
        mesh_index: Option<usize>,
        transform: MeshTransform,
    ) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(ObjSceneNode {
            name: name.to_string(),
            transform,
            mesh_index,
            children: Vec::new(),
        });
        idx
    }
    /// Set a node as a child of another node.
    pub fn add_child(&mut self, parent_idx: usize, child_idx: usize) {
        if parent_idx < self.nodes.len() {
            self.nodes[parent_idx].children.push(child_idx);
        }
    }
    /// Flatten the scene hierarchy into a single merged `ObjMesh`.
    ///
    /// Each node's mesh is transformed by the node's local transform.
    pub fn flatten(&self) -> ObjMesh {
        let mut result = ObjMesh::default();
        for node in &self.nodes {
            if let Some(mi) = node.mesh_index
                && let Some(mesh) = self.meshes.get(mi)
            {
                let inst = MeshInstance {
                    name: node.name.clone(),
                    transform: node.transform.clone(),
                };
                let xformed = instantiate_mesh(mesh, &inst);
                result = merge_obj_meshes(&result, &xformed);
            }
        }
        result
    }
    /// Total number of vertices across all meshes.
    pub fn total_vertices(&self) -> usize {
        self.meshes.iter().map(|m| m.vertices.len()).sum()
    }
    /// Total number of faces across all meshes.
    pub fn total_faces(&self) -> usize {
        self.meshes.iter().map(|m| m.faces.len()).sum()
    }
}
/// A simple OBJ curve (polyline or rational B-spline stub).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjCurve {
    /// Curve name.
    pub name: String,
    /// Degree of the curve (1 = polyline).
    pub degree: usize,
    /// Control point indices into the vertex array.
    pub control_points: Vec<usize>,
    /// Knot vector (empty for polylines).
    pub knots: Vec<f64>,
}
/// Reader for Wavefront OBJ files.
#[allow(dead_code)]
pub struct ObjReader;
#[allow(dead_code)]
impl ObjReader {
    /// Parse an OBJ string into an [`ObjMesh`].
    pub fn from_str(data: &str) -> Result<ObjMesh> {
        let mut mesh = ObjMesh::default();
        let mut current_group_name: Option<String> = None;
        let mut current_group_start: usize = 0;
        let mut current_smoothing_group: u32 = 0;
        let mut current_material: Option<String> = None;
        for raw in data.lines() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if line.starts_with("g ") || line.starts_with("o ") {
                if let Some(ref name) = current_group_name {
                    let count = mesh.faces.len() - current_group_start;
                    if count > 0 {
                        mesh.groups.push(ObjGroup {
                            name: name.clone(),
                            face_start: current_group_start,
                            face_count: count,
                        });
                    }
                }
                let name = line[2..].trim().to_string();
                current_group_name = Some(name);
                current_group_start = mesh.faces.len();
            } else if line.starts_with("s ") {
                let val = line[2..].trim();
                current_smoothing_group = if val == "off" || val == "0" {
                    0
                } else {
                    val.parse::<u32>().unwrap_or(0)
                };
            } else if line.starts_with("usemtl ") {
                current_material = Some(line[7..].trim().to_string());
            } else {
                Self::parse_line_extended(
                    line,
                    &mut mesh,
                    current_smoothing_group,
                    &current_material,
                )?;
            }
        }
        if let Some(ref name) = current_group_name {
            let count = mesh.faces.len() - current_group_start;
            if count > 0 {
                mesh.groups.push(ObjGroup {
                    name: name.clone(),
                    face_start: current_group_start,
                    face_count: count,
                });
            }
        }
        Ok(mesh)
    }
    fn parse_line_extended(
        line: &str,
        mesh: &mut ObjMesh,
        smoothing_group: u32,
        material: &Option<String>,
    ) -> Result<()> {
        if line.starts_with("vn ") {
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() >= 4 {
                let x = p[1]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                let y = p[2]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                let z = p[3]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                mesh.normals.push([x, y, z]);
            }
        } else if line.starts_with("vt ") {
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() >= 3 {
                let u = p[1]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                let v = p[2]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                mesh.uvs.push([u, v]);
            }
        } else if line.starts_with("v ") {
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() >= 4 {
                let x = p[1]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                let y = p[2]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                let z = p[3]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                mesh.vertices.push([x, y, z]);
            }
        } else if line.starts_with("f ") {
            let p: Vec<&str> = line.split_whitespace().collect();
            let mut vis = Vec::new();
            let mut vts = Vec::new();
            let mut vns = Vec::new();
            let mut has_vt = false;
            let mut has_vn = false;
            for tok in &p[1..] {
                let parts: Vec<&str> = tok.split('/').collect();
                let vi = parts[0]
                    .parse::<usize>()
                    .map_err(|e| Error::Parse(e.to_string()))?
                    - 1;
                vis.push(vi);
                if parts.len() >= 2 && !parts[1].is_empty() {
                    has_vt = true;
                    let vt = parts[1]
                        .parse::<usize>()
                        .map_err(|e| Error::Parse(e.to_string()))?
                        - 1;
                    vts.push(vt);
                }
                if parts.len() >= 3 && !parts[2].is_empty() {
                    has_vn = true;
                    let vn = parts[2]
                        .parse::<usize>()
                        .map_err(|e| Error::Parse(e.to_string()))?
                        - 1;
                    vns.push(vn);
                }
            }
            mesh.faces.push(ObjFace {
                vertex_indices: vis,
                normal_indices: if has_vn { Some(vns) } else { None },
                uv_indices: if has_vt { Some(vts) } else { None },
                smoothing_group,
                material: material.clone(),
            });
        }
        Ok(())
    }
    /// Read an OBJ file and parse into an [`ObjMesh`].
    pub fn from_file(path: &str) -> Result<ObjMesh> {
        let file = File::open(Path::new(path))?;
        let reader = BufReader::new(file);
        let mut data = String::new();
        for raw in reader.lines() {
            let raw = raw?;
            data.push_str(&raw);
            data.push('\n');
        }
        Self::from_str(&data)
    }
    /// Read vertices and triangle faces from an OBJ file (legacy API).
    ///
    /// Returns `(vertices, faces)`. Face indices are converted to 0-based.
    pub fn read(path: &str) -> Result<(Vec<Vec3>, Vec<[usize; 3]>)> {
        let mesh = Self::from_file(path)?;
        let vertices: Vec<Vec3> = mesh
            .vertices
            .iter()
            .map(|v| Vec3::new(v[0], v[1], v[2]))
            .collect();
        let faces: Vec<[usize; 3]> = mesh
            .faces
            .iter()
            .filter(|f| f.vertex_indices.len() >= 3)
            .map(|f| {
                [
                    f.vertex_indices[0],
                    f.vertex_indices[1],
                    f.vertex_indices[2],
                ]
            })
            .collect();
        Ok((vertices, faces))
    }
}
/// A simple OBJ material definition.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjMaterial {
    /// Material name.
    pub name: String,
    /// Diffuse colour (Kd).
    pub kd: [f64; 3],
    /// Specular colour (Ks).
    pub ks: [f64; 3],
    /// Shininess exponent (Ns).
    pub ns: f64,
    /// Ambient colour (Ka).
    pub ka: [f64; 3],
    /// Dissolve / transparency (d, 1.0 = opaque).
    pub dissolve: f64,
    /// Diffuse texture map filename.
    pub map_kd: Option<String>,
}
#[allow(dead_code)]
impl ObjMaterial {
    /// Create a basic material with just a name and diffuse color.
    pub fn basic(name: &str, kd: [f64; 3]) -> Self {
        Self {
            name: name.to_string(),
            kd,
            ks: [0.0; 3],
            ns: 1.0,
            ka: [0.0; 3],
            dissolve: 1.0,
            map_kd: None,
        }
    }
}
/// A level-of-detail collection: multiple meshes at decreasing resolution.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ObjLod {
    /// LOD levels, index 0 = highest detail.
    pub levels: Vec<ObjMesh>,
    /// Distance thresholds at which to switch to the next LOD level.
    /// `thresholds[i]` is the max distance for `levels[i]`.
    pub thresholds: Vec<f64>,
}
#[allow(dead_code)]
impl ObjLod {
    /// Create an empty LOD set.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a LOD level with the given switch distance.
    pub fn push(&mut self, mesh: ObjMesh, threshold: f64) {
        self.levels.push(mesh);
        self.thresholds.push(threshold);
    }
    /// Select the appropriate LOD level for a viewing distance.
    pub fn select(&self, distance: f64) -> Option<&ObjMesh> {
        for (i, &t) in self.thresholds.iter().enumerate() {
            if distance <= t {
                return self.levels.get(i);
            }
        }
        self.levels.last()
    }
    /// Number of LOD levels.
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }
    /// Simplify a mesh to approximately `target_faces` triangles using
    /// Quadric Error Metric (QEM) edge-collapse.
    ///
    /// Implements the Garland–Heckbert (1997) algorithm:
    /// 1. Compute per-face error quadrics from plane equations.
    /// 2. Sum quadrics at each vertex.
    /// 3. For every edge, compute the optimal collapse target position and
    ///    its associated cost.
    /// 4. Greedily collapse the minimum-cost edge until `target_faces` is
    ///    reached, using a lazy-deletion min-heap.
    pub fn decimate(mesh: &ObjMesh, target_faces: usize) -> ObjMesh {
        if mesh.faces.len() <= target_faces {
            return mesh.clone();
        }
        decimate_qem_impl(mesh, target_faces)
    }
}
/// A transform applied when instancing a mesh.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MeshTransform {
    /// Translation vector.
    pub translation: [f64; 3],
    /// Uniform scale factor.
    pub scale: f64,
    /// Rotation axis (unit vector).
    pub axis: [f64; 3],
    /// Rotation angle in radians.
    pub angle: f64,
}
#[allow(dead_code)]
impl MeshTransform {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            translation: [0.0; 3],
            scale: 1.0,
            axis: [0.0, 0.0, 1.0],
            angle: 0.0,
        }
    }
    /// Translation-only transform.
    pub fn from_translation(tx: f64, ty: f64, tz: f64) -> Self {
        Self {
            translation: [tx, ty, tz],
            scale: 1.0,
            axis: [0.0, 0.0, 1.0],
            angle: 0.0,
        }
    }
    /// Apply this transform to a single point.
    ///
    /// Applies scale, then rotation (Rodrigues formula), then translation.
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let s = [p[0] * self.scale, p[1] * self.scale, p[2] * self.scale];
        let (ax, ay, az) = (self.axis[0], self.axis[1], self.axis[2]);
        let (sin_a, cos_a) = (self.angle.sin(), self.angle.cos());
        let dot = ax * s[0] + ay * s[1] + az * s[2];
        let cross = [
            ay * s[2] - az * s[1],
            az * s[0] - ax * s[2],
            ax * s[1] - ay * s[0],
        ];
        let rx = s[0] * cos_a + cross[0] * sin_a + ax * dot * (1.0 - cos_a);
        let ry = s[1] * cos_a + cross[1] * sin_a + ay * dot * (1.0 - cos_a);
        let rz = s[2] * cos_a + cross[2] * sin_a + az * dot * (1.0 - cos_a);
        [
            rx + self.translation[0],
            ry + self.translation[1],
            rz + self.translation[2],
        ]
    }
}
/// Statistics computed from an `ObjMesh`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjMeshStats {
    /// Number of vertices.
    pub vertex_count: usize,
    /// Number of faces.
    pub face_count: usize,
    /// Number of triangles (after fan-triangulation).
    pub triangle_count: usize,
    /// Number of unique materials.
    pub material_count: usize,
    /// Number of groups.
    pub group_count: usize,
    /// Number of faces with normals.
    pub faces_with_normals: usize,
    /// Number of faces with UVs.
    pub faces_with_uvs: usize,
    /// Surface area (sum of triangle areas).
    pub surface_area: f64,
    /// Axis-aligned bounding box (min, max).
    pub bbox: Option<([f64; 3], [f64; 3])>,
}
/// A Wavefront OBJ mesh with vertices, normals, UVs, and faces.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ObjMesh {
    /// 3D vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Per-vertex normals.
    pub normals: Vec<[f64; 3]>,
    /// UV (texture) coordinates.
    pub uvs: Vec<[f64; 2]>,
    /// Face list.
    pub faces: Vec<ObjFace>,
    /// Named groups.
    pub groups: Vec<ObjGroup>,
}
#[allow(dead_code)]
impl ObjMesh {
    /// Extract a flat triangle soup from the mesh.
    ///
    /// Each triangle is `[[v0\], [v1], [v2]]` in world-space coordinates.
    /// Quad faces are triangulated as (0,1,2) and (0,2,3).
    pub fn to_triangle_soup(&self) -> Vec<[[f64; 3]; 3]> {
        let mut soup = Vec::new();
        for face in &self.faces {
            let verts = &face.vertex_indices;
            if verts.len() < 3 {
                continue;
            }
            for i in 1..(verts.len() - 1) {
                let v0 = self.vertices[verts[0]];
                let v1 = self.vertices[verts[i]];
                let v2 = self.vertices[verts[i + 1]];
                soup.push([v0, v1, v2]);
            }
        }
        soup
    }
    /// Return faces belonging to a specific group.
    pub fn faces_in_group(&self, group_name: &str) -> Vec<&ObjFace> {
        if let Some(group) = self.groups.iter().find(|g| g.name == group_name) {
            let end = group.face_start + group.face_count;
            let end = end.min(self.faces.len());
            self.faces[group.face_start..end].iter().collect()
        } else {
            Vec::new()
        }
    }
    /// Return faces belonging to a specific smoothing group.
    pub fn faces_in_smoothing_group(&self, sg: u32) -> Vec<&ObjFace> {
        self.faces
            .iter()
            .filter(|f| f.smoothing_group == sg)
            .collect()
    }
    /// Return faces with a specific material.
    pub fn faces_with_material(&self, mat_name: &str) -> Vec<&ObjFace> {
        self.faces
            .iter()
            .filter(|f| f.material.as_deref() == Some(mat_name))
            .collect()
    }
    /// Total number of triangles after fan-triangulation.
    pub fn triangle_count(&self) -> usize {
        self.faces
            .iter()
            .map(|f| {
                if f.vertex_indices.len() >= 3 {
                    f.vertex_indices.len() - 2
                } else {
                    0
                }
            })
            .sum()
    }
    /// Compute flat face normal for a triangular face.
    pub fn face_normal(&self, face_idx: usize) -> Option<[f64; 3]> {
        let face = self.faces.get(face_idx)?;
        if face.vertex_indices.len() < 3 {
            return None;
        }
        let v0 = self.vertices[face.vertex_indices[0]];
        let v1 = self.vertices[face.vertex_indices[1]];
        let v2 = self.vertices[face.vertex_indices[2]];
        let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len < 1e-30 {
            return None;
        }
        Some([n[0] / len, n[1] / len, n[2] / len])
    }
    /// Compute axis-aligned bounding box of all vertices.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.vertices.is_empty() {
            return None;
        }
        let mut min = self.vertices[0];
        let mut max = self.vertices[0];
        for v in &self.vertices[1..] {
            for k in 0..3 {
                if v[k] < min[k] {
                    min[k] = v[k];
                }
                if v[k] > max[k] {
                    max[k] = v[k];
                }
            }
        }
        Some((min, max))
    }
}
/// An instance of an `ObjMesh` with a transform and optional name override.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MeshInstance {
    /// Name of this instance.
    pub name: String,
    /// Transform to apply.
    pub transform: MeshTransform,
}
/// A named group/object within an OBJ file.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjGroup {
    /// Group or object name.
    pub name: String,
    /// Index of the first face belonging to this group (into the mesh face list).
    pub face_start: usize,
    /// Number of faces in this group.
    pub face_count: usize,
}
/// A single face in a Wavefront OBJ mesh.
///
/// All index arrays are optional to support the various OBJ face formats:
/// `f v`, `f v/vt`, `f v//vn`, `f v/vt/vn`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjFace {
    /// Vertex indices (0-based).
    pub vertex_indices: Vec<usize>,
    /// Normal indices (0-based), if present.
    pub normal_indices: Option<Vec<usize>>,
    /// UV/texture-coordinate indices (0-based), if present.
    pub uv_indices: Option<Vec<usize>>,
    /// Smoothing group this face belongs to (0 = off).
    pub smoothing_group: u32,
    /// Material name assigned to this face.
    pub material: Option<String>,
}
/// Writer for `ObjMesh` structures.
#[allow(dead_code)]
pub struct ObjWriter;
#[allow(dead_code)]
impl ObjWriter {
    /// Serialize an `ObjMesh` to a Wavefront OBJ string.
    pub fn write(mesh: &ObjMesh) -> String {
        Self::write_with_groups(mesh, false)
    }
    /// Serialize an `ObjMesh` with optional group headers.
    pub fn write_with_groups(mesh: &ObjMesh, emit_groups: bool) -> String {
        let mut s = String::from("# OxiPhysics OBJ export\n");
        for v in &mesh.vertices {
            s.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
        }
        for vn in &mesh.normals {
            s.push_str(&format!("vn {} {} {}\n", vn[0], vn[1], vn[2]));
        }
        for vt in &mesh.uvs {
            s.push_str(&format!("vt {} {}\n", vt[0], vt[1]));
        }
        let mut current_group: Option<&str> = None;
        let mut current_material: Option<&str> = None;
        let mut current_sg: u32 = 0;
        for (fi, face) in mesh.faces.iter().enumerate() {
            if emit_groups {
                for group in &mesh.groups {
                    if fi == group.face_start && current_group != Some(&group.name) {
                        s.push_str(&format!("g {}\n", group.name));
                        current_group = Some(&group.name);
                    }
                }
            }
            if let Some(ref mat) = face.material
                && current_material != Some(mat.as_str())
            {
                s.push_str(&format!("usemtl {}\n", mat));
                current_material = Some(mat);
            }
            if face.smoothing_group != current_sg {
                current_sg = face.smoothing_group;
                if current_sg == 0 {
                    s.push_str("s off\n");
                } else {
                    s.push_str(&format!("s {}\n", current_sg));
                }
            }
            s.push('f');
            for i in 0..face.vertex_indices.len() {
                let vi = face.vertex_indices[i] + 1;
                let vt_idx = face.uv_indices.as_ref().map(|uvs| uvs[i] + 1);
                let vn_idx = face.normal_indices.as_ref().map(|ns| ns[i] + 1);
                let token = match (vt_idx, vn_idx) {
                    (Some(vt), Some(vn)) => format!(" {}/{}/{}", vi, vt, vn),
                    (None, Some(vn)) => format!(" {}//{}", vi, vn),
                    (Some(vt), None) => format!(" {}/{}", vi, vt),
                    (None, None) => format!(" {}", vi),
                };
                s.push_str(&token);
            }
            s.push('\n');
        }
        s
    }
    /// Write an `ObjMesh` directly to a file at `path`.
    pub fn write_to_file(path: &str, mesh: &ObjMesh) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        write!(w, "{}", Self::write(mesh))?;
        w.flush()?;
        Ok(())
    }
    /// Write a triangle mesh (legacy API -- vertices + triangle index arrays).
    ///
    /// Optionally includes per-vertex normals.
    pub fn write_legacy(
        path: &str,
        vertices: &[Vec3],
        triangles: &[[usize; 3]],
        normals: Option<&[Vec3]>,
    ) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "# OxiPhysics OBJ export")?;
        for v in vertices {
            writeln!(w, "v {} {} {}", v.x, v.y, v.z)?;
        }
        if let Some(norms) = normals {
            for n in norms {
                writeln!(w, "vn {} {} {}", n.x, n.y, n.z)?;
            }
            for t in triangles {
                writeln!(
                    w,
                    "f {}//{} {}//{} {}//{}",
                    t[0] + 1,
                    t[0] + 1,
                    t[1] + 1,
                    t[1] + 1,
                    t[2] + 1,
                    t[2] + 1,
                )?;
            }
        } else {
            for t in triangles {
                writeln!(w, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1)?;
            }
        }
        w.flush()?;
        Ok(())
    }
    /// Write a mesh with texture coordinates (legacy API).
    pub fn write_with_uvs(
        path: &str,
        vertices: &[Vec3],
        uvs: &[[f64; 2]],
        triangles: &[[usize; 3]],
    ) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "# OxiPhysics OBJ export")?;
        for v in vertices {
            writeln!(w, "v {} {} {}", v.x, v.y, v.z)?;
        }
        for uv in uvs {
            writeln!(w, "vt {} {}", uv[0], uv[1])?;
        }
        for t in triangles {
            writeln!(
                w,
                "f {}/{} {}/{} {}/{}",
                t[0] + 1,
                t[0] + 1,
                t[1] + 1,
                t[1] + 1,
                t[2] + 1,
                t[2] + 1,
            )?;
        }
        w.flush()?;
        Ok(())
    }
}
/// A node in an OBJ scene hierarchy.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjSceneNode {
    /// Name of this node.
    pub name: String,
    /// Local transform.
    pub transform: MeshTransform,
    /// Mesh index (if this node has geometry).
    pub mesh_index: Option<usize>,
    /// Child node indices.
    pub children: Vec<usize>,
}
/// Per-vertex RGBA colour.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObjVertexColor {
    /// Red channel \[0, 1\].
    pub r: f64,
    /// Green channel \[0, 1\].
    pub g: f64,
    /// Blue channel \[0, 1\].
    pub b: f64,
    /// Alpha channel \[0, 1\].
    pub a: f64,
}
#[allow(dead_code)]
impl ObjVertexColor {
    /// Construct a colour from (r, g, b, a).
    pub fn rgba(r: f64, g: f64, b: f64, a: f64) -> Self {
        Self { r, g, b, a }
    }
    /// Construct an opaque colour from (r, g, b).
    pub fn rgb(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b, a: 1.0 }
    }
    /// Return as `[r, g, b, a]` array.
    pub fn to_array(self) -> [f64; 4] {
        [self.r, self.g, self.b, self.a]
    }
    /// Blend two colours with weight `t ∈ [0,1]`.
    pub fn lerp(self, other: Self, t: f64) -> Self {
        Self {
            r: self.r + (other.r - self.r) * t,
            g: self.g + (other.g - self.g) * t,
            b: self.b + (other.b - self.b) * t,
            a: self.a + (other.a - self.a) * t,
        }
    }
}
/// A simple OBJ surface (tensor-product B-spline stub).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ObjSurface {
    /// Surface name.
    pub name: String,
    /// Degree in u direction.
    pub degree_u: usize,
    /// Degree in v direction.
    pub degree_v: usize,
    /// Control point indices (row-major, n_u x n_v).
    pub control_points: Vec<usize>,
    /// Number of control points in u direction.
    pub n_u: usize,
    /// Knot vector in u direction.
    pub knots_u: Vec<f64>,
    /// Knot vector in v direction.
    pub knots_v: Vec<f64>,
}
/// Writer for Wavefront MTL material library files.
#[allow(dead_code)]
pub struct MtlWriter;
#[allow(dead_code)]
impl MtlWriter {
    /// Generate an MTL file string for the given materials.
    pub fn write(materials: &[ObjMaterial]) -> String {
        let mut s = String::from("# OxiPhysics MTL export\n");
        for mat in materials {
            s.push_str(&format!("\nnewmtl {}\n", mat.name));
            s.push_str(&format!("Ka {} {} {}\n", mat.ka[0], mat.ka[1], mat.ka[2]));
            s.push_str(&format!("Kd {} {} {}\n", mat.kd[0], mat.kd[1], mat.kd[2]));
            s.push_str(&format!("Ks {} {} {}\n", mat.ks[0], mat.ks[1], mat.ks[2]));
            s.push_str(&format!("Ns {}\n", mat.ns));
            s.push_str(&format!("d {}\n", mat.dissolve));
            if let Some(ref tex) = mat.map_kd {
                s.push_str(&format!("map_Kd {}\n", tex));
            }
        }
        s
    }
}
/// OBJ mesh extended with per-vertex colours (non-standard OBJ extension).
///
/// Some exporters write vertex colours as extra columns on `v` lines:
/// `v x y z r g b` or `v x y z r g b a`.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ObjVertexColorMesh {
    /// Underlying geometry.
    pub mesh: ObjMesh,
    /// Per-vertex colours (same length as `mesh.vertices` if present).
    pub colors: Vec<ObjVertexColor>,
}
#[allow(dead_code)]
impl ObjVertexColorMesh {
    /// Parse an OBJ string that may contain vertex colour data.
    pub fn from_str(data: &str) -> std::result::Result<Self, String> {
        let mut vcmesh = ObjVertexColorMesh::default();
        let mut smoothing_group: u32 = 0;
        let mut current_material: Option<String> = None;
        let mut current_group_name: Option<String> = None;
        let mut current_group_start: usize = 0;
        for raw in data.lines() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if line.starts_with("g ") || line.starts_with("o ") {
                if let Some(ref name) = current_group_name {
                    let count = vcmesh.mesh.faces.len() - current_group_start;
                    if count > 0 {
                        vcmesh.mesh.groups.push(ObjGroup {
                            name: name.clone(),
                            face_start: current_group_start,
                            face_count: count,
                        });
                    }
                }
                current_group_name = Some(line[2..].trim().to_string());
                current_group_start = vcmesh.mesh.faces.len();
            } else if line.starts_with("usemtl ") {
                current_material = Some(line[7..].trim().to_string());
            } else if line.starts_with("s ") {
                let val = line[2..].trim();
                smoothing_group = if val == "off" || val == "0" {
                    0
                } else {
                    val.parse::<u32>().unwrap_or(0)
                };
            } else if line.starts_with("v ") {
                let p: Vec<&str> = line.split_whitespace().collect();
                if p.len() < 4 {
                    return Err(format!("Vertex line too short: {}", line));
                }
                let x: f64 = p[1]
                    .parse()
                    .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                let y: f64 = p[2]
                    .parse()
                    .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                let z: f64 = p[3]
                    .parse()
                    .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                vcmesh.mesh.vertices.push([x, y, z]);
                let r: f64 = p.get(4).and_then(|v| v.parse().ok()).unwrap_or(1.0);
                let g: f64 = p.get(5).and_then(|v| v.parse().ok()).unwrap_or(1.0);
                let b: f64 = p.get(6).and_then(|v| v.parse().ok()).unwrap_or(1.0);
                let a: f64 = p.get(7).and_then(|v| v.parse().ok()).unwrap_or(1.0);
                vcmesh.colors.push(ObjVertexColor { r, g, b, a });
            } else if line.starts_with("vn ") {
                let p: Vec<&str> = line.split_whitespace().collect();
                if p.len() >= 4 {
                    let x: f64 = p[1]
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                    let y: f64 = p[2]
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                    let z: f64 = p[3]
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                    vcmesh.mesh.normals.push([x, y, z]);
                }
            } else if line.starts_with("vt ") {
                let p: Vec<&str> = line.split_whitespace().collect();
                if p.len() >= 3 {
                    let u: f64 = p[1]
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                    let v: f64 = p[2]
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| e.to_string())?;
                    vcmesh.mesh.uvs.push([u, v]);
                }
            } else if line.starts_with("f ") {
                let p: Vec<&str> = line.split_whitespace().collect();
                let mut vis = Vec::new();
                let mut vts = Vec::new();
                let mut vns = Vec::new();
                let mut has_vt = false;
                let mut has_vn = false;
                for tok in &p[1..] {
                    let parts: Vec<&str> = tok.split('/').collect();
                    let vi: usize = parts[0].parse::<usize>().map_err(|e| e.to_string())? - 1;
                    vis.push(vi);
                    if parts.len() >= 2 && !parts[1].is_empty() {
                        has_vt = true;
                        let vt: usize = parts[1].parse::<usize>().map_err(|e| e.to_string())? - 1;
                        vts.push(vt);
                    }
                    if parts.len() >= 3 && !parts[2].is_empty() {
                        has_vn = true;
                        let vn: usize = parts[2].parse::<usize>().map_err(|e| e.to_string())? - 1;
                        vns.push(vn);
                    }
                }
                vcmesh.mesh.faces.push(ObjFace {
                    vertex_indices: vis,
                    normal_indices: if has_vn { Some(vns) } else { None },
                    uv_indices: if has_vt { Some(vts) } else { None },
                    smoothing_group,
                    material: current_material.clone(),
                });
            }
        }
        if let Some(ref name) = current_group_name {
            let count = vcmesh.mesh.faces.len() - current_group_start;
            if count > 0 {
                vcmesh.mesh.groups.push(ObjGroup {
                    name: name.clone(),
                    face_start: current_group_start,
                    face_count: count,
                });
            }
        }
        Ok(vcmesh)
    }
    /// Serialize this coloured mesh to an OBJ string (with colour extension).
    pub fn to_obj_str(&self) -> String {
        let mut s = String::from("# OxiPhysics OBJ export (vertex colours)\n");
        for (i, v) in self.mesh.vertices.iter().enumerate() {
            if let Some(c) = self.colors.get(i) {
                s.push_str(&format!(
                    "v {} {} {} {} {} {}\n",
                    v[0], v[1], v[2], c.r, c.g, c.b
                ));
            } else {
                s.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
            }
        }
        for vn in &self.mesh.normals {
            s.push_str(&format!("vn {} {} {}\n", vn[0], vn[1], vn[2]));
        }
        for vt in &self.mesh.uvs {
            s.push_str(&format!("vt {} {}\n", vt[0], vt[1]));
        }
        for face in &self.mesh.faces {
            s.push('f');
            for i in 0..face.vertex_indices.len() {
                let vi = face.vertex_indices[i] + 1;
                let vt_idx = face.uv_indices.as_ref().map(|uvs| uvs[i] + 1);
                let vn_idx = face.normal_indices.as_ref().map(|ns| ns[i] + 1);
                let tok = match (vt_idx, vn_idx) {
                    (Some(vt), Some(vn)) => format!(" {}/{}/{}", vi, vt, vn),
                    (None, Some(vn)) => format!(" {}//{}", vi, vn),
                    (Some(vt), None) => format!(" {}/{}", vi, vt),
                    (None, None) => format!(" {}", vi),
                };
                s.push_str(&tok);
            }
            s.push('\n');
        }
        s
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// QEM (Quadric Error Metric) mesh decimation — Garland & Heckbert 1997
// ─────────────────────────────────────────────────────────────────────────────

/// Symmetric 4×4 error quadric stored as the 10-element upper triangle.
///
/// Layout: `[a2, ab, ac, ad, b2, bc, bd, c2, cd, d2]` where the plane is
/// `ax + by + cz + d = 0`.
#[derive(Clone, Copy)]
struct Quadric {
    m: [f64; 10],
}

impl Quadric {
    const ZERO: Self = Self { m: [0.0; 10] };

    /// Build the fundamental-error quadric `K_p = pp^T` for plane `(a,b,c,d)`.
    fn from_plane(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self {
            m: [
                a * a,
                a * b,
                a * c,
                a * d,
                b * b,
                b * c,
                b * d,
                c * c,
                c * d,
                d * d,
            ],
        }
    }

    /// Add two quadrics (element-wise).
    fn add(&self, other: &Self) -> Self {
        let mut m = [0.0f64; 10];
        for (i, val) in m.iter_mut().enumerate() {
            *val = self.m[i] + other.m[i];
        }
        Self { m }
    }

    /// Evaluate the quadric error at position `p`.
    fn eval(&self, p: [f64; 3]) -> f64 {
        let [x, y, z] = p;
        let m = &self.m;
        m[0] * x * x
            + 2.0 * m[1] * x * y
            + 2.0 * m[2] * x * z
            + 2.0 * m[3] * x
            + m[4] * y * y
            + 2.0 * m[5] * y * z
            + 2.0 * m[6] * y
            + m[7] * z * z
            + 2.0 * m[8] * z
            + m[9]
    }

    /// Solve for the optimal collapse position using Cramer's rule on the
    /// upper-left 3×3 block.  Falls back to `midpoint` for singular systems.
    fn optimal_pos(&self, midpoint: [f64; 3]) -> [f64; 3] {
        let m = &self.m;
        let a00 = m[0];
        let a01 = m[1];
        let a02 = m[2];
        let a11 = m[4];
        let a12 = m[5];
        let a22 = m[7];
        let bx = -m[3];
        let by = -m[6];
        let bz = -m[8];

        let det = a00 * (a11 * a22 - a12 * a12) - a01 * (a01 * a22 - a12 * a02)
            + a02 * (a01 * a12 - a11 * a02);

        if det.abs() < 1.0e-12 {
            return midpoint;
        }
        let inv = 1.0 / det;
        let x = inv
            * (bx * (a11 * a22 - a12 * a12) - a01 * (by * a22 - a12 * bz)
                + a02 * (by * a12 - a11 * bz));
        let y = inv
            * (a00 * (by * a22 - a12 * bz) - bx * (a01 * a22 - a12 * a02)
                + a02 * (a01 * bz - by * a02));
        let z = inv
            * (a00 * (a11 * bz - by * a12) - a01 * (a01 * bz - by * a02)
                + bx * (a01 * a12 - a11 * a02));
        [x, y, z]
    }
}

/// Entry in the collapse priority queue.
#[derive(Clone)]
struct EdgeEntry {
    cost: f64,
    v0: usize,
    v1: usize,
    pos: [f64; 3],
    /// Epoch stamp for lazy deletion: entry is stale when it no longer
    /// matches `vertex_epochs[v0] ^ (vertex_epochs[v1] << 32)`.
    epoch: u64,
}

impl PartialEq for EdgeEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost
    }
}
impl Eq for EdgeEntry {}
impl PartialOrd for EdgeEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for EdgeEntry {
    /// Reverse order so `BinaryHeap` acts as a min-heap.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Core QEM decimation implementation.
fn decimate_qem_impl(mesh: &ObjMesh, target_faces: usize) -> ObjMesh {
    use std::collections::{BinaryHeap, HashSet};

    let nv = mesh.vertices.len();
    if nv == 0 || mesh.faces.is_empty() {
        return mesh.clone();
    }

    // ── 1. Collect triangles ─────────────────────────────────────────────────
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for face in &mesh.faces {
        let vi = &face.vertex_indices;
        if vi.len() < 3 {
            continue;
        }
        for k in 1..(vi.len() - 1) {
            if vi[0] < nv && vi[k] < nv && vi[k + 1] < nv {
                triangles.push([vi[0], vi[k], vi[k + 1]]);
            }
        }
    }
    if triangles.is_empty() {
        return mesh.clone();
    }

    // ── 2. Compute per-vertex quadrics ───────────────────────────────────────
    let mut quadrics = vec![Quadric::ZERO; nv];
    for &[i0, i1, i2] in &triangles {
        let v0 = mesh.vertices[i0];
        let v1 = mesh.vertices[i1];
        let v2 = mesh.vertices[i2];
        let ex = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let fy = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let nx = ex[1] * fy[2] - ex[2] * fy[1];
        let ny = ex[2] * fy[0] - ex[0] * fy[2];
        let nz = ex[0] * fy[1] - ex[1] * fy[0];
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len < 1.0e-15 {
            continue;
        }
        let (a, b, c) = (nx / len, ny / len, nz / len);
        let d = -(a * v0[0] + b * v0[1] + c * v0[2]);
        let q = Quadric::from_plane(a, b, c, d);
        quadrics[i0] = quadrics[i0].add(&q);
        quadrics[i1] = quadrics[i1].add(&q);
        quadrics[i2] = quadrics[i2].add(&q);
    }

    // ── 3. Collect unique undirected edges ───────────────────────────────────
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    for &[i0, i1, i2] in &triangles {
        for (a, b) in [(i0, i1), (i1, i2), (i0, i2)] {
            edge_set.insert((a.min(b), a.max(b)));
        }
    }

    // ── 4. State arrays ──────────────────────────────────────────────────────
    let mut vertices: Vec<Option<[f64; 3]>> = mesh.vertices.iter().map(|&v| Some(v)).collect();
    let mut vertex_epochs: Vec<u64> = vec![0u64; nv];
    let mut remap: Vec<usize> = (0..nv).collect();
    let mut heap: BinaryHeap<EdgeEntry> = BinaryHeap::new();
    let mut faces: Vec<Option<[usize; 3]>> = triangles.iter().map(|&t| Some(t)).collect();
    let mut active_face_count = faces.len();

    // ── Helper: resolve remap chain ──────────────────────────────────────────
    fn resolve(remap: &[usize], mut v: usize) -> usize {
        while remap[v] != v {
            v = remap[v];
        }
        v
    }

    // ── Helper: push an edge into the heap ───────────────────────────────────
    let push_edge = |heap: &mut BinaryHeap<EdgeEntry>,
                     v0: usize,
                     v1: usize,
                     qv: &[Quadric],
                     verts: &[Option<[f64; 3]>],
                     ve: &[u64]| {
        let Some(p0) = verts[v0] else { return };
        let Some(p1) = verts[v1] else { return };
        let combined = qv[v0].add(&qv[v1]);
        let mid = [
            (p0[0] + p1[0]) * 0.5,
            (p0[1] + p1[1]) * 0.5,
            (p0[2] + p1[2]) * 0.5,
        ];
        let pos = combined.optimal_pos(mid);
        let cost = combined.eval(pos);
        heap.push(EdgeEntry {
            cost,
            v0,
            v1,
            pos,
            epoch: ve[v0] ^ (ve[v1] << 32),
        });
    };

    for &(a, b) in &edge_set {
        push_edge(&mut heap, a, b, &quadrics, &vertices, &vertex_epochs);
    }

    // ── 5. Greedy collapse loop ───────────────────────────────────────────────
    while active_face_count > target_faces {
        let entry = match heap.pop() {
            Some(e) => e,
            None => break,
        };
        let v0 = resolve(&remap, entry.v0);
        let v1 = resolve(&remap, entry.v1);
        if v0 == v1 {
            continue;
        }
        if entry.epoch != vertex_epochs[v0] ^ (vertex_epochs[v1] << 32) {
            continue;
        }
        if vertices[v0].is_none() || vertices[v1].is_none() {
            continue;
        }

        // Collapse v1 into v0.
        vertices[v0] = Some(entry.pos);
        vertices[v1] = None;
        quadrics[v0] = quadrics[v0].add(&quadrics[v1]);
        remap[v1] = v0;
        vertex_epochs[v0] += 1;
        vertex_epochs[v1] += 1;

        // Degenerate face removal.
        for face in faces.iter_mut() {
            let tri = match face {
                Some(t) => t,
                None => continue,
            };
            let mut changed = false;
            for idx in tri.iter_mut() {
                let r = resolve(&remap, *idx);
                if r != *idx {
                    *idx = r;
                    changed = true;
                }
            }
            if changed && (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) {
                *face = None;
                active_face_count -= 1;
            }
        }

        // Re-queue edges incident to the surviving vertex.
        for &(a, b) in &edge_set {
            let ra = resolve(&remap, a);
            let rb = resolve(&remap, b);
            if ra == rb {
                continue;
            }
            if ra == v0 || rb == v0 {
                push_edge(&mut heap, ra, rb, &quadrics, &vertices, &vertex_epochs);
            }
        }
    }

    // ── 6. Rebuild ObjMesh ───────────────────────────────────────────────────
    let mut new_vertices: Vec<[f64; 3]> = Vec::new();
    let mut old_to_new: Vec<Option<usize>> = vec![None; nv];
    for (i, v) in vertices.iter().enumerate() {
        if let Some(pos) = v {
            old_to_new[i] = Some(new_vertices.len());
            new_vertices.push(*pos);
        }
    }

    let mut new_faces: Vec<ObjFace> = Vec::new();
    for face in faces.iter().flatten() {
        let a = resolve(&remap, face[0]);
        let b = resolve(&remap, face[1]);
        let c = resolve(&remap, face[2]);
        if a == b || b == c || a == c {
            continue;
        }
        let Some(na) = old_to_new[a] else { continue };
        let Some(nb) = old_to_new[b] else { continue };
        let Some(nc) = old_to_new[c] else { continue };
        new_faces.push(ObjFace {
            vertex_indices: vec![na, nb, nc],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        if new_faces.len() >= target_faces {
            break;
        }
    }

    ObjMesh {
        vertices: new_vertices,
        normals: Vec::new(),
        uvs: Vec::new(),
        faces: new_faces,
        groups: Vec::new(),
    }
}
