//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// A triangulated parametric mesh (f64 precision).
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ParamMesh {
    /// Vertices.
    pub vertices: Vec<ParamVertex>,
    /// Triangle indices (triples).
    pub indices: Vec<usize>,
}
impl ParamMesh {
    /// Create an empty mesh.
    #[allow(dead_code)]
    pub fn empty() -> Self {
        ParamMesh {
            vertices: Vec::new(),
            indices: Vec::new(),
        }
    }
    /// Number of triangles.
    #[allow(dead_code)]
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
    /// Number of vertices.
    #[allow(dead_code)]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    /// Compute the axis-aligned bounding box.
    #[allow(dead_code)]
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            for k in 0..3 {
                if v.pos[k] < lo[k] {
                    lo[k] = v.pos[k];
                }
                if v.pos[k] > hi[k] {
                    hi[k] = v.pos[k];
                }
            }
        }
        (lo, hi)
    }
    /// Surface area (sum of triangle areas).
    #[allow(dead_code)]
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0;
        let n = self.indices.len() / 3;
        for t in 0..n {
            let i0 = self.indices[t * 3];
            let i1 = self.indices[t * 3 + 1];
            let i2 = self.indices[t * 3 + 2];
            let a = self.vertices[i0].pos;
            let b = self.vertices[i1].pos;
            let c = self.vertices[i2].pos;
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cross = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];
            let len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
            area += 0.5 * len;
        }
        area
    }
    /// Translate all vertex positions by `delta`.
    #[allow(dead_code)]
    pub fn translate(&mut self, delta: [f64; 3]) {
        for v in &mut self.vertices {
            v.pos[0] += delta[0];
            v.pos[1] += delta[1];
            v.pos[2] += delta[2];
        }
    }
    /// Scale all vertex positions by `s` relative to the origin.
    #[allow(dead_code)]
    pub fn scale(&mut self, s: f64) {
        for v in &mut self.vertices {
            v.pos[0] *= s;
            v.pos[1] *= s;
            v.pos[2] *= s;
        }
    }
}
/// One vertex on a parametric surface: position, normal, UV.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ParamVertex {
    /// Position in world space.
    pub pos: [f64; 3],
    /// Outward unit normal.
    pub normal: [f64; 3],
    /// UV texture coordinates.
    pub uv: [f64; 2],
}
/// A complete mesh with positions, normals, UVs, and triangle indices.
///
/// Suitable for tangent-space computation, GPU upload, and mesh operations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MeshData {
    /// Vertex positions (world space).
    pub positions: Vec<[f64; 3]>,
    /// Per-vertex normals (unit length).
    pub normals: Vec<[f64; 3]>,
    /// Per-vertex UV coordinates (u, v).
    pub uvs: Vec<[f64; 2]>,
    /// Triangle indices (triples).
    pub indices: Vec<[usize; 3]>,
}
#[allow(dead_code)]
impl MeshData {
    /// Create an empty `MeshData`.
    pub fn empty() -> Self {
        Self {
            positions: Vec::new(),
            normals: Vec::new(),
            uvs: Vec::new(),
            indices: Vec::new(),
        }
    }
    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }
    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len()
    }
    /// Returns `true` if positions, normals, and uvs all have the same length.
    pub fn is_consistent(&self) -> bool {
        let n = self.positions.len();
        self.normals.len() == n && self.uvs.len() == n
    }
    /// Apply a uniform scale to all vertex positions in place.
    pub fn scale(&mut self, s: f64) {
        for p in &mut self.positions {
            p[0] *= s;
            p[1] *= s;
            p[2] *= s;
        }
    }
    /// Translate all vertices by `offset`.
    pub fn translate(&mut self, offset: [f64; 3]) {
        for p in &mut self.positions {
            p[0] += offset[0];
            p[1] += offset[1];
            p[2] += offset[2];
        }
    }
    /// Merge another `MeshData` into this one, adjusting index offsets.
    pub fn merge_from(&mut self, other: &MeshData) {
        let base = self.positions.len();
        self.positions.extend_from_slice(&other.positions);
        self.normals.extend_from_slice(&other.normals);
        self.uvs.extend_from_slice(&other.uvs);
        for &[a, b, c] in &other.indices {
            self.indices.push([a + base, b + base, c + base]);
        }
    }
    /// Compute the axis-aligned bounding box of all vertices.
    ///
    /// Returns `(min, max)` or `None` if the mesh is empty.
    pub fn aabb(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.positions.is_empty() {
            return None;
        }
        let mut mn = self.positions[0];
        let mut mx = self.positions[0];
        for &p in &self.positions[1..] {
            mn[0] = mn[0].min(p[0]);
            mn[1] = mn[1].min(p[1]);
            mn[2] = mn[2].min(p[2]);
            mx[0] = mx[0].max(p[0]);
            mx[1] = mx[1].max(p[1]);
            mx[2] = mx[2].max(p[2]);
        }
        Some((mn, mx))
    }
    /// Recalculate per-vertex normals by averaging adjacent face normals.
    pub fn recompute_normals(&mut self) {
        use crate::wireframe::triangle_normal;
        let n = self.positions.len();
        let mut accum = vec![[0.0_f64; 3]; n];
        for &[a, b, c] in &self.indices {
            if a >= n || b >= n || c >= n {
                continue;
            }
            let fn_ = triangle_normal(self.positions[a], self.positions[b], self.positions[c]);
            for idx in [a, b, c] {
                accum[idx][0] += fn_[0];
                accum[idx][1] += fn_[1];
                accum[idx][2] += fn_[2];
            }
        }
        self.normals = accum
            .into_iter()
            .map(|mut v| {
                let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-30);
                v[0] /= mag;
                v[1] /= mag;
                v[2] /= mag;
                v
            })
            .collect();
    }
    /// Flatten triangle indices to a `Vec`usize` (for interop with flat index buffers).
    pub fn flat_indices(&self) -> Vec<usize> {
        self.indices
            .iter()
            .flat_map(|&[a, b, c]| [a, b, c])
            .collect()
    }
}
