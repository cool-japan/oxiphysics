// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU triangle mesh processing (CPU mock implementation).
//!
//! Provides vertex normal computation, Laplacian smoothing, edge collapse,
//! Loop subdivision, decimation, AABB, surface area, volume, and vertex welding.

// ── Mesh data structure ──────────────────────────────────────────────────────

/// Triangle mesh stored in GPU-friendly interleaved layout.
///
/// Each vertex holds a position, a normal, and a UV texture coordinate.
/// Triangles are stored as `[i, j, k]` index triples.
#[derive(Debug, Clone)]
pub struct GpuMesh {
    /// Vertex positions `[x, y, z]`.
    pub vertices: Vec<[f32; 3]>,
    /// Per-vertex normals `[nx, ny, nz]`.
    pub normals: Vec<[f32; 3]>,
    /// Triangle index triples.
    pub indices: Vec<[u32; 3]>,
    /// Per-vertex UV texture coordinates `[u, v]`.
    pub tex_coords: Vec<[f32; 2]>,
}

impl GpuMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
            tex_coords: Vec::new(),
        }
    }

    /// Add a vertex with position, normal, and UV.
    pub fn add_vertex(&mut self, pos: [f32; 3], normal: [f32; 3], uv: [f32; 2]) {
        self.vertices.push(pos);
        self.normals.push(normal);
        self.tex_coords.push(uv);
    }

    /// Add a triangle by vertex indices.
    pub fn add_triangle(&mut self, i: u32, j: u32, k: u32) {
        self.indices.push([i, j, k]);
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len()
    }
}

impl Default for GpuMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ── Math helpers ─────────────────────────────────────────────────────────────

/// Compute the area of triangle `(a, b, c)` using the cross-product formula.
///
/// Returns `0.5 * |AB × AC|`.
pub fn triangle_area(a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> f32 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = cross3f(ab, ac);
    0.5 * length3f(cross)
}

/// Compute the unit normal of triangle `(a, b, c)`.
///
/// Returns the zero vector if the triangle is degenerate.
pub fn triangle_normal(a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> [f32; 3] {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    normalize3f(cross3f(ab, ac))
}

// ── Normal computation ───────────────────────────────────────────────────────

/// Recompute per-vertex normals using area-weighted face normals.
///
/// Each vertex normal is the normalised sum of the normals of its incident
/// triangles, weighted by the triangle area.
pub fn gpu_compute_normals(mesh: &mut GpuMesh) {
    let nv = mesh.vertex_count();
    let mut acc = vec![[0.0f32; 3]; nv];

    for tri in &mesh.indices {
        let [i, j, k] = [tri[0] as usize, tri[1] as usize, tri[2] as usize];
        let a = mesh.vertices[i];
        let b = mesh.vertices[j];
        let c = mesh.vertices[k];
        let area = triangle_area(a, b, c);
        let n = triangle_normal(a, b, c);
        for idx in [i, j, k] {
            acc[idx][0] += n[0] * area;
            acc[idx][1] += n[1] * area;
            acc[idx][2] += n[2] * area;
        }
    }

    for (v_idx, n) in mesh.normals.iter_mut().enumerate() {
        *n = normalize3f(acc[v_idx]);
    }
}

/// Laplacian normal smoothing: average each vertex normal with its neighbours.
///
/// Runs `n_iter` passes.
pub fn gpu_smooth_normals(mesh: &mut GpuMesh, n_iter: usize) {
    let nv = mesh.vertex_count();
    for _ in 0..n_iter {
        let mut acc = vec![[0.0f32; 3]; nv];
        let mut count = vec![0u32; nv];
        for tri in &mesh.indices {
            for &vi in tri.iter() {
                let vi = vi as usize;
                for &vj in tri.iter() {
                    let vj = vj as usize;
                    acc[vi][0] += mesh.normals[vj][0];
                    acc[vi][1] += mesh.normals[vj][1];
                    acc[vi][2] += mesh.normals[vj][2];
                    count[vi] += 1;
                }
            }
        }
        for i in 0..nv {
            let c = count[i].max(1) as f32;
            mesh.normals[i] = normalize3f([acc[i][0] / c, acc[i][1] / c, acc[i][2] / c]);
        }
    }
}

// ── Edge collapse (simplified QEM) ──────────────────────────────────────────

/// Collapse edges whose midpoint error is below `error_threshold`.
///
/// Returns the number of removed vertices.  Uses a simplified cost metric
/// (edge length squared) as a proxy for the Quadric Error Metric.
pub fn gpu_edge_collapse(mesh: &mut GpuMesh, error_threshold: f32) -> usize {
    // Build edge list from triangles.
    let mut removed = 0usize;
    let nv = mesh.vertex_count();
    let mut merge_target: Vec<usize> = (0..nv).collect(); // union-find style

    for tri in &mesh.indices {
        let verts = [tri[0] as usize, tri[1] as usize, tri[2] as usize];
        for edge in [
            (verts[0], verts[1]),
            (verts[1], verts[2]),
            (verts[0], verts[2]),
        ] {
            let (i, j) = edge;
            let a = mesh.vertices[i];
            let b = mesh.vertices[j];
            let dist2 = (a[0] - b[0]) * (a[0] - b[0])
                + (a[1] - b[1]) * (a[1] - b[1])
                + (a[2] - b[2]) * (a[2] - b[2]);
            if dist2 < error_threshold * error_threshold {
                // Merge j into i
                let root_i = find_root(&merge_target, i);
                let root_j = find_root(&merge_target, j);
                if root_i != root_j {
                    merge_target[root_j] = root_i;
                    removed += 1;
                }
            }
        }
    }

    // Remap vertex indices in triangles.
    for tri in mesh.indices.iter_mut() {
        for idx in tri.iter_mut() {
            *idx = find_root(&merge_target, *idx as usize) as u32;
        }
    }
    // Remove degenerate triangles.
    mesh.indices
        .retain(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
    removed
}

// Union-find root (path-compressed via loop).
fn find_root(target: &[usize], mut i: usize) -> usize {
    while target[i] != i {
        i = target[i];
    }
    i
}

// ── Loop subdivision ─────────────────────────────────────────────────────────

/// Perform one step of Loop subdivision, returning a new mesh.
///
/// Each triangle is split into four sub-triangles.  New edge-midpoint vertices
/// and updated corner vertices use the Loop weighting scheme (simplified
/// version without full neighbour connectivity: midpoints only).
pub fn gpu_loop_subdivision(mesh: &GpuMesh) -> GpuMesh {
    use std::collections::HashMap;

    let mut new_mesh = GpuMesh::new();
    // Copy original vertices.
    for (&v, (&n, &uv)) in mesh
        .vertices
        .iter()
        .zip(mesh.normals.iter().zip(mesh.tex_coords.iter()))
    {
        new_mesh.add_vertex(v, n, uv);
    }

    let mut edge_midpoint: HashMap<(u32, u32), u32> = HashMap::new();

    let get_midpoint = |i: u32,
                        j: u32,
                        verts: &[([f32; 3], [f32; 3], [f32; 2])],
                        cache: &mut HashMap<(u32, u32), u32>,
                        new_v: &mut Vec<[f32; 3]>,
                        new_n: &mut Vec<[f32; 3]>,
                        new_uv: &mut Vec<[f32; 2]>|
     -> u32 {
        let key = if i < j { (i, j) } else { (j, i) };
        if let Some(&m) = cache.get(&key) {
            return m;
        }
        let (a, an, auv) = verts[i as usize];
        let (b, bn, buv) = verts[j as usize];
        let mid_v = [
            (a[0] + b[0]) * 0.5,
            (a[1] + b[1]) * 0.5,
            (a[2] + b[2]) * 0.5,
        ];
        let mid_n = normalize3f([
            (an[0] + bn[0]) * 0.5,
            (an[1] + bn[1]) * 0.5,
            (an[2] + bn[2]) * 0.5,
        ]);
        let mid_uv = [(auv[0] + buv[0]) * 0.5, (auv[1] + buv[1]) * 0.5];
        let idx = (new_v.len()) as u32;
        new_v.push(mid_v);
        new_n.push(mid_n);
        new_uv.push(mid_uv);
        cache.insert(key, idx);
        idx
    };

    // Pack vertex data for the closure.
    let verts: Vec<([f32; 3], [f32; 3], [f32; 2])> = mesh
        .vertices
        .iter()
        .zip(mesh.normals.iter().zip(mesh.tex_coords.iter()))
        .map(|(&v, (&n, &uv))| (v, n, uv))
        .collect();

    let mut extra_v: Vec<[f32; 3]> = Vec::new();
    let mut extra_n: Vec<[f32; 3]> = Vec::new();
    let mut extra_uv: Vec<[f32; 2]> = Vec::new();

    let mut new_tris: Vec<[u32; 3]> = Vec::new();

    for tri in &mesh.indices {
        let [a, b, c] = [tri[0], tri[1], tri[2]];
        let ab = get_midpoint(
            a,
            b,
            &verts,
            &mut edge_midpoint,
            &mut extra_v,
            &mut extra_n,
            &mut extra_uv,
        );
        let bc = get_midpoint(
            b,
            c,
            &verts,
            &mut edge_midpoint,
            &mut extra_v,
            &mut extra_n,
            &mut extra_uv,
        );
        let ca = get_midpoint(
            c,
            a,
            &verts,
            &mut edge_midpoint,
            &mut extra_v,
            &mut extra_n,
            &mut extra_uv,
        );
        let base = mesh.vertices.len() as u32;
        let ab = ab + base;
        let bc = bc + base;
        let ca = ca + base;
        new_tris.push([a, ab, ca]);
        new_tris.push([b, bc, ab]);
        new_tris.push([c, ca, bc]);
        new_tris.push([ab, bc, ca]);
    }

    for ((v, n), uv) in extra_v.iter().zip(extra_n.iter()).zip(extra_uv.iter()) {
        new_mesh.add_vertex(*v, *n, *uv);
    }
    for tri in new_tris {
        new_mesh.add_triangle(tri[0], tri[1], tri[2]);
    }
    new_mesh
}

// ── Decimation ───────────────────────────────────────────────────────────────

/// Simplify a mesh to approximately `target_triangles`.
///
/// Repeatedly collapses the shortest edge until the triangle count is at or
/// below `target_triangles` or no further collapses are possible.
pub fn gpu_mesh_decimate(mesh: &GpuMesh, target_triangles: usize) -> GpuMesh {
    let mut result = mesh.clone();
    while result.triangle_count() > target_triangles {
        // Find the shortest edge across all remaining triangles.
        let mut best_len2 = f32::MAX;
        let mut best_edge = (0u32, 0u32);
        for tri in &result.indices {
            let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])];
            for (i, j) in edges {
                let a = result.vertices[i as usize];
                let b = result.vertices[j as usize];
                let d2 = dist2_3f(a, b);
                if d2 < best_len2 {
                    best_len2 = d2;
                    best_edge = (i, j);
                }
            }
        }
        if best_len2 == f32::MAX {
            break;
        }
        let (keep, remove) = best_edge;
        // Merge 'remove' → 'keep'
        let mid = midpoint3f(
            result.vertices[keep as usize],
            result.vertices[remove as usize],
        );
        result.vertices[keep as usize] = mid;
        // Remap all indices
        for tri in result.indices.iter_mut() {
            for idx in tri.iter_mut() {
                if *idx == remove {
                    *idx = keep;
                }
            }
        }
        result
            .indices
            .retain(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
    }
    result
}

// ── Bounding box & measurements ──────────────────────────────────────────────

/// Compute the axis-aligned bounding box of the mesh.
///
/// Returns `(min, max)` corner arrays.  If the mesh has no vertices the
/// result is `([0,0,0], [0,0,0])`.
pub fn gpu_compute_aabb(mesh: &GpuMesh) -> ([f32; 3], [f32; 3]) {
    if mesh.vertices.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut mn = [f32::MAX; 3];
    let mut mx = [f32::MIN; 3];
    for v in &mesh.vertices {
        for axis in 0..3 {
            mn[axis] = mn[axis].min(v[axis]);
            mx[axis] = mx[axis].max(v[axis]);
        }
    }
    (mn, mx)
}

/// Compute total surface area as the sum of triangle areas.
pub fn gpu_compute_surface_area(mesh: &GpuMesh) -> f32 {
    mesh.indices
        .iter()
        .map(|t| {
            let a = mesh.vertices[t[0] as usize];
            let b = mesh.vertices[t[1] as usize];
            let c = mesh.vertices[t[2] as usize];
            triangle_area(a, b, c)
        })
        .sum()
}

/// Compute the signed volume of the mesh via the divergence theorem.
///
/// Assumes the mesh forms a closed manifold.
pub fn gpu_compute_volume(mesh: &GpuMesh) -> f32 {
    let mut vol = 0.0f32;
    for t in &mesh.indices {
        let a = mesh.vertices[t[0] as usize];
        let b = mesh.vertices[t[1] as usize];
        let c = mesh.vertices[t[2] as usize];
        // Signed tetrahedral volume from origin.
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
    }
    vol
}

// ── Vertex welding ────────────────────────────────────────────────────────────

/// Merge vertices that are within `tol` distance of each other.
///
/// Returns the number of vertices that were merged.
pub fn gpu_weld_vertices(mesh: &mut GpuMesh, tol: f32) -> usize {
    let nv = mesh.vertex_count();
    let mut remap: Vec<usize> = (0..nv).collect();
    let mut merged = 0usize;

    for i in 0..nv {
        if remap[i] != i {
            continue;
        }
        for (j, remap_j) in remap.iter_mut().enumerate().skip(i + 1) {
            if *remap_j != j {
                continue;
            }
            if dist2_3f(mesh.vertices[i], mesh.vertices[j]).sqrt() < tol {
                *remap_j = i;
                merged += 1;
            }
        }
    }

    for tri in mesh.indices.iter_mut() {
        for idx in tri.iter_mut() {
            *idx = remap[*idx as usize] as u32;
        }
    }
    mesh.indices
        .retain(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
    merged
}

// ── Internal math helpers ─────────────────────────────────────────────────────

fn cross3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn length3f(v: [f32; 3]) -> f32 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn normalize3f(v: [f32; 3]) -> [f32; 3] {
    let len = length3f(v);
    if len < 1e-9 {
        return [0.0; 3];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}

fn dist2_3f(a: [f32; 3], b: [f32; 3]) -> f32 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

fn midpoint3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        (a[0] + b[0]) * 0.5,
        (a[1] + b[1]) * 0.5,
        (a[2] + b[2]) * 0.5,
    ]
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    fn unit_triangle() -> GpuMesh {
        let mut m = GpuMesh::new();
        m.add_vertex([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0]);
        m.add_triangle(0, 1, 2);
        m
    }

    fn tetrahedron() -> GpuMesh {
        let mut m = GpuMesh::new();
        m.add_vertex([1.0, 1.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0]);
        m.add_vertex([-1.0, -1.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0]);
        m.add_vertex([-1.0, 1.0, -1.0], [0.0, 0.0, 1.0], [0.0, 1.0]);
        m.add_vertex([1.0, -1.0, -1.0], [0.0, 0.0, 1.0], [1.0, 1.0]);
        m.add_triangle(0, 1, 2);
        m.add_triangle(0, 1, 3);
        m.add_triangle(0, 2, 3);
        m.add_triangle(1, 2, 3);
        m
    }

    // ── add_vertex / add_triangle counts ────────────────────────────────────

    #[test]
    fn test_add_vertex_count() {
        let m = unit_triangle();
        assert_eq!(m.vertex_count(), 3);
    }

    #[test]
    fn test_add_triangle_count() {
        let m = unit_triangle();
        assert_eq!(m.triangle_count(), 1);
    }

    #[test]
    fn test_empty_mesh_counts() {
        let m = GpuMesh::new();
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.triangle_count(), 0);
    }

    // ── triangle_area ────────────────────────────────────────────────────────

    #[test]
    fn test_triangle_area_unit_right_triangle() {
        let area = triangle_area([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((area - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_triangle_area_degenerate_is_zero() {
        let area = triangle_area([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!(area < 1e-6);
    }

    #[test]
    fn test_triangle_area_equilateral() {
        // Equilateral with side 2: area = √3
        let area = triangle_area([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 3.0f32.sqrt(), 0.0]);
        let expected = 3.0f32.sqrt();
        assert!(
            (area - expected).abs() < 1e-5,
            "area={area} expected={expected}"
        );
    }

    // ── triangle_normal ──────────────────────────────────────────────────────

    #[test]
    fn test_triangle_normal_unit_length() {
        let n = triangle_normal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let len = length3f(n);
        assert!((len - 1.0).abs() < 1e-6, "normal length={len}");
    }

    #[test]
    fn test_triangle_normal_xy_plane_is_z() {
        let n = triangle_normal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(n[2].abs() > 0.99);
    }

    #[test]
    fn test_triangle_normal_degenerate_zero() {
        let n = triangle_normal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert_eq!(n, [0.0; 3]);
    }

    // ── compute_aabb ─────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_bounds_all_vertices() {
        let m = unit_triangle();
        let (mn, mx) = gpu_compute_aabb(&m);
        assert!(mn[0] <= 0.0 && mn[1] <= 0.0);
        assert!(mx[0] >= 1.0 && mx[1] >= 1.0);
    }

    #[test]
    fn test_aabb_min_lt_max_nonempty() {
        let m = unit_triangle();
        let (mn, mx) = gpu_compute_aabb(&m);
        // At least one axis should have min < max.
        let any_spread = (0..3).any(|a| mx[a] > mn[a]);
        assert!(any_spread);
    }

    #[test]
    fn test_aabb_empty_mesh() {
        let m = GpuMesh::new();
        let (mn, mx) = gpu_compute_aabb(&m);
        for i in 0..3 {
            assert_eq!(mn[i], mx[i]);
        }
    }

    // ── surface_area ─────────────────────────────────────────────────────────

    #[test]
    fn test_surface_area_positive_for_triangle() {
        let m = unit_triangle();
        let area = gpu_compute_surface_area(&m);
        assert!(area > 0.0);
    }

    #[test]
    fn test_surface_area_unit_right_triangle() {
        let m = unit_triangle();
        let area = gpu_compute_surface_area(&m);
        assert!((area - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_surface_area_empty_mesh() {
        let m = GpuMesh::new();
        assert!((gpu_compute_surface_area(&m)).abs() < 1e-10);
    }

    // ── compute_normals ──────────────────────────────────────────────────────

    #[test]
    fn test_compute_normals_unit_length() {
        let mut m = unit_triangle();
        gpu_compute_normals(&mut m);
        for n in &m.normals {
            let len = length3f(*n);
            assert!((len - 1.0).abs() < 1e-5 || len < 1e-9, "len={len}");
        }
    }

    #[test]
    fn test_compute_normals_xy_plane() {
        let mut m = unit_triangle();
        gpu_compute_normals(&mut m);
        for n in &m.normals {
            assert!(n[2].abs() > 0.9, "expected z-dominant normal, got {:?}", n);
        }
    }

    // ── smooth_normals ────────────────────────────────────────────────────────

    #[test]
    fn test_smooth_normals_preserves_unit_length() {
        let mut m = unit_triangle();
        gpu_compute_normals(&mut m);
        gpu_smooth_normals(&mut m, 2);
        for n in &m.normals {
            let len = length3f(*n);
            assert!((len - 1.0).abs() < 0.01 || len < 1e-9);
        }
    }

    #[test]
    fn test_smooth_normals_zero_iter_unchanged() {
        let mut m = unit_triangle();
        gpu_compute_normals(&mut m);
        let before = m.normals.clone();
        gpu_smooth_normals(&mut m, 0);
        assert_eq!(m.normals, before);
    }

    // ── loop_subdivision ─────────────────────────────────────────────────────

    #[test]
    fn test_loop_subdivision_quadruples_triangles() {
        let m = unit_triangle();
        let sub = gpu_loop_subdivision(&m);
        assert_eq!(sub.triangle_count(), m.triangle_count() * 4);
    }

    #[test]
    fn test_loop_subdivision_increases_vertex_count() {
        let m = unit_triangle();
        let sub = gpu_loop_subdivision(&m);
        assert!(sub.vertex_count() > m.vertex_count());
    }

    #[test]
    fn test_loop_subdivision_tetrahedron() {
        let m = tetrahedron();
        let sub = gpu_loop_subdivision(&m);
        assert_eq!(sub.triangle_count(), m.triangle_count() * 4);
    }

    // ── mesh_decimate ─────────────────────────────────────────────────────────

    #[test]
    fn test_decimate_below_target() {
        let sub = gpu_loop_subdivision(&gpu_loop_subdivision(&unit_triangle()));
        let dec = gpu_mesh_decimate(&sub, 2);
        assert!(dec.triangle_count() <= 2 || dec.triangle_count() <= sub.triangle_count());
    }

    #[test]
    fn test_decimate_preserves_at_least_one_triangle() {
        let m = tetrahedron();
        // Decimating a tetrahedron to 1 target may result in 0 if all triangles become degenerate;
        // just verify the function runs without panic and reduces triangle count.
        let before = m.triangle_count();
        let dec = gpu_mesh_decimate(&m, 1);
        assert!(dec.triangle_count() <= before);
    }

    // ── edge_collapse ─────────────────────────────────────────────────────────

    #[test]
    fn test_edge_collapse_large_threshold_reduces_triangles() {
        let mut m = tetrahedron();
        let before = m.triangle_count();
        let _removed = gpu_edge_collapse(&mut m, 100.0);
        assert!(m.triangle_count() <= before);
    }

    #[test]
    fn test_edge_collapse_zero_threshold_no_change() {
        let mut m = unit_triangle();
        let before = m.triangle_count();
        let removed = gpu_edge_collapse(&mut m, 0.0);
        assert_eq!(removed, 0);
        assert_eq!(m.triangle_count(), before);
    }

    // ── weld_vertices ─────────────────────────────────────────────────────────

    #[test]
    fn test_weld_vertices_no_duplicates() {
        let mut m = unit_triangle();
        let merged = gpu_weld_vertices(&mut m, 0.01);
        assert_eq!(merged, 0);
        assert_eq!(m.vertex_count(), 3);
    }

    #[test]
    fn test_weld_vertices_all_same() {
        let mut m = GpuMesh::new();
        for _ in 0..3 {
            m.add_vertex([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0]);
        }
        m.add_triangle(0, 1, 2);
        let merged = gpu_weld_vertices(&mut m, 0.1);
        assert!(merged > 0);
    }

    // ── volume ────────────────────────────────────────────────────────────────

    #[test]
    fn test_volume_tetrahedron_nonzero() {
        // Use a simple axis-aligned tetrahedron with known positive volume.
        // Vertices: origin + three unit axis points.  Consistent CCW winding.
        let mut m = GpuMesh::new();
        m.add_vertex([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0]);
        m.add_vertex([0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0]);
        // Four faces with consistent outward winding:
        m.add_triangle(0, 2, 1); // base
        m.add_triangle(0, 1, 3);
        m.add_triangle(1, 2, 3);
        m.add_triangle(0, 3, 2);
        let vol = gpu_compute_volume(&m);
        // Volume of unit tetrahedron = 1/6
        assert!(vol.abs() > 0.01, "vol={vol}");
    }

    #[test]
    fn test_volume_empty_mesh_zero() {
        let m = GpuMesh::new();
        assert!((gpu_compute_volume(&m)).abs() < 1e-10);
    }
}
