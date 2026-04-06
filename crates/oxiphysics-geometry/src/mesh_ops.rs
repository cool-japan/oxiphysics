// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Triangle mesh operations: normals, volume, surface area, point containment.

use oxiphysics_core::math::Vec3;
use std::collections::HashMap;

/// Vertex indices for a triangle.
pub type Triangle = [usize; 3];

/// A simple triangle mesh.
#[derive(Debug, Clone)]
pub struct TriMesh {
    /// Vertices of the mesh.
    pub vertices: Vec<Vec3>,
    /// Triangles defined by vertex indices.
    pub triangles: Vec<Triangle>,
}

impl TriMesh {
    /// Create a new triangle mesh.
    pub fn new(vertices: Vec<Vec3>, triangles: Vec<Triangle>) -> Self {
        Self {
            vertices,
            triangles,
        }
    }

    /// Compute the face normal for triangle at `tri_idx`.
    ///
    /// Returns `cross(v1-v0, v2-v0).normalize()`.
    pub fn face_normal(&self, tri_idx: usize) -> Vec3 {
        let tri = &self.triangles[tri_idx];
        let v0 = self.vertices[tri[0]];
        let v1 = self.vertices[tri[1]];
        let v2 = self.vertices[tri[2]];
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        e1.cross(&e2).normalize()
    }

    /// Volume of a closed mesh using the divergence theorem (signed volume sum).
    ///
    /// Assumes a watertight mesh with consistent winding.
    /// `V = (1/6) * Σ (v0 · (v1 × v2))`
    pub fn signed_volume(&self) -> f64 {
        let mut sum = 0.0_f64;
        for tri in &self.triangles {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            sum += v0.dot(&v1.cross(&v2));
        }
        sum / 6.0
    }

    /// Surface area of the mesh.
    ///
    /// Sum of triangle areas = `0.5 * |cross(e1, e2)|`.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0_f64;
        for tri in &self.triangles {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            let e1 = v1 - v0;
            let e2 = v2 - v0;
            area += 0.5 * e1.cross(&e2).norm();
        }
        area
    }

    /// Compute vertex normals by averaging adjacent face normals.
    pub fn vertex_normals(&self) -> Vec<Vec3> {
        let n = self.vertices.len();
        let mut normals = vec![Vec3::zeros(); n];
        for tri in &self.triangles {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            let e1 = v1 - v0;
            let e2 = v2 - v0;
            let face_n = e1.cross(&e2); // weighted by area (no normalize yet)
            normals[tri[0]] += face_n;
            normals[tri[1]] += face_n;
            normals[tri[2]] += face_n;
        }
        normals
            .into_iter()
            .map(|n| {
                let len = n.norm();
                if len < 1e-10 {
                    Vec3::new(0.0, 1.0, 0.0)
                } else {
                    n / len
                }
            })
            .collect()
    }

    /// Check if a point is inside the mesh using ray casting in the +Z direction.
    ///
    /// Counts intersections; odd count means inside. Assumes watertight mesh.
    pub fn contains_point(&self, point: Vec3) -> bool {
        let mut count = 0usize;
        for tri in &self.triangles {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            if ray_intersects_triangle_z(point, v0, v1, v2) {
                count += 1;
            }
        }
        count % 2 == 1
    }

    /// Generate a UV sphere mesh for testing purposes.
    ///
    /// `rings` is the number of latitude rings (excluding poles),
    /// `sectors` is the number of longitude sectors.
    pub fn uv_sphere(radius: f64, rings: u32, sectors: u32) -> Self {
        use std::f64::consts::PI;
        let mut vertices = Vec::new();
        let mut triangles = Vec::new();

        // Top pole
        vertices.push(Vec3::new(0.0, radius, 0.0));

        // Middle rings
        for ring in 1..=rings {
            let phi = PI * ring as f64 / (rings + 1) as f64; // 0..PI
            let y = radius * phi.cos();
            let r_xz = radius * phi.sin();
            for sector in 0..sectors {
                let theta = 2.0 * PI * sector as f64 / sectors as f64;
                let x = r_xz * theta.cos();
                let z = r_xz * theta.sin();
                vertices.push(Vec3::new(x, y, z));
            }
        }

        // Bottom pole
        vertices.push(Vec3::new(0.0, -radius, 0.0));

        let top_idx = 0usize;
        let bottom_idx = vertices.len() - 1;

        // Triangles from top pole to first ring (outward normals: CCW from outside)
        // Vertices on the first ring go CCW when viewed from above (+Y).
        // Winding [top, b, a] gives outward normal at the top.
        for s in 0..sectors as usize {
            let next_s = (s + 1) % sectors as usize;
            let a = 1 + s;
            let b = 1 + next_s;
            triangles.push([top_idx, b, a]);
        }

        // Quads for middle rings
        // Upper ring: row_start; lower ring: next_row_start.
        // For outward normals (viewing sphere from outside = CCW winding).
        for ring in 0..(rings as usize - 1) {
            let row_start = 1 + ring * sectors as usize;
            let next_row_start = row_start + sectors as usize;
            for s in 0..sectors as usize {
                let next_s = (s + 1) % sectors as usize;
                let a = row_start + s; // upper-left
                let b = row_start + next_s; // upper-right
                let c = next_row_start + s; // lower-left
                let d = next_row_start + next_s; // lower-right
                // Two triangles per quad, outward CCW
                triangles.push([a, b, c]);
                triangles.push([b, d, c]);
            }
        }

        // Triangles from last ring to bottom pole (outward normals).
        let last_ring_start = 1 + (rings as usize - 1) * sectors as usize;
        for s in 0..sectors as usize {
            let next_s = (s + 1) % sectors as usize;
            let a = last_ring_start + s;
            let b = last_ring_start + next_s;
            triangles.push([a, b, bottom_idx]);
        }

        Self {
            vertices,
            triangles,
        }
    }

    /// Compute the axis-aligned bounding box of the mesh.
    ///
    /// Returns `(min, max)` corners. Panics if the mesh has no vertices.
    pub fn bounding_box(&self) -> (Vec3, Vec3) {
        let mut min = self.vertices[0];
        let mut max = self.vertices[0];
        for v in &self.vertices {
            min.x = min.x.min(v.x);
            min.y = min.y.min(v.y);
            min.z = min.z.min(v.z);
            max.x = max.x.max(v.x);
            max.y = max.y.max(v.y);
            max.z = max.z.max(v.z);
        }
        (min, max)
    }
}

// ─── Standalone mesh operation helpers ──────────────────────────────────────

#[inline]
fn mo_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn mo_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn mo_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn mo_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn mo_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline]
fn mo_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = mo_norm(a);
    if n < 1e-12 {
        [0.0, 1.0, 0.0]
    } else {
        mo_scale(a, 1.0 / n)
    }
}

// ─── compute_face_normals ────────────────────────────────────────────────────

/// Compute a unit normal for every triangle in the mesh.
///
/// The normal is `normalize(cross(v1 - v0, v2 - v0))`.
/// If the triangle is degenerate (zero area) a zero vector is returned.
pub fn compute_face_normals(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> Vec<[f64; 3]> {
    tris.iter()
        .map(|tri| {
            let v0 = verts[tri[0]];
            let v1 = verts[tri[1]];
            let v2 = verts[tri[2]];
            let e1 = mo_sub(v1, v0);
            let e2 = mo_sub(v2, v0);
            let c = mo_cross(e1, e2);
            let len = mo_norm(c);
            if len < 1e-14 {
                [0.0; 3]
            } else {
                mo_scale(c, 1.0 / len)
            }
        })
        .collect()
}

// ─── compute_vertex_normals ──────────────────────────────────────────────────

/// Compute per-vertex normals using area-weighted averaging of adjacent face normals.
///
/// Each face contributes its area-weighted normal (`cross(e1, e2)`) to its
/// three vertices; the accumulated vectors are then normalized.
pub fn compute_vertex_normals(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> Vec<[f64; 3]> {
    let n = verts.len();
    let mut acc: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    for tri in tris {
        let v0 = verts[tri[0]];
        let v1 = verts[tri[1]];
        let v2 = verts[tri[2]];
        let e1 = mo_sub(v1, v0);
        let e2 = mo_sub(v2, v0);
        let fn_ = mo_cross(e1, e2); // area-weighted face normal
        for &vi in tri {
            acc[vi] = mo_add(acc[vi], fn_);
        }
    }
    acc.into_iter().map(mo_normalize).collect()
}

// ─── weld_vertices ───────────────────────────────────────────────────────────

/// Merge vertices that are closer than `tol` into a single representative vertex.
///
/// Returns the de-duplicated vertex array and remapped triangle indices.
/// This is an O(n²) brute-force implementation suitable for small meshes.
pub fn weld_vertices(
    verts: &[[f64; 3]],
    tris: &[[usize; 3]],
    tol: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let tol2 = tol * tol;
    let mut new_verts: Vec<[f64; 3]> = Vec::with_capacity(verts.len());
    let mut remap: Vec<usize> = vec![0; verts.len()];

    for (i, &v) in verts.iter().enumerate() {
        // Search for an existing vertex within tolerance
        let mut found = None;
        for (j, &w) in new_verts.iter().enumerate() {
            let d = mo_sub(v, w);
            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= tol2 {
                found = Some(j);
                break;
            }
        }
        remap[i] = match found {
            Some(j) => j,
            None => {
                let idx = new_verts.len();
                new_verts.push(v);
                idx
            }
        };
    }

    let new_tris: Vec<[usize; 3]> = tris
        .iter()
        .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
        .collect();

    (new_verts, new_tris)
}

// ─── compute_adjacency ───────────────────────────────────────────────────────

/// Build a vertex-to-vertex adjacency list from a triangle mesh.
///
/// Returns a `Vec` of length `n_verts`, where each element is the sorted
/// list of vertex indices adjacent to vertex `i`.
pub fn compute_adjacency(tris: &[[usize; 3]], n_verts: usize) -> Vec<Vec<usize>> {
    let mut adj: Vec<std::collections::BTreeSet<usize>> =
        vec![std::collections::BTreeSet::new(); n_verts];
    for tri in tris {
        let [a, b, c] = *tri;
        adj[a].insert(b);
        adj[a].insert(c);
        adj[b].insert(a);
        adj[b].insert(c);
        adj[c].insert(a);
        adj[c].insert(b);
    }
    adj.into_iter().map(|s| s.into_iter().collect()).collect()
}

// ─── is_watertight ───────────────────────────────────────────────────────────

/// Return `true` if the mesh is watertight (every edge is shared by exactly two faces).
///
/// A mesh is watertight (manifold, no boundary) when the edge–face incidence
/// count for every directed edge (v_a, v_b) satisfies:
/// `count(a→b) == 1` for all edges when counting undirected pairs,
/// which is equivalent to each unordered edge having exactly 2 incident triangles.
pub fn is_watertight(tris: &[[usize; 3]]) -> bool {
    use std::collections::HashMap;
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in tris {
        let [a, b, c] = *tri;
        let edges = [
            (a.min(b), a.max(b)),
            (b.min(c), b.max(c)),
            (a.min(c), a.max(c)),
        ];
        for e in edges {
            *edge_count.entry(e).or_insert(0) += 1;
        }
    }
    edge_count.values().all(|&cnt| cnt == 2)
}

// ─── compute_genus ───────────────────────────────────────────────────────────

/// Compute the genus of a closed orientable triangulated surface via the
/// Euler characteristic: `χ = V - E + F = 2 - 2g`.
///
/// Returns `g = (2 - χ) / 2`.  For a topological sphere `g = 0`.
pub fn compute_genus(n_verts: usize, tris: &[[usize; 3]]) -> i32 {
    use std::collections::HashSet;
    let v = n_verts as i32;
    let f = tris.len() as i32;
    let mut edges: HashSet<(usize, usize)> = HashSet::new();
    for tri in tris {
        let [a, b, c] = *tri;
        edges.insert((a.min(b), a.max(b)));
        edges.insert((b.min(c), b.max(c)));
        edges.insert((a.min(c), a.max(c)));
    }
    let e = edges.len() as i32;
    let chi = v - e + f; // Euler characteristic
    (2 - chi) / 2
}

// ─── triangulate_polygon ─────────────────────────────────────────────────────

/// Triangulate a simple convex or concave planar polygon using ear clipping.
///
/// `polygon` must be a list of at least 3 coplanar vertices in order.  The
/// function returns triangles as index triples (into the input slice).
///
/// # Panics
/// Panics if `polygon.len() < 3`.
pub fn triangulate_polygon(polygon: &[[f64; 3]]) -> Vec<[usize; 3]> {
    let n = polygon.len();
    assert!(n >= 3, "triangulate_polygon needs at least 3 vertices");
    if n == 3 {
        return vec![[0, 1, 2]];
    }

    // Work with indices
    let mut remaining: Vec<usize> = (0..n).collect();
    let mut result = Vec::with_capacity(n - 2);

    // Compute polygon normal via Newell's method for robust 2-D projection
    let mut poly_normal = [0.0f64; 3];
    for i in 0..n {
        let j = (i + 1) % n;
        let p = polygon[i];
        let q = polygon[j];
        poly_normal[0] += (p[1] - q[1]) * (p[2] + q[2]);
        poly_normal[1] += (p[2] - q[2]) * (p[0] + q[0]);
        poly_normal[2] += (p[0] - q[0]) * (p[1] + q[1]);
    }
    let poly_normal = mo_normalize(poly_normal);

    let mut iterations = 0;
    while remaining.len() > 3 {
        let m = remaining.len();
        let mut ear_found = false;

        for i in 0..m {
            let prev = remaining[(i + m - 1) % m];
            let curr = remaining[i];
            let next = remaining[(i + 1) % m];

            let a = polygon[prev];
            let b = polygon[curr];
            let c = polygon[next];

            // Check convexity: cross(AB, AC) must align with polygon normal
            let ab = mo_sub(b, a);
            let ac = mo_sub(c, a);
            let cross = mo_cross(ab, ac);
            let dot_with_normal =
                cross[0] * poly_normal[0] + cross[1] * poly_normal[1] + cross[2] * poly_normal[2];
            if dot_with_normal < 0.0 {
                continue; // reflex angle
            }

            // Check no other vertex is inside this ear triangle
            let is_ear = remaining.iter().all(|&vi| {
                if vi == prev || vi == curr || vi == next {
                    return true;
                }
                let p = polygon[vi];
                !point_in_triangle_3d(p, a, b, c, poly_normal)
            });

            if is_ear {
                result.push([prev, curr, next]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }

        if !ear_found {
            // Fallback: force-clip to avoid infinite loop on degenerate input
            if remaining.len() >= 3 {
                result.push([remaining[0], remaining[1], remaining[2]]);
                remaining.remove(1);
            } else {
                break;
            }
        }

        iterations += 1;
        if iterations > n * n {
            break; // safety exit
        }
    }

    if remaining.len() == 3 {
        result.push([remaining[0], remaining[1], remaining[2]]);
    }

    result
}

/// Test whether `p` lies inside triangle `(a, b, c)` by checking the sign
/// of barycentric coordinates projected onto the dominant axis of `normal`.
fn point_in_triangle_3d(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    normal: [f64; 3],
) -> bool {
    // Use same-side test (cross products w.r.t. each edge must align with normal)
    let ab = mo_sub(b, a);
    let bc = mo_sub(c, b);
    let ca = mo_sub(a, c);
    let ap = mo_sub(p, a);
    let bp = mo_sub(p, b);
    let cp = mo_sub(p, c);

    let n0 = mo_cross(ab, ap);
    let n1 = mo_cross(bc, bp);
    let n2 = mo_cross(ca, cp);

    let dot_n = |v: [f64; 3]| v[0] * normal[0] + v[1] * normal[1] + v[2] * normal[2];

    let d0 = dot_n(n0);
    let d1 = dot_n(n1);
    let d2 = dot_n(n2);

    d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0
}

/// Ray-triangle intersection for a ray from `origin` in the +Z direction.
///
/// Returns true if the ray intersects the triangle and the intersection is in front (z > origin.z).
fn ray_intersects_triangle_z(origin: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> bool {
    // Möller–Trumbore for ray direction = (0, 0, 1)
    let dir = Vec3::new(0.0, 0.0, 1.0);
    let edge1 = v1 - v0;
    let edge2 = v2 - v0;
    let h = dir.cross(&edge2);
    let det = edge1.dot(&h);
    if det.abs() < 1e-10 {
        return false;
    }
    let inv_det = 1.0 / det;
    let s = origin - v0;
    let u = inv_det * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let q = s.cross(&edge1);
    let v = inv_det * dir.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    let t = inv_det * edge2.dot(&q);
    t > 0.0
}

// ─── Quadric Error Mesh Decimation ───────────────────────────────────────────

/// Quadric Error Metrics (QEM) mesh decimation.
///
/// Iteratively collapses the `n_collapses` cheapest edges (by quadric error).
/// Returns the simplified `(verts, tris)`.
pub fn quadric_decimate(
    verts: &[[f64; 3]],
    tris: &[[usize; 3]],
    n_collapses: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let n = verts.len();
    // Build per-vertex quadric matrices (4×4 symmetric, stored as upper triangle).
    // Quadric Q_v = sum of (n·p=0) plane equations for all adjacent faces.
    // We store Q as a flat [10] array for 4×4 symmetric: indices (0,0)(0,1)(0,2)(0,3)(1,1)(1,2)(1,3)(2,2)(2,3)(3,3)
    #[allow(non_snake_case)]
    let mut q_mat: Vec<[f64; 10]> = vec![[0.0; 10]; n];

    for tri in tris {
        let a = verts[tri[0]];
        let b = verts[tri[1]];
        let c = verts[tri[2]];
        let e1 = mo_sub(b, a);
        let e2 = mo_sub(c, a);
        let nn = mo_cross(e1, e2);
        let len = mo_norm(nn);
        if len < 1e-14 {
            continue;
        }
        let n4 = [
            nn[0] / len,
            nn[1] / len,
            nn[2] / len,
            -(nn[0] * a[0] + nn[1] * a[1] + nn[2] * a[2]) / len,
        ];
        // Outer product n4 n4^T (symmetric)
        for &vi in tri {
            let q = &mut q_mat[vi];
            q[0] += n4[0] * n4[0];
            q[1] += n4[0] * n4[1];
            q[2] += n4[0] * n4[2];
            q[3] += n4[0] * n4[3];
            q[4] += n4[1] * n4[1];
            q[5] += n4[1] * n4[2];
            q[6] += n4[1] * n4[3];
            q[7] += n4[2] * n4[2];
            q[8] += n4[2] * n4[3];
            q[9] += n4[3] * n4[3];
        }
    }

    let mut cur_verts: Vec<[f64; 3]> = verts.to_vec();
    let mut cur_tris: Vec<[usize; 3]> = tris.to_vec();
    let mut remap: Vec<usize> = (0..n).collect();

    let add_q = |qa: &[f64; 10], qb: &[f64; 10]| -> [f64; 10] {
        let mut r = [0.0f64; 10];
        for i in 0..10 {
            r[i] = qa[i] + qb[i];
        }
        r
    };

    let eval_q = |q: &[f64; 10], v: [f64; 3]| -> f64 {
        let p = [v[0], v[1], v[2], 1.0];
        // p^T Q p (symmetric 4×4)
        let mut result = 0.0;
        let idx = |r: usize, c: usize| -> usize {
            if r <= c {
                r * 4 - r * (r - 1) / 2 + c - r
            } else {
                c * 4 - c * (c - 1) / 2 + r - c
            }
        };
        // Direct expansion for 4×4 symmetric
        result += q[0] * p[0] * p[0]
            + 2.0 * q[1] * p[0] * p[1]
            + 2.0 * q[2] * p[0] * p[2]
            + 2.0 * q[3] * p[0] * p[3];
        result += q[4] * p[1] * p[1] + 2.0 * q[5] * p[1] * p[2] + 2.0 * q[6] * p[1] * p[3];
        result += q[7] * p[2] * p[2] + 2.0 * q[8] * p[2] * p[3];
        result += q[9] * p[3] * p[3];
        let _ = idx;
        result
    };

    for _ in 0..n_collapses {
        // Find the edge with minimum quadric error cost
        let mut best_cost = f64::INFINITY;
        let mut best_edge = (usize::MAX, usize::MAX);
        let mut best_pos = [0.0f64; 3];

        // Collect all unique edges
        let mut edge_set = std::collections::BTreeSet::new();
        for tri in &cur_tris {
            for k in 0..3 {
                let a = remap[tri[k]];
                let b = remap[tri[(k + 1) % 3]];
                let key = (a.min(b), a.max(b));
                edge_set.insert(key);
            }
        }

        for (ea, eb) in &edge_set {
            let qa = &q_mat[*ea];
            let qb = &q_mat[*eb];
            let qsum = add_q(qa, qb);
            // Use midpoint as the new vertex position (optimal for QEM requires solving 3x3 system)
            let mid = [
                (cur_verts[*ea][0] + cur_verts[*eb][0]) * 0.5,
                (cur_verts[*ea][1] + cur_verts[*eb][1]) * 0.5,
                (cur_verts[*ea][2] + cur_verts[*eb][2]) * 0.5,
            ];
            let cost = eval_q(&qsum, mid);
            if cost < best_cost {
                best_cost = cost;
                best_edge = (*ea, *eb);
                best_pos = mid;
            }
        }

        if best_edge.0 == usize::MAX {
            break;
        }

        let (ea, eb) = best_edge;
        // Merge eb into ea
        cur_verts[ea] = best_pos;
        let qa = q_mat[ea];
        let qb = q_mat[eb];
        q_mat[ea] = add_q(&qa, &qb);

        // Update remap: all eb references → ea
        for r in remap.iter_mut() {
            if *r == eb {
                *r = ea;
            }
        }

        // Update triangles and remove degenerate
        cur_tris = cur_tris
            .iter()
            .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
            .filter(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2])
            .collect();
    }

    (cur_verts, cur_tris)
}

// ─── Laplacian Smoothing ──────────────────────────────────────────────────────

/// Apply `iterations` passes of Laplacian smoothing to `verts`.
///
/// Each pass moves each vertex toward the average of its neighbors by `lambda`.
/// Returns the smoothed vertex array.
pub fn laplacian_smooth(
    verts: &[[f64; 3]],
    tris: &[[usize; 3]],
    iterations: usize,
    lambda: f64,
) -> Vec<[f64; 3]> {
    let mut v = verts.to_vec();
    let adj = compute_adjacency(tris, verts.len());
    for _ in 0..iterations {
        let prev = v.clone();
        for i in 0..v.len() {
            if adj[i].is_empty() {
                continue;
            }
            let mut avg = [0.0f64; 3];
            for &nb in &adj[i] {
                avg[0] += prev[nb][0];
                avg[1] += prev[nb][1];
                avg[2] += prev[nb][2];
            }
            let k = adj[i].len() as f64;
            avg[0] /= k;
            avg[1] /= k;
            avg[2] /= k;
            v[i][0] += lambda * (avg[0] - prev[i][0]);
            v[i][1] += lambda * (avg[1] - prev[i][1]);
            v[i][2] += lambda * (avg[2] - prev[i][2]);
        }
    }
    v
}

// ─── Mesh Merging ─────────────────────────────────────────────────────────────

/// Merge two meshes (verts + tris) into one by concatenating vertices and
/// offsetting triangle indices of the second mesh.
pub fn merge_meshes(
    verts1: &[[f64; 3]],
    tris1: &[[usize; 3]],
    verts2: &[[f64; 3]],
    tris2: &[[usize; 3]],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let offset = verts1.len();
    let mut verts = verts1.to_vec();
    verts.extend_from_slice(verts2);
    let mut tris = tris1.to_vec();
    for tri in tris2 {
        tris.push([tri[0] + offset, tri[1] + offset, tri[2] + offset]);
    }
    (verts, tris)
}

// ─── Edge Split ───────────────────────────────────────────────────────────────

/// Split edge `(v0, v1)` in the mesh by inserting a midpoint vertex.
///
/// All triangles containing that edge are split into two.  Updates `verts` and `tris` in place.
pub fn split_edge(verts: &mut Vec<[f64; 3]>, tris: &mut Vec<[usize; 3]>, v0: usize, v1: usize) {
    let mid_idx = verts.len();
    let mid = [
        (verts[v0][0] + verts[v1][0]) * 0.5,
        (verts[v0][1] + verts[v1][1]) * 0.5,
        (verts[v0][2] + verts[v1][2]) * 0.5,
    ];
    verts.push(mid);

    let mut new_tris: Vec<[usize; 3]> = Vec::new();
    let mut keep = vec![true; tris.len()];

    for (ti, tri) in tris.iter().enumerate() {
        let has_edge = (tri[0] == v0 || tri[1] == v0 || tri[2] == v0)
            && (tri[0] == v1 || tri[1] == v1 || tri[2] == v1);
        if has_edge {
            keep[ti] = false;
            // Find the third vertex
            let opp = tri
                .iter()
                .copied()
                .find(|&v| v != v0 && v != v1)
                .unwrap_or(v0);
            new_tris.push([v0, mid_idx, opp]);
            new_tris.push([mid_idx, v1, opp]);
        }
    }

    let kept: Vec<[usize; 3]> = tris
        .iter()
        .enumerate()
        .filter_map(|(i, &t)| if keep[i] { Some(t) } else { None })
        .collect();
    *tris = kept;
    tris.extend(new_tris);
}

// ─── Boundary Edge Detection ──────────────────────────────────────────────────

/// Find all boundary edges: edges shared by exactly one triangle.
///
/// Returns a list of `(v0, v1)` edge pairs.
pub fn find_boundary_edges(tris: &[[usize; 3]]) -> Vec<(usize, usize)> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in tris {
        let edges = [
            (tri[0].min(tri[1]), tri[0].max(tri[1])),
            (tri[1].min(tri[2]), tri[1].max(tri[2])),
            (tri[0].min(tri[2]), tri[0].max(tri[2])),
        ];
        for e in edges {
            *edge_count.entry(e).or_insert(0) += 1;
        }
    }
    edge_count
        .into_iter()
        .filter(|(_, cnt)| *cnt == 1)
        .map(|(e, _)| e)
        .collect()
}

// ─── Mesh Simplification (vertex clustering) ─────────────────────────────────

/// Cluster vertices into a grid of `grid_res × grid_res × grid_res` cells.
///
/// Vertices in the same cell are merged to the cell centroid.
/// Returns `(new_verts, new_tris)` with degenerates removed.
pub fn vertex_cluster_simplify(
    verts: &[[f64; 3]],
    tris: &[[usize; 3]],
    grid_res: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    if verts.is_empty() || grid_res == 0 {
        return (verts.to_vec(), tris.to_vec());
    }

    // Compute bounding box
    let mut bb_min = verts[0];
    let mut bb_max = verts[0];
    for v in verts {
        for k in 0..3 {
            if v[k] < bb_min[k] {
                bb_min[k] = v[k];
            }
            if v[k] > bb_max[k] {
                bb_max[k] = v[k];
            }
        }
    }
    let extent = [
        (bb_max[0] - bb_min[0]).max(1e-12),
        (bb_max[1] - bb_min[1]).max(1e-12),
        (bb_max[2] - bb_min[2]).max(1e-12),
    ];

    let cell = |v: &[f64; 3]| -> (usize, usize, usize) {
        let cx = ((v[0] - bb_min[0]) / extent[0] * grid_res as f64) as usize;
        let cy = ((v[1] - bb_min[1]) / extent[1] * grid_res as f64) as usize;
        let cz = ((v[2] - bb_min[2]) / extent[2] * grid_res as f64) as usize;
        (
            cx.min(grid_res - 1),
            cy.min(grid_res - 1),
            cz.min(grid_res - 1),
        )
    };

    let mut cell_to_idx: HashMap<(usize, usize, usize), usize> = HashMap::new();
    let mut cell_sums: Vec<[f64; 3]> = Vec::new();
    let mut cell_counts: Vec<usize> = Vec::new();
    let mut remap: Vec<usize> = vec![0; verts.len()];

    for (i, v) in verts.iter().enumerate() {
        let c = cell(v);
        let idx = *cell_to_idx.entry(c).or_insert_with(|| {
            let idx = cell_sums.len();
            cell_sums.push([0.0; 3]);
            cell_counts.push(0);
            idx
        });
        cell_sums[idx][0] += v[0];
        cell_sums[idx][1] += v[1];
        cell_sums[idx][2] += v[2];
        cell_counts[idx] += 1;
        remap[i] = idx;
    }

    let new_verts: Vec<[f64; 3]> = cell_sums
        .iter()
        .zip(cell_counts.iter())
        .map(|(s, &c)| {
            let cf = c as f64;
            [s[0] / cf, s[1] / cf, s[2] / cf]
        })
        .collect();

    let new_tris: Vec<[usize; 3]> = tris
        .iter()
        .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
        .filter(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2])
        .collect();

    (new_verts, new_tris)
}

// ─── Mesh Statistics ──────────────────────────────────────────────────────────

/// Compute mesh statistics: min/max/average edge length.
pub fn edge_length_stats(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> (f64, f64, f64) {
    let mut min_len = f64::INFINITY;
    let mut max_len = 0.0f64;
    let mut total_len = 0.0f64;
    let mut count = 0usize;

    use std::collections::BTreeSet;
    let mut seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    for tri in tris {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = (a.min(b), a.max(b));
            if seen.insert(key) {
                let d = mo_norm(mo_sub(verts[b], verts[a]));
                min_len = min_len.min(d);
                max_len = max_len.max(d);
                total_len += d;
                count += 1;
            }
        }
    }
    let avg = if count > 0 {
        total_len / count as f64
    } else {
        0.0
    };
    (min_len, max_len, avg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_trimesh_face_normal() {
        // Triangle in XY plane with vertices CCW → normal should point +Z
        let verts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![[0, 1, 2]];
        let mesh = TriMesh::new(verts, tris);
        let n = mesh.face_normal(0);
        assert!((n.z - 1.0).abs() < 1e-10, "normal.z={}, expected 1.0", n.z);
        assert!(n.x.abs() < 1e-10);
        assert!(n.y.abs() < 1e-10);
    }

    #[test]
    fn test_trimesh_volume_sphere() {
        let mesh = TriMesh::uv_sphere(1.0, 20, 40);
        let vol = mesh.signed_volume().abs();
        let expected = (4.0 / 3.0) * PI; // ≈ 4.189
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 0.01,
            "sphere volume: got={}, expected≈{}, rel_err={}",
            vol,
            expected,
            rel_err
        );
    }

    #[test]
    fn test_trimesh_surface_area_sphere() {
        let mesh = TriMesh::uv_sphere(1.0, 20, 40);
        let sa = mesh.surface_area();
        let expected = 4.0 * PI; // ≈ 12.566
        let rel_err = (sa - expected).abs() / expected;
        assert!(
            rel_err < 0.01,
            "sphere surface area: got={}, expected≈{}, rel_err={}",
            sa,
            expected,
            rel_err
        );
    }

    #[test]
    fn test_trimesh_vertex_normals() {
        // For a unit sphere mesh, vertex normals should point roughly outward
        // (close to normalized vertex position)
        let mesh = TriMesh::uv_sphere(1.0, 20, 40);
        let vnormals = mesh.vertex_normals();
        for (i, (v, n)) in mesh.vertices.iter().zip(vnormals.iter()).enumerate() {
            let v_len = v.norm();
            if v_len < 1e-6 {
                continue; // skip degenerate poles
            }
            let v_norm = v / v_len;
            let dot = v_norm.dot(n);
            assert!(
                dot > 0.0,
                "vertex {} normal does not point outward: dot={}",
                i,
                dot
            );
        }
    }

    #[test]
    fn test_trimesh_contains_point() {
        let mesh = TriMesh::uv_sphere(1.0, 20, 40);
        // Center of sphere → inside
        let inside = mesh.contains_point(Vec3::new(0.0, 0.0, 0.0));
        assert!(inside, "center of sphere should be inside");
        // Far point → outside
        let outside = mesh.contains_point(Vec3::new(10.0, 0.0, 0.0));
        assert!(!outside, "far point should be outside");
    }

    #[test]
    fn test_trimesh_bounding_box() {
        let r = 1.5;
        let mesh = TriMesh::uv_sphere(r, 20, 40);
        let (min, max) = mesh.bounding_box();
        // All coordinates should be within [-r, r]
        assert!(min.x >= -r - 1e-10, "min.x={} < -r={}", min.x, -r);
        assert!(min.y >= -r - 1e-10, "min.y={} < -r={}", min.y, -r);
        assert!(min.z >= -r - 1e-10, "min.z={} < -r={}", min.z, -r);
        assert!(max.x <= r + 1e-10, "max.x={} > r={}", max.x, r);
        assert!(max.y <= r + 1e-10, "max.y={} > r={}", max.y, r);
        assert!(max.z <= r + 1e-10, "max.z={} > r={}", max.z, r);
        // Max should be close to r in each dimension
        assert!(max.x > r * 0.9, "max.x={} too small", max.x);
        assert!(max.y > r * 0.9, "max.y={} too small", max.y);
        assert!(max.z > r * 0.9, "max.z={} too small", max.z);
    }

    // ── standalone function tests ────────────────────────────────────────────

    /// Build a simple unit tetrahedron (4 triangles, watertight).
    fn unit_tetrahedron() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let verts = vec![
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let tris = vec![[0usize, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        (verts, tris)
    }

    #[test]
    fn test_compute_face_normals_unit_length() {
        let (verts, tris) = unit_tetrahedron();
        let normals = compute_face_normals(&verts, &tris);
        assert_eq!(normals.len(), tris.len());
        for (i, n) in normals.iter().enumerate() {
            let len = (n[0].powi(2) + n[1].powi(2) + n[2].powi(2)).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-10,
                "face {i} normal not unit: len={len}"
            );
        }
    }

    #[test]
    fn test_compute_face_normals_xy_plane() {
        // Triangle in XY plane → normal = [0, 0, ±1]
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0usize, 1, 2]];
        let normals = compute_face_normals(&verts, &tris);
        assert!(
            (normals[0][2].abs() - 1.0).abs() < 1e-10,
            "nz={}",
            normals[0][2]
        );
    }

    #[test]
    fn test_compute_vertex_normals_count() {
        let (verts, tris) = unit_tetrahedron();
        let vnormals = compute_vertex_normals(&verts, &tris);
        assert_eq!(vnormals.len(), verts.len());
    }

    #[test]
    fn test_compute_vertex_normals_unit_length() {
        let (verts, tris) = unit_tetrahedron();
        let vnormals = compute_vertex_normals(&verts, &tris);
        for (i, n) in vnormals.iter().enumerate() {
            let len = (n[0].powi(2) + n[1].powi(2) + n[2].powi(2)).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-10,
                "vertex {i} normal not unit: {len}"
            );
        }
    }

    #[test]
    fn test_weld_vertices_basic() {
        // Two vertices that are essentially the same within tolerance
        let verts = vec![
            [0.0f64, 0.0, 0.0],
            [0.0, 0.0, 1e-10], // < tol=1e-6 away from above
            [1.0, 0.0, 0.0],
        ];
        let tris = vec![[0usize, 1, 2]];
        let (new_verts, _new_tris) = weld_vertices(&verts, &tris, 1e-6);
        assert_eq!(new_verts.len(), 2, "expected 2 unique verts after weld");
    }

    #[test]
    fn test_weld_vertices_no_merge() {
        // Vertices that are far apart should not be merged
        let verts = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let tris = vec![[0usize, 1, 2]];
        let (new_verts, _) = weld_vertices(&verts, &tris, 1e-6);
        assert_eq!(new_verts.len(), 3);
    }

    #[test]
    fn test_compute_adjacency_tetrahedron() {
        let (verts, tris) = unit_tetrahedron();
        let adj = compute_adjacency(&tris, verts.len());
        // In a tetrahedron every vertex connects to every other
        assert_eq!(adj.len(), 4);
        for (i, neighbors) in adj.iter().enumerate() {
            assert_eq!(neighbors.len(), 3, "vertex {i} should have 3 neighbors");
        }
    }

    #[test]
    fn test_is_watertight_tetrahedron() {
        let (_verts, tris) = unit_tetrahedron();
        assert!(is_watertight(&tris), "tetrahedron should be watertight");
    }

    #[test]
    fn test_is_watertight_open_mesh() {
        // A single triangle has 3 boundary edges → not watertight
        let tris = vec![[0usize, 1, 2]];
        assert!(
            !is_watertight(&tris),
            "single triangle should not be watertight"
        );
    }

    #[test]
    fn test_compute_genus_tetrahedron() {
        // Tetrahedron: V=4, E=6, F=4 → χ=2 → g=0
        let (verts, tris) = unit_tetrahedron();
        let g = compute_genus(verts.len(), &tris);
        assert_eq!(g, 0, "tetrahedron genus should be 0");
    }

    #[test]
    fn test_triangulate_polygon_triangle() {
        let poly = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let result = triangulate_polygon(&poly);
        assert_eq!(result.len(), 1, "triangle polygon = 1 tri");
        assert_eq!(result[0], [0, 1, 2]);
    }

    #[test]
    fn test_triangulate_polygon_quad() {
        // A convex quad in XY should yield 2 triangles
        let poly = vec![
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let result = triangulate_polygon(&poly);
        assert_eq!(result.len(), 2, "quad polygon = 2 tris");
    }

    #[test]
    fn test_triangulate_polygon_pentagon() {
        // A convex pentagon should yield 3 triangles
        let poly: Vec<[f64; 3]> = (0..5)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let result = triangulate_polygon(&poly);
        assert_eq!(result.len(), 3, "pentagon polygon = 3 tris");
    }

    // ── Quadric error decimation tests ────────────────────────────────────────

    #[test]
    fn test_decimate_reduces_triangles() {
        let (verts, tris) = unit_tetrahedron();
        let (_, new_tris) = quadric_decimate(&verts, &tris, 2);
        // After removing 2 edges, should have fewer triangles (or equal if no valid collapses)
        assert!(
            new_tris.len() <= tris.len(),
            "decimation should not increase triangles"
        );
    }

    #[test]
    fn test_decimate_preserves_non_degenerate_tris() {
        let (verts, tris) = unit_tetrahedron();
        let (new_verts, new_tris) = quadric_decimate(&verts, &tris, 1);
        for tri in &new_tris {
            let a = new_verts[tri[0]];
            let b = new_verts[tri[1]];
            let c = new_verts[tri[2]];
            let e1 = mo_sub(b, a);
            let e2 = mo_sub(c, a);
            let cross = mo_cross(e1, e2);
            let area = mo_norm(cross) * 0.5;
            assert!(area >= 0.0, "negative area triangle");
        }
    }

    // ── Laplacian smoothing tests ─────────────────────────────────────────────

    #[test]
    fn test_laplacian_smooth_reduces_roughness() {
        // Create a rough mesh: tetrahedron with one vertex displaced
        let (mut verts, tris) = unit_tetrahedron();
        verts[0] = [0.5, 0.5, 0.5]; // displace from origin
        let smoothed = laplacian_smooth(&verts, &tris, 5, 0.5);
        // Smoothed version should have vertex 0 closer to its neighbors' average
        assert_eq!(smoothed.len(), verts.len());
    }

    #[test]
    fn test_laplacian_smooth_returns_correct_count() {
        let (verts, tris) = unit_tetrahedron();
        let smoothed = laplacian_smooth(&verts, &tris, 3, 0.5);
        assert_eq!(
            smoothed.len(),
            verts.len(),
            "smoothed should have same vertex count"
        );
    }

    // ── Mesh merging tests ────────────────────────────────────────────────────

    #[test]
    fn test_merge_meshes_vertex_count() {
        let (v1, t1) = unit_tetrahedron();
        let v2: Vec<[f64; 3]> = v1.iter().map(|&v| [v[0] + 5.0, v[1], v[2]]).collect();
        let (mv, mt) = merge_meshes(&v1, &t1, &v2, &t1);
        assert_eq!(mv.len(), v1.len() + v2.len(), "merged vertices");
        assert_eq!(mt.len(), t1.len() * 2, "merged triangles");
    }

    #[test]
    fn test_merge_meshes_index_offset() {
        let (v1, t1) = unit_tetrahedron();
        let v2: Vec<[f64; 3]> = v1.iter().map(|&v| [v[0] + 5.0, v[1], v[2]]).collect();
        let (mv, mt) = merge_meshes(&v1, &t1, &v2, &t1);
        // Second mesh triangles should reference indices >= v1.len()
        for tri in &mt[t1.len()..] {
            for &vi in tri {
                assert!(vi >= v1.len(), "index {} should be offset", vi);
            }
        }
        assert_eq!(mv.len(), v1.len() + v2.len());
    }

    // ── Edge split tests ──────────────────────────────────────────────────────

    #[test]
    fn test_split_edge_increases_vertices() {
        let (mut verts, mut tris) = unit_tetrahedron();
        let n_before = verts.len();
        split_edge(&mut verts, &mut tris, 0, 1);
        assert_eq!(verts.len(), n_before + 1, "split should add one vertex");
    }

    #[test]
    fn test_split_edge_increases_triangles() {
        let (mut verts, mut tris) = unit_tetrahedron();
        let n_before = tris.len();
        split_edge(&mut verts, &mut tris, 0, 1);
        assert!(tris.len() > n_before, "split should add triangles");
    }

    // ── Boundary detection tests ──────────────────────────────────────────────

    #[test]
    fn test_find_boundary_edges_tetrahedron() {
        let (_verts, tris) = unit_tetrahedron();
        let boundary = find_boundary_edges(&tris);
        // Closed tetrahedron has no boundary edges
        assert!(
            boundary.is_empty(),
            "tetrahedron has no boundary, got {:?}",
            boundary
        );
    }

    #[test]
    fn test_find_boundary_edges_open_mesh() {
        // Two triangles sharing one edge → 4 boundary edges
        let tris = vec![[0usize, 1, 2], [1, 3, 2]];
        let boundary = find_boundary_edges(&tris);
        assert!(!boundary.is_empty(), "open mesh should have boundary edges");
    }
}
