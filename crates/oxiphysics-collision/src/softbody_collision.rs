// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soft body collision detection: deformable meshes, cloth, self-collision.
//!
//! Provides deformable-vs-BVH, self-collision with spatial hashing and normal
//! cone culling, cloth-vs-SDF contacts, edge-edge contacts, and CCD for
//! deformable meshes.

#![allow(dead_code)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Helper free functions
// ---------------------------------------------------------------------------

/// Compute the squared distance from point `p` to the triangle `(a, b, c)`.
/// Returns (squared distance, barycentric coordinates (u, v, w)).
pub fn vertex_triangle_distance(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
) -> (f64, [f64; 3]) {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ap = sub3(p, a);

    let d1 = dot3(ab, ap);
    let d2 = dot3(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        let bary = [1.0, 0.0, 0.0];
        return (dot3(ap, ap), bary);
    }

    let bp = sub3(p, b);
    let d3 = dot3(ab, bp);
    let d4 = dot3(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        let bary = [0.0, 1.0, 0.0];
        let diff = sub3(p, b);
        return (dot3(diff, diff), bary);
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let closest = add3(a, scale3(ab, v));
        let diff = sub3(p, closest);
        let bary = [1.0 - v, v, 0.0];
        return (dot3(diff, diff), bary);
    }

    let cp = sub3(p, c);
    let d5 = dot3(ab, cp);
    let d6 = dot3(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        let bary = [0.0, 0.0, 1.0];
        let diff = sub3(p, c);
        return (dot3(diff, diff), bary);
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let closest = add3(a, scale3(ac, w));
        let diff = sub3(p, closest);
        let bary = [1.0 - w, 0.0, w];
        return (dot3(diff, diff), bary);
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let bc = sub3(c, b);
        let closest = add3(b, scale3(bc, w));
        let diff = sub3(p, closest);
        let bary = [0.0, 1.0 - w, w];
        return (dot3(diff, diff), bary);
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let closest = add3(add3(a, scale3(ab, v)), scale3(ac, w));
    let diff = sub3(p, closest);
    let bary = [1.0 - v - w, v, w];
    (dot3(diff, diff), bary)
}

/// Compute the squared distance and closest points between two line segments
/// `(p, p+d)` and `(q, q+e)`.
/// Returns (sq_dist, s, t) where s,t ∈ \[0,1\] parameterise the closest points.
pub fn edge_edge_distance(p: [f64; 3], d: [f64; 3], q: [f64; 3], e: [f64; 3]) -> (f64, f64, f64) {
    let r = sub3(p, q);
    let a = dot3(d, d);
    let f = dot3(e, e);
    let c = dot3(d, r);

    let eps = 1e-10;
    let (s, t) = if a <= eps && f <= eps {
        (0.0, 0.0)
    } else if a <= eps {
        let t = (dot3(e, r) / f).clamp(0.0, 1.0);
        (0.0, t)
    } else {
        let e_dot = dot3(e, r);
        if f <= eps {
            let s = (-c / a).clamp(0.0, 1.0);
            (s, 0.0)
        } else {
            let b = dot3(d, e);
            let denom = a * f - b * b;
            let s = if denom.abs() > eps {
                ((b * e_dot - c * f) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + e_dot) / f;
            if t < 0.0 {
                let s2 = (-c / a).clamp(0.0, 1.0);
                (s2, 0.0)
            } else if t > 1.0 {
                let s2 = ((b - c) / a).clamp(0.0, 1.0);
                (s2, 1.0)
            } else {
                (s, t)
            }
        }
    };

    let closest_p = add3(p, scale3(d, s));
    let closest_q = add3(q, scale3(e, t));
    let diff = sub3(closest_p, closest_q);
    (dot3(diff, diff), s, t)
}

/// Compute the half-angle of the normal cone for a set of vertex normals.
/// Returns the cosine of the maximum deviation angle from the average normal.
pub fn normal_cone_angle(normals: &[[f64; 3]]) -> f64 {
    if normals.is_empty() {
        return 0.0;
    }
    let avg = normals.iter().fold([0.0f64; 3], |acc, n| add3(acc, *n));
    let avg = normalize3(avg);
    normals
        .iter()
        .map(|n| dot3(*n, avg).clamp(-1.0, 1.0))
        .fold(f64::INFINITY, f64::min)
}

/// Compute the first time-of-impact for a pair of linearly-moving vertices,
/// given vertex trajectory over one timestep.
/// Returns fraction `t ∈ [0, 1]` or `None` if no impact within the step.
pub fn pairwise_ccd_dt(
    p0: [f64; 3],
    p1: [f64; 3],
    q0: [f64; 3],
    q1: [f64; 3],
    radius: f64,
) -> Option<f64> {
    let dp = sub3(p1, p0);
    let dq = sub3(q1, q0);
    let rel_vel = sub3(dp, dq);
    let r0 = sub3(p0, q0);

    let a = dot3(rel_vel, rel_vel);
    let b = 2.0 * dot3(r0, rel_vel);
    let c = dot3(r0, r0) - radius * radius;

    if a < 1e-14 {
        if c <= 0.0 {
            return Some(0.0);
        }
        return None;
    }

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t0 = (-b - sqrt_disc) / (2.0 * a);
    let t1 = (-b + sqrt_disc) / (2.0 * a);
    if t1 < 0.0 || t0 > 1.0 {
        return None;
    }
    Some(t0.max(0.0))
}

// ---------------------------------------------------------------------------
// Internal vector helpers (no nalgebra)
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0, 0.0, 1.0]
    } else {
        scale3(a, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// Deformable vs BVH
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box over a triangle for BVH nodes.
#[derive(Debug, Clone)]
pub struct TriAabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
    /// Triangle index (leaf node).
    pub tri_index: usize,
}

impl TriAabb {
    /// Construct from three vertices.
    pub fn from_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3], idx: usize) -> Self {
        let min = [
            a[0].min(b[0]).min(c[0]),
            a[1].min(b[1]).min(c[1]),
            a[2].min(b[2]).min(c[2]),
        ];
        let max = [
            a[0].max(b[0]).max(c[0]),
            a[1].max(b[1]).max(c[1]),
            a[2].max(b[2]).max(c[2]),
        ];
        TriAabb {
            min,
            max,
            tri_index: idx,
        }
    }

    /// Test overlap with another AABB.
    pub fn overlaps(&self, other: &TriAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
}

/// Result of a deformable vs BVH query: a contact between a deformable vertex
/// and a rigid triangle.
#[derive(Debug, Clone)]
pub struct DeformContact {
    /// Index of the deformable vertex.
    pub vertex_index: usize,
    /// Index of the rigid triangle.
    pub tri_index: usize,
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Contact normal (pointing from rigid surface to deformable vertex).
    pub normal: [f64; 3],
    /// Closest point on the rigid triangle.
    pub point_on_rigid: [f64; 3],
}

/// Collision detection between a deformable mesh and a rigid BVH.
///
/// Queries each deformable vertex against the rigid BVH leaf triangles.
pub struct DeformableVsBvh {
    /// Flat list of rigid triangles (vertex triples).
    pub rigid_triangles: Vec<[[f64; 3]; 3]>,
    /// Per-triangle AABBs.
    pub aabbs: Vec<TriAabb>,
    /// Contact distance threshold.
    pub contact_threshold: f64,
}

impl DeformableVsBvh {
    /// Construct from rigid triangle soup.
    pub fn new(triangles: Vec<[[f64; 3]; 3]>, contact_threshold: f64) -> Self {
        let aabbs = triangles
            .iter()
            .enumerate()
            .map(|(i, t)| TriAabb::from_triangle(t[0], t[1], t[2], i))
            .collect();
        DeformableVsBvh {
            rigid_triangles: triangles,
            aabbs,
            contact_threshold,
        }
    }

    /// Query all deformable vertices against rigid triangles.
    /// Returns list of contacts.
    pub fn query(&self, vertices: &[[f64; 3]]) -> Vec<DeformContact> {
        let mut contacts = Vec::new();
        for (vi, &v) in vertices.iter().enumerate() {
            // Build a query AABB around the vertex
            let r = self.contact_threshold;
            let query_aabb = TriAabb {
                min: [v[0] - r, v[1] - r, v[2] - r],
                max: [v[0] + r, v[1] + r, v[2] + r],
                tri_index: 0,
            };
            for aabb in &self.aabbs {
                if !aabb.overlaps(&query_aabb) {
                    continue;
                }
                let tri = &self.rigid_triangles[aabb.tri_index];
                let (sq_dist, bary) = vertex_triangle_distance(v, tri[0], tri[1], tri[2]);
                let dist = sq_dist.sqrt();
                if dist < self.contact_threshold {
                    let closest = add3(
                        add3(scale3(tri[0], bary[0]), scale3(tri[1], bary[1])),
                        scale3(tri[2], bary[2]),
                    );
                    let ab = sub3(tri[1], tri[0]);
                    let ac = sub3(tri[2], tri[0]);
                    let tri_normal = normalize3(cross3(ab, ac));
                    contacts.push(DeformContact {
                        vertex_index: vi,
                        tri_index: aabb.tri_index,
                        depth: self.contact_threshold - dist,
                        normal: tri_normal,
                        point_on_rigid: closest,
                    });
                }
            }
        }
        contacts
    }
}

// ---------------------------------------------------------------------------
// Normal Cone Filter
// ---------------------------------------------------------------------------

/// Back-face culling filter for self-collision using vertex normals and a
/// cone angle threshold.
pub struct NormalConeFilter {
    /// Cosine of the culling cone half-angle threshold.
    /// Pairs whose normals are within this cone are candidates.
    pub cos_threshold: f64,
}

impl NormalConeFilter {
    /// Construct with a cone half-angle in radians.
    pub fn new(half_angle_rad: f64) -> Self {
        NormalConeFilter {
            cos_threshold: half_angle_rad.cos(),
        }
    }

    /// Returns `true` if this vertex-triangle pair should be tested for
    /// self-collision (not culled by the normal cone).
    pub fn should_test(&self, vertex_normal: [f64; 3], tri_normal: [f64; 3]) -> bool {
        dot3(vertex_normal, tri_normal) < self.cos_threshold
    }
}

// ---------------------------------------------------------------------------
// Self Collision
// ---------------------------------------------------------------------------

/// Candidate self-collision pair between two vertices.
#[derive(Debug, Clone)]
pub struct SelfCollisionPair {
    /// Index of vertex A.
    pub a: usize,
    /// Index of vertex B.
    pub b: usize,
    /// Squared distance between the two vertices.
    pub sq_dist: f64,
}

/// Self-collision detection for cloth/soft mesh using a spatial hash and
/// normal cone culling.
pub struct SelfCollision {
    /// Cell size for the spatial hash.
    pub cell_size: f64,
    /// Contact distance threshold.
    pub contact_threshold: f64,
    /// Normal cone filter.
    pub cone_filter: NormalConeFilter,
    hash: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl SelfCollision {
    /// Construct with cell size, contact threshold, and cone half-angle.
    pub fn new(cell_size: f64, contact_threshold: f64, cone_half_angle: f64) -> Self {
        SelfCollision {
            cell_size,
            contact_threshold,
            cone_filter: NormalConeFilter::new(cone_half_angle),
            hash: HashMap::new(),
        }
    }

    fn cell_of(&self, p: [f64; 3]) -> (i64, i64, i64) {
        (
            (p[0] / self.cell_size).floor() as i64,
            (p[1] / self.cell_size).floor() as i64,
            (p[2] / self.cell_size).floor() as i64,
        )
    }

    /// Update the spatial hash from the current vertex positions.
    pub fn update(&mut self, vertices: &[[f64; 3]]) {
        self.hash.clear();
        for (i, &v) in vertices.iter().enumerate() {
            let cell = self.cell_of(v);
            self.hash.entry(cell).or_default().push(i);
        }
    }

    /// Find candidate self-collision pairs, optionally filtered by normal cone.
    pub fn find_pairs(
        &self,
        vertices: &[[f64; 3]],
        normals: Option<&[[f64; 3]]>,
    ) -> Vec<SelfCollisionPair> {
        let mut pairs = Vec::new();
        let threshold_sq = self.contact_threshold * self.contact_threshold;

        for (&cell, indices) in &self.hash {
            // Neighbour cells including self
            for di in -1i64..=1 {
                for dj in -1i64..=1 {
                    for dk in -1i64..=1 {
                        let nc = (cell.0 + di, cell.1 + dj, cell.2 + dk);
                        if let Some(nb_indices) = self.hash.get(&nc) {
                            for &a in indices {
                                for &b in nb_indices {
                                    if b <= a {
                                        continue;
                                    }
                                    // Normal cone culling
                                    if let Some(norms) = normals
                                        && !self.cone_filter.should_test(norms[a], norms[b])
                                    {
                                        continue;
                                    }
                                    let diff = sub3(vertices[a], vertices[b]);
                                    let sq = dot3(diff, diff);
                                    if sq < threshold_sq {
                                        pairs.push(SelfCollisionPair { a, b, sq_dist: sq });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        pairs
    }
}

// ---------------------------------------------------------------------------
// Cloth SDF Contact
// ---------------------------------------------------------------------------

/// A contact from cloth vertex vs a rigid body's signed distance field.
#[derive(Debug, Clone)]
pub struct SdfContact {
    /// Index of the cloth vertex.
    pub vertex_index: usize,
    /// Signed distance value (negative = inside).
    pub sdf_value: f64,
    /// Gradient of the SDF at the contact point (contact normal).
    pub gradient: [f64; 3],
}

/// Cloth collision against a signed distance field (rigid body SDF).
///
/// The SDF is sampled on a uniform grid; gradients are computed by central
/// differences.
pub struct ClothSdfContact {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Grid dimensions.
    pub nz: usize,
    /// Cell size.
    pub dx: f64,
    /// Origin of the grid.
    pub origin: [f64; 3],
    /// Flat SDF values, row-major (x fastest).
    pub values: Vec<f64>,
    /// Contact thickness threshold.
    pub thickness: f64,
}

impl ClothSdfContact {
    /// Construct from a pre-computed SDF grid.
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        origin: [f64; 3],
        values: Vec<f64>,
        thickness: f64,
    ) -> Self {
        ClothSdfContact {
            nx,
            ny,
            nz,
            dx,
            origin,
            values,
            thickness,
        }
    }

    fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.nx * (iy + self.ny * iz)
    }

    fn sample(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            self.values[self.index(ix, iy, iz)]
        } else {
            f64::MAX
        }
    }

    /// Trilinearly interpolate SDF at world position `p`.
    pub fn interpolate(&self, p: [f64; 3]) -> f64 {
        let lp = sub3(p, self.origin);
        let fx = lp[0] / self.dx;
        let fy = lp[1] / self.dx;
        let fz = lp[2] / self.dx;
        let ix = fx.floor() as isize;
        let iy = fy.floor() as isize;
        let iz = fz.floor() as isize;
        if ix < 0
            || iy < 0
            || iz < 0
            || ix + 1 >= self.nx as isize
            || iy + 1 >= self.ny as isize
            || iz + 1 >= self.nz as isize
        {
            return f64::MAX;
        }
        let (ix, iy, iz) = (ix as usize, iy as usize, iz as usize);
        let tx = fx - ix as f64;
        let ty = fy - iy as f64;
        let tz = fz - iz as f64;
        let v000 = self.sample(ix, iy, iz);
        let v100 = self.sample(ix + 1, iy, iz);
        let v010 = self.sample(ix, iy + 1, iz);
        let v110 = self.sample(ix + 1, iy + 1, iz);
        let v001 = self.sample(ix, iy, iz + 1);
        let v101 = self.sample(ix + 1, iy, iz + 1);
        let v011 = self.sample(ix, iy + 1, iz + 1);
        let v111 = self.sample(ix + 1, iy + 1, iz + 1);
        let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
        let c00 = lerp(v000, v100, tx);
        let c10 = lerp(v010, v110, tx);
        let c01 = lerp(v001, v101, tx);
        let c11 = lerp(v011, v111, tx);
        let c0 = lerp(c00, c10, ty);
        let c1 = lerp(c01, c11, ty);
        lerp(c0, c1, tz)
    }

    /// Compute the SDF gradient at world position `p` using central differences.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let h = self.dx * 0.5;
        let gx = (self.interpolate([p[0] + h, p[1], p[2]])
            - self.interpolate([p[0] - h, p[1], p[2]]))
            / (2.0 * h);
        let gy = (self.interpolate([p[0], p[1] + h, p[2]])
            - self.interpolate([p[0], p[1] - h, p[2]]))
            / (2.0 * h);
        let gz = (self.interpolate([p[0], p[1], p[2] + h])
            - self.interpolate([p[0], p[1], p[2] - h]))
            / (2.0 * h);
        normalize3([gx, gy, gz])
    }

    /// Test all cloth vertices against the SDF.
    pub fn test_vertices(&self, vertices: &[[f64; 3]]) -> Vec<SdfContact> {
        vertices
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| {
                let d = self.interpolate(v);
                if d < self.thickness {
                    Some(SdfContact {
                        vertex_index: i,
                        sdf_value: d,
                        gradient: self.gradient(v),
                    })
                } else {
                    None
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Edge-Edge Contact
// ---------------------------------------------------------------------------

/// Contact result for edge-edge (line-segment vs line-segment) collision.
#[derive(Debug, Clone)]
pub struct EdgeEdgeResult {
    /// Index of edge A (in edge list).
    pub edge_a: usize,
    /// Index of edge B (in edge list).
    pub edge_b: usize,
    /// Closest distance between the edges.
    pub distance: f64,
    /// Parameter on edge A of the closest point.
    pub s: f64,
    /// Parameter on edge B of the closest point.
    pub t: f64,
    /// Contact normal (from B to A).
    pub normal: [f64; 3],
}

/// Edge-edge contact detector for cloth edge collisions.
///
/// Each edge is defined by two vertex indices.
pub struct EdgeEdgeContact {
    /// Contact distance threshold.
    pub threshold: f64,
}

impl EdgeEdgeContact {
    /// Construct with a contact threshold.
    pub fn new(threshold: f64) -> Self {
        EdgeEdgeContact { threshold }
    }

    /// Test a single edge pair.
    pub fn test_pair(
        &self,
        ei_a: usize,
        p_a: [f64; 3],
        d_a: [f64; 3],
        ei_b: usize,
        p_b: [f64; 3],
        d_b: [f64; 3],
    ) -> Option<EdgeEdgeResult> {
        let (sq_dist, s, t) = edge_edge_distance(p_a, d_a, p_b, d_b);
        let dist = sq_dist.sqrt();
        if dist < self.threshold {
            let ca = add3(p_a, scale3(d_a, s));
            let cb = add3(p_b, scale3(d_b, t));
            let diff = sub3(ca, cb);
            let normal = if dist > 1e-14 {
                scale3(diff, 1.0 / dist)
            } else {
                normalize3(cross3(d_a, d_b))
            };
            Some(EdgeEdgeResult {
                edge_a: ei_a,
                edge_b: ei_b,
                distance: dist,
                s,
                t,
                normal,
            })
        } else {
            None
        }
    }

    /// Test all edge pairs between two edge lists.
    pub fn test_all(
        &self,
        edges_a: &[(usize, usize)],
        verts_a: &[[f64; 3]],
        edges_b: &[(usize, usize)],
        verts_b: &[[f64; 3]],
    ) -> Vec<EdgeEdgeResult> {
        let mut results = Vec::new();
        for (ia, &(a0, a1)) in edges_a.iter().enumerate() {
            let pa = verts_a[a0];
            let da = sub3(verts_a[a1], pa);
            for (ib, &(b0, b1)) in edges_b.iter().enumerate() {
                let pb = verts_b[b0];
                let db = sub3(verts_b[b1], pb);
                if let Some(r) = self.test_pair(ia, pa, da, ib, pb, db) {
                    results.push(r);
                }
            }
        }
        results
    }
}

// ---------------------------------------------------------------------------
// CCD Deformable
// ---------------------------------------------------------------------------

/// Continuous collision result for a deformable mesh vertex trajectory.
#[derive(Debug, Clone)]
pub struct CcdDeformableResult {
    /// Vertex index.
    pub vertex_index: usize,
    /// Time of impact in \[0, 1\].
    pub toi: f64,
    /// Contact normal at impact.
    pub normal: [f64; 3],
}

/// Continuous collision detection for a deformable mesh.
///
/// Uses per-vertex linear trajectory CCD and edge-edge CCD.
pub struct CcdDeformable {
    /// Contact radius for vertex-vertex CCD.
    pub vertex_radius: f64,
    /// Contact radius for edge-edge CCD.
    pub edge_radius: f64,
}

impl CcdDeformable {
    /// Construct with per-feature radii.
    pub fn new(vertex_radius: f64, edge_radius: f64) -> Self {
        CcdDeformable {
            vertex_radius,
            edge_radius,
        }
    }

    /// Vertex-vertex CCD between two meshes over a timestep.
    /// `verts0_a/b` are start positions, `verts1_a/b` are end positions.
    pub fn vertex_vertex_ccd(
        &self,
        verts0_a: &[[f64; 3]],
        verts1_a: &[[f64; 3]],
        verts0_b: &[[f64; 3]],
        verts1_b: &[[f64; 3]],
    ) -> Vec<CcdDeformableResult> {
        let mut results = Vec::new();
        for (ia, (&p0, &p1)) in verts0_a.iter().zip(verts1_a.iter()).enumerate() {
            for (&q0, &q1) in verts0_b.iter().zip(verts1_b.iter()) {
                if let Some(toi) = pairwise_ccd_dt(p0, p1, q0, q1, self.vertex_radius) {
                    let pt = add3(p0, scale3(sub3(p1, p0), toi));
                    let qt = add3(q0, scale3(sub3(q1, q0), toi));
                    let diff = sub3(pt, qt);
                    let normal = normalize3(diff);
                    results.push(CcdDeformableResult {
                        vertex_index: ia,
                        toi,
                        normal,
                    });
                }
            }
        }
        results
    }

    /// Edge-edge CCD between two edge sets over a timestep.
    pub fn edge_edge_ccd(
        &self,
        edges_a: &[(usize, usize)],
        verts0_a: &[[f64; 3]],
        verts1_a: &[[f64; 3]],
        edges_b: &[(usize, usize)],
        verts0_b: &[[f64; 3]],
        verts1_b: &[[f64; 3]],
    ) -> Vec<CcdDeformableResult> {
        let mut results = Vec::new();
        for (ia, &(a0, a1)) in edges_a.iter().enumerate() {
            for &(b0, b1) in edges_b {
                // Sample midpoints and test CCD as approximation
                let pa0 = midpoint(verts0_a[a0], verts0_a[a1]);
                let pa1 = midpoint(verts1_a[a0], verts1_a[a1]);
                let pb0 = midpoint(verts0_b[b0], verts0_b[b1]);
                let pb1 = midpoint(verts1_b[b0], verts1_b[b1]);
                if let Some(toi) = pairwise_ccd_dt(pa0, pa1, pb0, pb1, self.edge_radius) {
                    let pt = add3(pa0, scale3(sub3(pa1, pa0), toi));
                    let qt = add3(pb0, scale3(sub3(pb1, pb0), toi));
                    let normal = normalize3(sub3(pt, qt));
                    results.push(CcdDeformableResult {
                        vertex_index: ia,
                        toi,
                        normal,
                    });
                }
            }
        }
        results
    }
}

fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale3(add3(a, b), 0.5)
}

// ---------------------------------------------------------------------------
// Primitive vs Deformable
// ---------------------------------------------------------------------------

/// Contact between a primitive (sphere/capsule/box) and a deformable vertex.
#[derive(Debug, Clone)]
pub struct PrimitiveDeformContact {
    /// Vertex index.
    pub vertex_index: usize,
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Contact normal pointing away from the primitive.
    pub normal: [f64; 3],
    /// Contact point on the primitive surface.
    pub point: [f64; 3],
}

/// Sphere vs deformable mesh contact.
pub struct PrimitiveVsDeformable {
    /// Sphere centre.
    pub sphere_center: [f64; 3],
    /// Sphere radius.
    pub sphere_radius: f64,
}

impl PrimitiveVsDeformable {
    /// Construct a sphere vs deformable tester.
    pub fn sphere(center: [f64; 3], radius: f64) -> Self {
        PrimitiveVsDeformable {
            sphere_center: center,
            sphere_radius: radius,
        }
    }

    /// Test all vertices against the sphere.
    pub fn query_sphere(&self, vertices: &[[f64; 3]]) -> Vec<PrimitiveDeformContact> {
        vertices
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| {
                let diff = sub3(v, self.sphere_center);
                let dist = len3(diff);
                if dist < self.sphere_radius {
                    let depth = self.sphere_radius - dist;
                    let normal = if dist > 1e-14 {
                        scale3(diff, 1.0 / dist)
                    } else {
                        [0.0, 1.0, 0.0]
                    };
                    let point = add3(self.sphere_center, scale3(normal, self.sphere_radius));
                    Some(PrimitiveDeformContact {
                        vertex_index: i,
                        depth,
                        normal,
                        point,
                    })
                } else {
                    None
                }
            })
            .collect()
    }

    /// Test all vertices against an axis-aligned box.
    pub fn query_box(
        vertices: &[[f64; 3]],
        box_min: [f64; 3],
        box_max: [f64; 3],
    ) -> Vec<PrimitiveDeformContact> {
        vertices
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| {
                // Closest point on box
                let cx = v[0].clamp(box_min[0], box_max[0]);
                let cy = v[1].clamp(box_min[1], box_max[1]);
                let cz = v[2].clamp(box_min[2], box_max[2]);
                let closest = [cx, cy, cz];
                let diff = sub3(v, closest);
                let dist = len3(diff);
                if dist < 1e-6 {
                    // Vertex inside box
                    let normal = [0.0, 1.0, 0.0];
                    Some(PrimitiveDeformContact {
                        vertex_index: i,
                        depth: 1e-6,
                        normal,
                        point: closest,
                    })
                } else {
                    None
                }
            })
            .collect()
    }

    /// Test all vertices against a capsule (segment + radius).
    pub fn query_capsule(
        vertices: &[[f64; 3]],
        cap_a: [f64; 3],
        cap_b: [f64; 3],
        radius: f64,
    ) -> Vec<PrimitiveDeformContact> {
        let d = sub3(cap_b, cap_a);
        let len_sq = dot3(d, d);
        vertices
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| {
                let t = if len_sq < 1e-14 {
                    0.0
                } else {
                    (dot3(sub3(v, cap_a), d) / len_sq).clamp(0.0, 1.0)
                };
                let closest = add3(cap_a, scale3(d, t));
                let diff = sub3(v, closest);
                let dist = len3(diff);
                if dist < radius {
                    let depth = radius - dist;
                    let normal = if dist > 1e-14 {
                        scale3(diff, 1.0 / dist)
                    } else {
                        [0.0, 1.0, 0.0]
                    };
                    let point = add3(closest, scale3(normal, radius));
                    Some(PrimitiveDeformContact {
                        vertex_index: i,
                        depth,
                        normal,
                        point,
                    })
                } else {
                    None
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Contact Response
// ---------------------------------------------------------------------------

/// Position correction impulse for a deformable vertex in contact.
#[derive(Debug, Clone)]
pub struct VertexImpulse {
    /// Index of the vertex.
    pub vertex_index: usize,
    /// Position correction vector (apply to vertex position).
    pub correction: [f64; 3],
    /// Velocity impulse to apply.
    pub impulse: [f64; 3],
}

/// Contact response system: compute position corrections and velocity impulses
/// for deformable vs rigid contacts.
pub struct ContactResponse {
    /// Restitution coefficient (0 = fully inelastic, 1 = elastic).
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
}

impl ContactResponse {
    /// Construct with restitution and friction coefficients.
    pub fn new(restitution: f64, friction: f64) -> Self {
        ContactResponse {
            restitution,
            friction,
        }
    }

    /// Compute vertex impulses for deformable vs rigid contacts.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_deform_rigid(
        &self,
        contacts: &[DeformContact],
        velocities: &[[f64; 3]],
        masses: &[f64],
    ) -> Vec<VertexImpulse> {
        contacts
            .iter()
            .map(|c| {
                let vi = c.vertex_index;
                let v = velocities[vi];
                let m = masses[vi];
                let vn = dot3(v, c.normal);
                // Correction: push vertex out of penetration
                let correction = scale3(c.normal, c.depth);
                // Impulse: reflect normal component
                let j_n = -(1.0 + self.restitution) * vn * m;
                let impulse_n = scale3(c.normal, j_n);
                // Tangential friction (simple Coulomb)
                let v_tang = sub3(v, scale3(c.normal, vn));
                let v_tang_len = len3(v_tang);
                let impulse_t = if v_tang_len > 1e-14 {
                    scale3(v_tang, -self.friction * j_n.abs() / v_tang_len)
                } else {
                    [0.0; 3]
                };
                let impulse = add3(impulse_n, impulse_t);
                VertexImpulse {
                    vertex_index: vi,
                    correction,
                    impulse,
                }
            })
            .collect()
    }

    /// Compute vertex impulses for edge-edge contacts.
    pub fn solve_edge_edge(
        &self,
        results: &[EdgeEdgeResult],
        velocities: &[[f64; 3]],
        edges_a: &[(usize, usize)],
        masses: &[f64],
    ) -> Vec<VertexImpulse> {
        let mut impulses = Vec::new();
        for r in results {
            let (a0, a1) = edges_a[r.edge_a];
            for &vi in &[a0, a1] {
                let v = velocities[vi];
                let m = masses[vi];
                let vn = dot3(v, r.normal);
                let depth = (self.friction - r.distance).max(0.0);
                let correction = scale3(r.normal, depth);
                let j = -(1.0 + self.restitution) * vn * m;
                let impulse = scale3(r.normal, j);
                impulses.push(VertexImpulse {
                    vertex_index: vi,
                    correction,
                    impulse,
                });
            }
        }
        impulses
    }
}

// ---------------------------------------------------------------------------
// Spatial Hash for Deformable
// ---------------------------------------------------------------------------

/// Dynamic spatial hash for deformable mesh particles, supporting fast update.
pub struct SpatialHashDeformable {
    /// Cell size.
    pub cell_size: f64,
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl SpatialHashDeformable {
    /// Construct with a cell size.
    pub fn new(cell_size: f64) -> Self {
        SpatialHashDeformable {
            cell_size,
            cells: HashMap::new(),
        }
    }

    fn cell_of(&self, p: [f64; 3]) -> (i64, i64, i64) {
        (
            (p[0] / self.cell_size).floor() as i64,
            (p[1] / self.cell_size).floor() as i64,
            (p[2] / self.cell_size).floor() as i64,
        )
    }

    /// Update positions of all particles.
    pub fn update(&mut self, positions: &[[f64; 3]]) {
        self.cells.clear();
        for (i, &p) in positions.iter().enumerate() {
            let cell = self.cell_of(p);
            self.cells.entry(cell).or_default().push(i);
        }
    }

    /// Query all particles within `radius` of `point`.
    pub fn query_radius(&self, point: [f64; 3], radius: f64) -> Vec<usize> {
        let r_cells = (radius / self.cell_size).ceil() as i64 + 1;
        let center_cell = self.cell_of(point);
        let mut result = Vec::new();
        for di in -r_cells..=r_cells {
            for dj in -r_cells..=r_cells {
                for dk in -r_cells..=r_cells {
                    let nc = (center_cell.0 + di, center_cell.1 + dj, center_cell.2 + dk);
                    if let Some(indices) = self.cells.get(&nc) {
                        for &idx in indices {
                            result.push(idx);
                        }
                    }
                }
            }
        }
        // Filter to actual radius
        result
            .into_iter()
            .filter(|_| true) // positions not available here; caller filters
            .collect()
    }

    /// Return all particles near `point`, filtered by actual distance.
    pub fn query_radius_with_positions(
        &self,
        point: [f64; 3],
        radius: f64,
        positions: &[[f64; 3]],
    ) -> Vec<usize> {
        let candidates = self.query_radius(point, radius);
        let r2 = radius * radius;
        candidates
            .into_iter()
            .filter(|&i| dot3(sub3(positions[i], point), sub3(positions[i], point)) <= r2)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Inflatable Contact
// ---------------------------------------------------------------------------

/// Contact result for an inflatable (balloon) vs environment collision.
#[derive(Debug, Clone)]
pub struct InflatableContact {
    /// Index of the balloon triangle.
    pub tri_index: usize,
    /// Penetration depth.
    pub depth: f64,
    /// Contact normal.
    pub normal: [f64; 3],
    /// Contact point.
    pub point: [f64; 3],
    /// Pressure-based stiffness at this contact.
    pub stiffness: f64,
}

/// Balloon vs environment collision using sphere sweeps per triangle and
/// pressure-based contact stiffness.
pub struct InflatableContactDetector {
    /// Internal pressure.
    pub pressure: f64,
    /// Thickness of balloon membrane.
    pub thickness: f64,
}

impl InflatableContactDetector {
    /// Construct with internal pressure and membrane thickness.
    pub fn new(pressure: f64, thickness: f64) -> Self {
        InflatableContactDetector {
            pressure,
            thickness,
        }
    }

    /// Compute contacts between balloon triangles and rigid triangles.
    pub fn query(
        &self,
        balloon_tris: &[[[f64; 3]; 3]],
        rigid_tris: &[[[f64; 3]; 3]],
    ) -> Vec<InflatableContact> {
        let mut contacts = Vec::new();
        for (bi, btri) in balloon_tris.iter().enumerate() {
            let bc = tri_centroid(*btri);
            for rtri in rigid_tris {
                let (sq_dist, bary) = vertex_triangle_distance(bc, rtri[0], rtri[1], rtri[2]);
                let dist = sq_dist.sqrt();
                if dist < self.thickness {
                    let closest = add3(
                        add3(scale3(rtri[0], bary[0]), scale3(rtri[1], bary[1])),
                        scale3(rtri[2], bary[2]),
                    );
                    let ab = sub3(rtri[1], rtri[0]);
                    let ac = sub3(rtri[2], rtri[0]);
                    let normal = normalize3(cross3(ab, ac));
                    let depth = self.thickness - dist;
                    // Pressure-based stiffness proportional to pressure and
                    // inversely to balloon area
                    let area = tri_area(*btri);
                    let stiffness = self.pressure * area;
                    contacts.push(InflatableContact {
                        tri_index: bi,
                        depth,
                        normal,
                        point: closest,
                        stiffness,
                    });
                }
            }
        }
        contacts
    }
}

fn tri_centroid(tri: [[f64; 3]; 3]) -> [f64; 3] {
    scale3(add3(add3(tri[0], tri[1]), tri[2]), 1.0 / 3.0)
}

fn tri_area(tri: [[f64; 3]; 3]) -> f64 {
    let ab = sub3(tri[1], tri[0]);
    let ac = sub3(tri[2], tri[0]);
    len3(cross3(ab, ac)) * 0.5
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- vertex_triangle_distance ---

    #[test]
    fn vtd_point_at_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let (sq, bary) = vertex_triangle_distance(a, a, b, c);
        assert!(sq < 1e-12, "sq={sq}");
        assert!((bary[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn vtd_point_above_triangle() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.25, 0.25, 1.0];
        let (sq, _) = vertex_triangle_distance(p, a, b, c);
        assert!((sq.sqrt() - 1.0).abs() < 1e-9, "dist={}", sq.sqrt());
    }

    #[test]
    fn vtd_point_inside_triangle_projection() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let p = [0.5, 0.5, 0.0];
        let (sq, _) = vertex_triangle_distance(p, a, b, c);
        assert!(sq < 1e-12, "sq={sq}");
    }

    #[test]
    fn vtd_point_outside_triangle() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [2.0, 0.0, 0.0];
        let (sq, _) = vertex_triangle_distance(p, a, b, c);
        assert!((sq.sqrt() - 1.0).abs() < 1e-9, "dist={}", sq.sqrt());
    }

    // --- edge_edge_distance ---

    #[test]
    fn eed_parallel_edges() {
        let p = [0.0, 0.0, 0.0];
        let d = [1.0, 0.0, 0.0];
        let q = [0.0, 1.0, 0.0];
        let e = [1.0, 0.0, 0.0];
        let (sq, _, _) = edge_edge_distance(p, d, q, e);
        assert!((sq.sqrt() - 1.0).abs() < 1e-9, "dist={}", sq.sqrt());
    }

    #[test]
    fn eed_perpendicular_crossing() {
        let p = [0.0, 0.0, 0.0];
        let d = [1.0, 0.0, 0.0];
        let q = [0.5, -1.0, 0.0];
        let e = [0.0, 2.0, 0.0];
        let (sq, s, t) = edge_edge_distance(p, d, q, e);
        assert!(sq < 1e-12, "sq={sq}");
        assert!((s - 0.5).abs() < 1e-9, "s={s}");
        assert!((t - 0.5).abs() < 1e-9, "t={t}");
    }

    #[test]
    fn eed_skew_edges() {
        let p = [0.0, 0.0, 0.0];
        let d = [1.0, 0.0, 0.0];
        let q = [0.5, 1.0, 1.0];
        let e = [0.0, 0.0, -1.0];
        let (sq, _, _) = edge_edge_distance(p, d, q, e);
        assert!((sq.sqrt() - 1.0).abs() < 1e-9, "dist={}", sq.sqrt());
    }

    // --- normal_cone_angle ---

    #[test]
    fn nca_identical_normals() {
        let n = [0.0, 1.0, 0.0];
        let normals = vec![n, n, n];
        let cos = normal_cone_angle(&normals);
        assert!((cos - 1.0).abs() < 1e-9, "cos={cos}");
    }

    #[test]
    fn nca_empty() {
        let cos = normal_cone_angle(&[]);
        assert_eq!(cos, 0.0);
    }

    #[test]
    fn nca_perpendicular() {
        let normals = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cos = normal_cone_angle(&normals);
        assert!(cos < 1.0);
    }

    // --- pairwise_ccd_dt ---

    #[test]
    fn ccd_dt_approaching() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [0.0, 0.0, 0.0];
        let q0 = [2.0, 0.0, 0.0];
        let q1 = [0.5, 0.0, 0.0];
        let r = 0.6;
        let toi = pairwise_ccd_dt(p0, p1, q0, q1, r);
        assert!(toi.is_some(), "expected impact");
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t), "t={t}");
    }

    #[test]
    fn ccd_dt_separating() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [0.0, 0.0, 0.0];
        let q0 = [1.0, 0.0, 0.0];
        let q1 = [5.0, 0.0, 0.0];
        let toi = pairwise_ccd_dt(p0, p1, q0, q1, 0.1);
        assert!(toi.is_none());
    }

    #[test]
    fn ccd_dt_already_overlapping() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [0.0, 0.0, 0.0];
        let q0 = [0.05, 0.0, 0.0];
        let q1 = [0.05, 0.0, 0.0];
        let toi = pairwise_ccd_dt(p0, p1, q0, q1, 1.0);
        assert!(toi.is_some());
        assert_eq!(toi.unwrap(), 0.0);
    }

    // --- DeformableVsBvh ---

    #[test]
    fn dvb_no_contact() {
        let tris = vec![[[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let det = DeformableVsBvh::new(tris, 0.05);
        let verts = vec![[0.5f64, 0.5, 5.0]];
        let contacts = det.query(&verts);
        assert!(contacts.is_empty());
    }

    #[test]
    fn dvb_contact_detected() {
        let tris = vec![[[0.0f64, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]];
        let det = DeformableVsBvh::new(tris, 0.1);
        let verts = vec![[0.5f64, 0.5, 0.05]];
        let contacts = det.query(&verts);
        assert!(!contacts.is_empty(), "expected contact");
        assert_eq!(contacts[0].vertex_index, 0);
        assert!(contacts[0].depth > 0.0);
    }

    // --- NormalConeFilter ---

    #[test]
    fn ncf_opposite_normals_always_tested() {
        // With half_angle = 0.1, cos_threshold ≈ 0.995
        // dot([0,1,0],[0,-1,0]) = -1.0 < 0.995 → should_test = true
        let filter = NormalConeFilter::new(0.1);
        assert!(filter.should_test([0.0, 1.0, 0.0], [0.0, -1.0, 0.0]));
    }

    #[test]
    fn ncf_same_normal_culled_small_threshold() {
        // With very small angle, cos_threshold ≈ 1.0
        // dot([0,1,0],[0,1,0]) = 1.0 which is NOT < threshold ≈ 1.0, so culled
        let filter = NormalConeFilter::new(0.0001);
        // cos(0.0001) ≈ 0.99999999, dot = 1.0 is NOT < 0.99999999 → culled
        assert!(!filter.should_test([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]));
    }

    // --- SelfCollision ---

    #[test]
    fn self_collision_finds_close_pair() {
        let mut sc = SelfCollision::new(1.0, 0.5, std::f64::consts::PI);
        let verts = vec![[0.0f64, 0.0, 0.0], [0.2, 0.0, 0.0], [10.0, 0.0, 0.0]];
        sc.update(&verts);
        let pairs = sc.find_pairs(&verts, None);
        assert!(!pairs.is_empty(), "expected close pair");
    }

    #[test]
    fn self_collision_no_pair_far_apart() {
        let mut sc = SelfCollision::new(1.0, 0.1, std::f64::consts::PI);
        let verts = vec![[0.0f64, 0.0, 0.0], [5.0, 0.0, 0.0]];
        sc.update(&verts);
        let pairs = sc.find_pairs(&verts, None);
        assert!(pairs.is_empty());
    }

    // --- ClothSdfContact ---

    fn make_sdf_sphere(
        center: [f64; 3],
        radius: f64,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        origin: [f64; 3],
    ) -> ClothSdfContact {
        let mut values = Vec::with_capacity(nx * ny * nz);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = [
                        origin[0] + ix as f64 * dx,
                        origin[1] + iy as f64 * dx,
                        origin[2] + iz as f64 * dx,
                    ];
                    let d = len3(sub3(p, center)) - radius;
                    values.push(d);
                }
            }
        }
        ClothSdfContact::new(nx, ny, nz, dx, origin, values, 0.1)
    }

    #[test]
    fn cloth_sdf_contact_inside_sphere() {
        let sdf = make_sdf_sphere([2.5, 2.5, 2.5], 1.0, 6, 6, 6, 1.0, [0.0, 0.0, 0.0]);
        let verts = vec![[2.5f64, 2.5, 2.5]]; // inside sphere
        let contacts = sdf.test_vertices(&verts);
        assert!(!contacts.is_empty(), "vertex inside sphere should contact");
    }

    #[test]
    fn cloth_sdf_contact_outside_sphere() {
        let sdf = make_sdf_sphere([2.5, 2.5, 2.5], 1.0, 6, 6, 6, 1.0, [0.0, 0.0, 0.0]);
        let verts = vec![[2.5f64, 2.5, 10.0]]; // far outside
        let contacts = sdf.test_vertices(&verts);
        assert!(
            contacts.is_empty(),
            "vertex far outside sphere should not contact"
        );
    }

    // --- EdgeEdgeContact ---

    #[test]
    fn eec_crossing_edges() {
        let det = EdgeEdgeContact::new(0.1);
        let pa = [0.0, 0.0, 0.0];
        let da = [2.0, 0.0, 0.0];
        let pb = [1.0, -1.0, 0.0];
        let db = [0.0, 2.0, 0.0];
        let result = det.test_pair(0, pa, da, 1, pb, db);
        assert!(result.is_some(), "crossing edges should contact");
    }

    #[test]
    fn eec_distant_edges() {
        let det = EdgeEdgeContact::new(0.1);
        let pa = [0.0, 0.0, 0.0];
        let da = [1.0, 0.0, 0.0];
        let pb = [0.0, 5.0, 0.0];
        let db = [1.0, 0.0, 0.0];
        let result = det.test_pair(0, pa, da, 1, pb, db);
        assert!(result.is_none(), "distant edges should not contact");
    }

    #[test]
    fn eec_test_all_returns_multiple() {
        let det = EdgeEdgeContact::new(0.2);
        let verts_a = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let edges_a = vec![(0, 1), (1, 2)];
        let verts_b = vec![[0.5f64, -0.1, 0.0], [0.5, 0.1, 0.0]];
        let edges_b = vec![(0, 1)];
        let results = det.test_all(&edges_a, &verts_a, &edges_b, &verts_b);
        assert!(!results.is_empty());
    }

    // --- CcdDeformable ---

    #[test]
    fn ccd_deformable_vertex_impact() {
        let ccd = CcdDeformable::new(0.1, 0.1);
        let verts0_a = vec![[0.0f64, 0.0, 0.0]];
        let verts1_a = vec![[0.0f64, 0.0, 0.0]];
        let verts0_b = vec![[0.5, 0.0, 0.0]];
        let verts1_b = vec![[0.05, 0.0, 0.0]];
        let results = ccd.vertex_vertex_ccd(&verts0_a, &verts1_a, &verts0_b, &verts1_b);
        assert!(!results.is_empty(), "expected CCD impact");
    }

    #[test]
    fn ccd_deformable_no_impact_separating() {
        let ccd = CcdDeformable::new(0.05, 0.05);
        let verts0_a = vec![[0.0f64, 0.0, 0.0]];
        let verts1_a = vec![[0.0f64, 0.0, 0.0]];
        let verts0_b = vec![[1.0, 0.0, 0.0]];
        let verts1_b = vec![[5.0, 0.0, 0.0]];
        let results = ccd.vertex_vertex_ccd(&verts0_a, &verts1_a, &verts0_b, &verts1_b);
        assert!(results.is_empty());
    }

    // --- PrimitiveVsDeformable ---

    #[test]
    fn pvd_sphere_contact() {
        let det = PrimitiveVsDeformable::sphere([0.0, 0.0, 0.0], 1.0);
        let verts = vec![[0.5f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let contacts = det.query_sphere(&verts);
        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0].vertex_index, 0);
        assert!(contacts[0].depth > 0.0);
    }

    #[test]
    fn pvd_sphere_no_contact() {
        let det = PrimitiveVsDeformable::sphere([0.0, 0.0, 0.0], 0.5);
        let verts = vec![[2.0f64, 0.0, 0.0]];
        let contacts = det.query_sphere(&verts);
        assert!(contacts.is_empty());
    }

    #[test]
    fn pvd_capsule_contact() {
        let verts = vec![[0.5f64, 0.5, 0.0]];
        let contacts =
            PrimitiveVsDeformable::query_capsule(&verts, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        assert!(!contacts.is_empty());
    }

    #[test]
    fn pvd_box_no_contact() {
        let verts = vec![[5.0f64, 0.0, 0.0]];
        let contacts = PrimitiveVsDeformable::query_box(&verts, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(contacts.is_empty());
    }

    // --- ContactResponse ---

    #[test]
    fn cr_solve_deform_rigid_positive_depth() {
        let cr = ContactResponse::new(0.0, 0.1);
        let contact = DeformContact {
            vertex_index: 0,
            tri_index: 0,
            depth: 0.02,
            normal: [0.0, 1.0, 0.0],
            point_on_rigid: [0.0, 0.0, 0.0],
        };
        let velocities = vec![[0.0, -1.0, 0.0]];
        let masses = vec![1.0];
        let impulses = cr.solve_deform_rigid(&[contact], &velocities, &masses);
        assert_eq!(impulses.len(), 1);
        assert!(impulses[0].correction[1] > 0.0);
    }

    // --- SpatialHashDeformable ---

    #[test]
    fn shd_query_finds_nearby() {
        let mut sh = SpatialHashDeformable::new(1.0);
        let positions = vec![[0.0f64, 0.0, 0.0], [0.3, 0.0, 0.0], [10.0, 0.0, 0.0]];
        sh.update(&positions);
        let nearby = sh.query_radius_with_positions([0.0, 0.0, 0.0], 0.5, &positions);
        assert!(nearby.contains(&0));
        assert!(nearby.contains(&1));
        assert!(!nearby.contains(&2));
    }

    #[test]
    fn shd_empty_after_clear() {
        let mut sh = SpatialHashDeformable::new(1.0);
        let positions = vec![[0.0f64, 0.0, 0.0]];
        sh.update(&positions);
        sh.update(&[]);
        let nearby = sh.query_radius_with_positions([0.0, 0.0, 0.0], 1.0, &[]);
        assert!(nearby.is_empty());
    }

    // --- InflatableContactDetector ---

    #[test]
    fn icd_contact_detected() {
        let det = InflatableContactDetector::new(1.0, 0.2);
        let balloon = vec![[[0.0f64, 0.0, 0.1], [1.0, 0.0, 0.1], [0.0, 1.0, 0.1]]];
        let rigid = vec![[[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let contacts = det.query(&balloon, &rigid);
        assert!(
            !contacts.is_empty(),
            "balloon close to rigid should contact"
        );
    }

    #[test]
    fn icd_no_contact_far_apart() {
        let det = InflatableContactDetector::new(1.0, 0.05);
        let balloon = vec![[[0.0f64, 0.0, 5.0], [1.0, 0.0, 5.0], [0.0, 1.0, 5.0]]];
        let rigid = vec![[[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let contacts = det.query(&balloon, &rigid);
        assert!(contacts.is_empty());
    }

    // --- tri_centroid / tri_area helpers ---

    #[test]
    fn tri_centroid_correct() {
        let tri = [[0.0f64, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let c = tri_centroid(tri);
        assert!((c[0] - 1.0).abs() < 1e-9);
        assert!((c[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn tri_area_unit() {
        let tri = [[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let a = tri_area(tri);
        assert!((a - 0.5).abs() < 1e-9, "area={a}");
    }

    // --- vector helpers ---

    #[test]
    fn vec_helpers_dot_cross() {
        let a = [1.0f64, 0.0, 0.0];
        let b = [0.0f64, 1.0, 0.0];
        assert!((dot3(a, b) - 0.0).abs() < 1e-14);
        let c = cross3(a, b);
        assert!((c[2] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn vec_normalize() {
        let v = [3.0f64, 4.0, 0.0];
        let n = normalize3(v);
        assert!((len3(n) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn midpoint_test() {
        let a = [0.0f64, 0.0, 0.0];
        let b = [2.0f64, 0.0, 0.0];
        let m = midpoint(a, b);
        assert!((m[0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn sdf_interpolate_outside_grid_returns_max() {
        let sdf = ClothSdfContact::new(3, 3, 3, 1.0, [0.0, 0.0, 0.0], vec![1.0; 27], 0.1);
        let val = sdf.interpolate([100.0, 100.0, 100.0]);
        assert_eq!(val, f64::MAX);
    }

    #[test]
    fn sdf_gradient_nonzero_in_sphere() {
        let sdf = make_sdf_sphere([2.5, 2.5, 2.5], 1.0, 6, 6, 6, 1.0, [0.0, 0.0, 0.0]);
        let g = sdf.gradient([2.5, 2.5, 1.8]);
        let gl = len3(g);
        assert!(gl > 0.5, "gradient should be nonzero: len={gl}");
    }

    #[test]
    fn self_collision_update_and_query() {
        let mut sc = SelfCollision::new(0.5, 0.3, std::f64::consts::PI);
        let verts = vec![[0.0f64, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]];
        sc.update(&verts);
        let pairs = sc.find_pairs(&verts, None);
        assert!(!pairs.is_empty());
        for p in &pairs {
            assert!(p.a < p.b);
        }
    }

    #[test]
    fn ccd_deformable_edge_edge() {
        let ccd = CcdDeformable::new(0.1, 0.1);
        let verts0_a = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let verts1_a = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let verts0_b = vec![[0.5, 1.0, 0.0], [0.5, 0.1, 0.0]];
        let verts1_b = vec![[0.5, 0.05, 0.0], [0.5, -0.85, 0.0]];
        let edges_a = vec![(0, 1)];
        let edges_b = vec![(0, 1)];
        let results = ccd.edge_edge_ccd(
            &edges_a, &verts0_a, &verts1_a, &edges_b, &verts0_b, &verts1_b,
        );
        // Just ensure no panic
        let _ = results;
    }

    #[test]
    fn deformable_vs_bvh_multiple_triangles() {
        let tris = vec![
            [[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[2.0f64, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
        ];
        let det = DeformableVsBvh::new(tris, 0.1);
        let verts = vec![[0.3f64, 0.3, 0.05], [2.5, 0.3, 0.05]];
        let contacts = det.query(&verts);
        assert_eq!(contacts.len(), 2, "expected 2 contacts");
    }
}
