// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polygon mesh construction from contour set, and conversion to NavMesh primitives.

use super::RecastConfig;
use super::contour::ContourSet;

/// The polygon mesh output of the Recast pipeline.
pub struct PolyMesh {
    /// Vertices in integer voxel coordinates: `[x, y, z]`.
    pub verts: Vec<[i32; 3]>,
    /// Polygon vertex indices (triangles after ear-clipping contours).
    /// Each inner `Vec` holds exactly 3 indices for a triangle.
    pub polys: Vec<Vec<u16>>,
    /// Region ID per polygon.
    pub regs: Vec<u16>,
    /// Number of vertices.
    pub nverts: usize,
    /// Number of polygons.
    pub npolys: usize,
    /// Maximum vertices per polygon (config-driven, capped at 6).
    pub nvp: usize,
    /// World-space bounding box minimum corner.
    pub bmin: [f64; 3],
    /// Cell size (XZ voxel size).
    pub cs: f64,
    /// Cell height (Y voxel size).
    pub ch: f64,
}

pub fn build_poly_mesh(cs: &ContourSet, cfg: &RecastConfig) -> Result<PolyMesh, String> {
    let nvp = cfg.max_verts_per_poly.min(6);
    let mut all_verts: Vec<[i32; 3]> = Vec::new();
    let mut all_polys: Vec<Vec<u16>> = Vec::new();
    let mut all_regs: Vec<u16> = Vec::new();

    for contour in &cs.contours {
        if contour.verts.len() < 3 {
            continue;
        }

        let start_v = all_verts.len();
        for v in &contour.verts {
            all_verts.push([v[0], v[1], v[2]]);
        }

        // Ear-clip (fan) triangulate this contour
        let tris = ear_clip_triangulate(&contour.verts);
        for tri in tris {
            let poly: Vec<u16> = tri.iter().map(|&i| (start_v + i) as u16).collect();
            all_polys.push(poly);
            all_regs.push(contour.region);
        }
    }

    let nverts = all_verts.len();
    let npolys = all_polys.len();

    Ok(PolyMesh {
        verts: all_verts,
        polys: all_polys,
        regs: all_regs,
        nverts,
        npolys,
        nvp,
        bmin: cs.bmin,
        cs: cs.cs,
        ch: cs.ch,
    })
}

/// Simple fan triangulation for convex polygons (v0.1 simplification).
/// Returns triangle indices into the contour vertex slice.
fn ear_clip_triangulate(verts: &[[i32; 4]]) -> Vec<[usize; 3]> {
    let n = verts.len();
    if n < 3 {
        return vec![];
    }
    if n == 3 {
        return vec![[0, 1, 2]];
    }
    // Fan from vertex 0
    (1..n - 1).map(|i| [0, i, i + 1]).collect()
}

/// Convert a `PolyMesh` to world-space triangles (flat `[f64;9]` each).
pub fn poly_mesh_to_triangles(pm: &PolyMesh) -> Vec<[f64; 9]> {
    let mut tris = Vec::new();
    for poly in &pm.polys {
        if poly.len() < 3 {
            continue;
        }
        let v0 = to_world(&pm.verts[poly[0] as usize], &pm.bmin, pm.cs, pm.ch);
        for i in 1..poly.len() - 1 {
            let v1 = to_world(&pm.verts[poly[i] as usize], &pm.bmin, pm.cs, pm.ch);
            let v2 = to_world(&pm.verts[poly[i + 1] as usize], &pm.bmin, pm.cs, pm.ch);
            tris.push([
                v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2],
            ]);
        }
    }
    tris
}

/// Convert a `PolyMesh` to deduplicated vertex + triangle index lists suitable for
/// `NavMesh::from_triangles(vertices, tris)`.
///
/// This function returns `(vertices, triangle_indices)` in world-space coordinates.
/// The caller can feed these directly into `NavMesh::from_triangles`.
pub fn poly_mesh_to_nav_mesh_primitives(pm: &PolyMesh) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    // Convert all PolyMesh voxel-space vertices to world-space
    let world_verts: Vec<[f64; 3]> = pm
        .verts
        .iter()
        .map(|v| to_world(v, &pm.bmin, pm.cs, pm.ch))
        .collect();

    // Each polygon is already a triangle (3-element Vec) after ear-clipping
    let mut out_verts: Vec<[f64; 3]> = Vec::new();
    let mut out_tris: Vec<[u32; 3]> = Vec::new();

    // Build a deduplication map: world vertex (quantized to avoid float keys) → index
    // Use a Vec-based linear search (mesh sizes are small enough for v0.1)
    let mut vert_map: Vec<([i32; 3], u32)> = Vec::new();

    let mut get_or_insert = |vkey: [i32; 3], world: [f64; 3]| -> u32 {
        for &(k, idx) in &vert_map {
            if k == vkey {
                return idx;
            }
        }
        let idx = out_verts.len() as u32;
        out_verts.push(world);
        vert_map.push((vkey, idx));
        idx
    };

    for poly in &pm.polys {
        if poly.len() < 3 {
            continue;
        }
        // Triangles are already fan-split into exactly 3 vertices
        let v0_key = pm.verts[poly[0] as usize];
        let v1_key = pm.verts[poly[1] as usize];
        let v2_key = pm.verts[poly[2] as usize];

        let w0 = world_verts[poly[0] as usize];
        let w1 = world_verts[poly[1] as usize];
        let w2 = world_verts[poly[2] as usize];

        let i0 = get_or_insert(v0_key, w0);
        let i1 = get_or_insert(v1_key, w1);
        let i2 = get_or_insert(v2_key, w2);

        out_tris.push([i0, i1, i2]);
    }

    (out_verts, out_tris)
}

#[inline]
fn to_world(v: &[i32; 3], bmin: &[f64; 3], cs: f64, ch: f64) -> [f64; 3] {
    [
        bmin[0] + v[0] as f64 * cs,
        bmin[1] + v[1] as f64 * ch,
        bmin[2] + v[2] as f64 * cs,
    ]
}
