// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polygon mesh construction from contour set, and conversion to NavMesh primitives.

use std::collections::HashMap;

use super::RecastConfig;
use super::contour::ContourSet;
use super::grid::{Pt, intersect, intersect_prop, left, left_on, next, prev};

/// The polygon mesh output of the Recast pipeline.
pub struct PolyMesh {
    /// Vertices in integer voxel coordinates: `[x, y, z]`.
    pub verts: Vec<[i32; 3]>,
    /// Convex polygons as vertex indices (3 to `nvp` indices each), wound
    /// counter-clockwise when viewed from above, so that fan triangles have
    /// upward (+Y) normals.
    pub polys: Vec<Vec<u16>>,
    /// Region ID per polygon.
    pub regs: Vec<u16>,
    /// Number of vertices.
    pub nverts: usize,
    /// Number of polygons.
    pub npolys: usize,
    /// Maximum vertices per polygon (`RecastConfig::max_verts_per_poly`, 3..=6).
    pub nvp: usize,
    /// World-space bounding box minimum corner.
    pub bmin: [f64; 3],
    /// Cell size (XZ voxel size).
    pub cs: f64,
    /// Cell height (Y voxel size).
    pub ch: f64,
}

/// Largest supported `max_verts_per_poly`.
pub const MAX_VERTS_PER_POLY: usize = 6;

/// Vertices with the same XZ corner whose heights differ by at most this many
/// voxels are welded into one mesh vertex.
const VERTEX_WELD_HEIGHT: i32 = 2;

/// Largest vertex count addressable by the `u16` polygon indices.
const MAX_MESH_VERTS: usize = 0xFFFE;

/// Build a polygon mesh from a contour set (port of Recast's `rcBuildPolyMesh`).
///
/// Each contour (a possibly non-convex polygon, holes already bridged in) is
/// ear-clipped into triangles, which are then greedily merged into convex
/// polygons of at most `cfg.max_verts_per_poly` vertices (longest shared edge
/// first).  Vertices are welded across contours, so neighbouring regions share
/// the vertices of their common portal edges.
pub fn build_poly_mesh(cs: &ContourSet, cfg: &RecastConfig) -> Result<PolyMesh, String> {
    let nvp = cfg.max_verts_per_poly;
    if !(3..=MAX_VERTS_PER_POLY).contains(&nvp) {
        return Err(format!(
            "max_verts_per_poly must be in 3..={MAX_VERTS_PER_POLY}, got {nvp}"
        ));
    }
    let mut all_verts: Vec<[i32; 3]> = Vec::new();
    let mut lookup: HashMap<(i32, i32), Vec<u16>> = HashMap::new();
    let mut all_polys: Vec<Vec<u16>> = Vec::new();
    let mut all_regs: Vec<u16> = Vec::new();

    for contour in &cs.contours {
        if contour.verts.len() < 3 {
            continue;
        }
        // Triangles that the ear clipper could produce; Recast keeps the
        // partial result when a contour is too badly simplified to finish.
        let (tris, _complete) = triangulate(&contour.verts);

        let mut local_to_mesh = Vec::with_capacity(contour.verts.len());
        for v in &contour.verts {
            local_to_mesh.push(add_vertex(&mut all_verts, &mut lookup, [v[0], v[1], v[2]])?);
        }

        let mut polys: Vec<Vec<u16>> = tris
            .iter()
            .filter_map(|t| {
                let (a, b, c) = (
                    local_to_mesh[t[0]],
                    local_to_mesh[t[1]],
                    local_to_mesh[t[2]],
                );
                (a != b && a != c && b != c).then(|| vec![a, b, c])
            })
            .collect();
        if nvp > 3 {
            merge_convex_polys(&mut polys, &all_verts, nvp);
        }
        for poly in polys {
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

/// Weld `v` into the vertex list, returning its index.
fn add_vertex(
    verts: &mut Vec<[i32; 3]>,
    lookup: &mut HashMap<(i32, i32), Vec<u16>>,
    v: [i32; 3],
) -> Result<u16, String> {
    let bucket = lookup.entry((v[0], v[2])).or_default();
    if let Some(&idx) = bucket
        .iter()
        .find(|&&idx| (verts[usize::from(idx)][1] - v[1]).abs() <= VERTEX_WELD_HEIGHT)
    {
        return Ok(idx);
    }
    if verts.len() >= MAX_MESH_VERTS {
        return Err(format!(
            "polygon mesh exceeds {MAX_MESH_VERTS} vertices (u16 index range)"
        ));
    }
    let idx = verts.len() as u16;
    verts.push(v);
    bucket.push(idx);
    Ok(idx)
}

#[inline]
fn cpt(v: &[i32; 4]) -> Pt {
    [i64::from(v[0]), i64::from(v[2])]
}

/// Triangle-ear test: the diagonal `(i, j)` is inside the polygon near `i`.
fn in_cone(i: usize, j: usize, idx: &[usize], verts: &[[i32; 4]], loose: bool) -> bool {
    let n = idx.len();
    let pi = cpt(&verts[idx[i]]);
    let pj = cpt(&verts[idx[j]]);
    let pi1 = cpt(&verts[idx[next(i, n)]]);
    let pin1 = cpt(&verts[idx[prev(i, n)]]);
    if left_on(pin1, pi, pi1) {
        // Convex vertex.
        if loose {
            return left_on(pi, pj, pin1) && left_on(pj, pi, pi1);
        }
        return left(pi, pj, pin1) && left(pj, pi, pi1);
    }
    // Reflex vertex.
    !(left_on(pi, pj, pi1) && left_on(pj, pi, pin1))
}

/// The diagonal `(i, j)` crosses no polygon edge (edges incident to `i` or
/// `j`, and edges touching the diagonal's endpoints, are ignored).
fn diagonalie(i: usize, j: usize, idx: &[usize], verts: &[[i32; 4]], loose: bool) -> bool {
    let n = idx.len();
    let d0 = cpt(&verts[idx[i]]);
    let d1 = cpt(&verts[idx[j]]);
    for k in 0..n {
        let k1 = next(k, n);
        if k == i || k1 == i || k == j || k1 == j {
            continue;
        }
        let p0 = cpt(&verts[idx[k]]);
        let p1 = cpt(&verts[idx[k1]]);
        if d0 == p0 || d1 == p0 || d0 == p1 || d1 == p1 {
            continue;
        }
        let hit = if loose {
            intersect_prop(d0, d1, p0, p1)
        } else {
            intersect(d0, d1, p0, p1)
        };
        if hit {
            return false;
        }
    }
    true
}

/// `(i, j)` is a proper internal diagonal of the polygon `idx`.
fn diagonal(i: usize, j: usize, idx: &[usize], verts: &[[i32; 4]], loose: bool) -> bool {
    in_cone(i, j, idx, verts, loose) && diagonalie(i, j, idx, verts, loose)
}

/// Ear-clipping triangulation of a contour polygon (port of Recast's
/// `triangulate`).  Shortest ears are clipped first, which yields well-shaped
/// triangles.  Returns contour-local triangle indices and whether the whole
/// polygon was triangulated.
fn triangulate(verts: &[[i32; 4]]) -> (Vec<[usize; 3]>, bool) {
    let mut idx: Vec<usize> = (0..verts.len()).collect();
    let mut tris: Vec<[usize; 3]> = Vec::new();
    if idx.len() < 3 {
        return (tris, false);
    }
    let len2 = |a: Pt, b: Pt| {
        let (dx, dz) = (b[0] - a[0], b[1] - a[1]);
        dx * dx + dz * dz
    };

    let n0 = idx.len();
    let mut ear: Vec<bool> = (0..n0)
        .map(|i1| {
            let i = prev(i1, n0);
            let i2 = next(i1, n0);
            diagonal(i, i2, &idx, verts, false)
        })
        .collect();

    while idx.len() > 3 {
        let n = idx.len();
        let mut best: Option<(i64, usize)> = None;
        for i in 0..n {
            let i1 = next(i, n);
            if ear[i1] {
                let l = len2(cpt(&verts[idx[i]]), cpt(&verts[idx[next(i1, n)]]));
                if best.is_none_or(|(bl, _)| l < bl) {
                    best = Some((l, i));
                }
            }
        }
        if best.is_none() {
            // Overlapping segments (e.g. a hole bridge) can hide every strict
            // ear; retry with a looser cone test.
            for i in 0..n {
                let i1 = next(i, n);
                let i2 = next(i1, n);
                if diagonal(i, i2, &idx, verts, true) {
                    let l = len2(cpt(&verts[idx[i]]), cpt(&verts[idx[i2]]));
                    if best.is_none_or(|(bl, _)| l < bl) {
                        best = Some((l, i));
                    }
                }
            }
        }
        let Some((_, i)) = best else {
            return (tris, false);
        };
        let i1 = next(i, n);
        let i2 = next(i1, n);
        tris.push([idx[i], idx[i1], idx[i2]]);

        idx.remove(i1);
        ear.remove(i1);
        let n = idx.len();
        let i1 = if i1 >= n { 0 } else { i1 };
        let i = prev(i1, n);
        ear[i] = diagonal(prev(i, n), i1, &idx, verts, false);
        ear[i1] = diagonal(i, next(i1, n), &idx, verts, false);
    }
    tris.push([idx[0], idx[1], idx[2]]);
    (tris, true)
}

/// `c` is strictly left of `a -> b` on mesh vertices (Recast convention).
#[inline]
fn uleft(a: [i32; 3], b: [i32; 3], c: [i32; 3]) -> bool {
    let (a, b, c) = (
        [i64::from(a[0]), i64::from(a[2])],
        [i64::from(b[0]), i64::from(b[2])],
        [i64::from(c[0]), i64::from(c[2])],
    );
    (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]) < 0
}

/// Merge value of polygons `pa` and `pb`: the squared length of their shared
/// edge when the union is convex and small enough, plus the shared edge
/// position in each polygon.
fn poly_merge_value(
    pa: &[u16],
    pb: &[u16],
    verts: &[[i32; 3]],
    nvp: usize,
) -> Option<(i64, usize, usize)> {
    let (na, nb) = (pa.len(), pb.len());
    if na + nb < 2 || na + nb - 2 > nvp {
        return None;
    }
    let mut shared: Option<(usize, usize)> = None;
    'outer: for i in 0..na {
        let (mut va0, mut va1) = (pa[i], pa[(i + 1) % na]);
        if va0 > va1 {
            std::mem::swap(&mut va0, &mut va1);
        }
        for j in 0..nb {
            let (mut vb0, mut vb1) = (pb[j], pb[(j + 1) % nb]);
            if vb0 > vb1 {
                std::mem::swap(&mut vb0, &mut vb1);
            }
            if va0 == vb0 && va1 == vb1 {
                shared = Some((i, j));
                break 'outer;
            }
        }
    }
    let (ea, eb) = shared?;
    let v = |k: u16| verts[usize::from(k)];

    // The merged polygon must stay convex at both ends of the shared edge.
    if !uleft(v(pa[(ea + na - 1) % na]), v(pa[ea]), v(pb[(eb + 2) % nb])) {
        return None;
    }
    if !uleft(v(pb[(eb + nb - 1) % nb]), v(pb[eb]), v(pa[(ea + 2) % na])) {
        return None;
    }
    let (a, b) = (v(pa[ea]), v(pa[(ea + 1) % na]));
    let dx = i64::from(a[0]) - i64::from(b[0]);
    let dz = i64::from(a[2]) - i64::from(b[2]);
    Some((dx * dx + dz * dz, ea, eb))
}

/// Polygon formed by joining `pa` and `pb` along edges `ea` / `eb`.
fn merge_poly_verts(pa: &[u16], pb: &[u16], ea: usize, eb: usize) -> Vec<u16> {
    let (na, nb) = (pa.len(), pb.len());
    let mut out = Vec::with_capacity(na + nb - 2);
    for i in 0..na - 1 {
        out.push(pa[(ea + 1 + i) % na]);
    }
    for i in 0..nb - 1 {
        out.push(pb[(eb + 1 + i) % nb]);
    }
    out
}

/// Greedily merge convex polygons that share an edge, longest edge first.
fn merge_convex_polys(polys: &mut Vec<Vec<u16>>, verts: &[[i32; 3]], nvp: usize) {
    loop {
        let mut best: Option<(i64, usize, usize, usize, usize)> = None;
        for j in 0..polys.len() {
            for k in j + 1..polys.len() {
                if let Some((val, ea, eb)) = poly_merge_value(&polys[j], &polys[k], verts, nvp)
                    && best.is_none_or(|b| val > b.0)
                {
                    best = Some((val, j, k, ea, eb));
                }
            }
        }
        let Some((_, pa, pb, ea, eb)) = best else {
            break;
        };
        let merged = merge_poly_verts(&polys[pa], &polys[pb], ea, eb);
        polys[pa] = merged;
        polys.swap_remove(pb);
    }
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
/// Every convex polygon is fan-triangulated; only vertices referenced by a
/// polygon are emitted.  The caller can feed these directly into
/// `NavMesh::from_triangles`.
pub fn poly_mesh_to_nav_mesh_primitives(pm: &PolyMesh) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut out_verts: Vec<[f64; 3]> = Vec::new();
    let mut out_tris: Vec<[u32; 3]> = Vec::new();
    // Mesh vertices are already welded; map each referenced one to its output slot.
    let mut remap: HashMap<u16, u32> = HashMap::new();

    let mut get_or_insert = |vi: u16, out_verts: &mut Vec<[f64; 3]>| -> Option<u32> {
        if let Some(&idx) = remap.get(&vi) {
            return Some(idx);
        }
        let v = pm.verts.get(usize::from(vi))?;
        let idx = out_verts.len() as u32;
        out_verts.push(to_world(v, &pm.bmin, pm.cs, pm.ch));
        remap.insert(vi, idx);
        Some(idx)
    };

    for poly in &pm.polys {
        if poly.len() < 3 {
            continue;
        }
        let Some(i0) = get_or_insert(poly[0], &mut out_verts) else {
            continue;
        };
        for k in 1..poly.len() - 1 {
            let (Some(i1), Some(i2)) = (
                get_or_insert(poly[k], &mut out_verts),
                get_or_insert(poly[k + 1], &mut out_verts),
            ) else {
                continue;
            };
            out_tris.push([i0, i1, i2]);
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recast::contour::Contour;

    /// Outline-wound L polygon (cells 1..13 x 1..5 plus 1..5 x 5..13).
    fn l_outline() -> Vec<[i32; 4]> {
        [(1, 1), (1, 13), (5, 13), (5, 5), (13, 5), (13, 1)]
            .iter()
            .map(|&(x, z)| [x, 0, z, 0])
            .collect()
    }

    fn in_l(x: f64, z: f64) -> bool {
        ((1.0..=13.0).contains(&x) && (1.0..=5.0).contains(&z))
            || ((1.0..=5.0).contains(&x) && (5.0..=13.0).contains(&z))
    }

    /// Twice the signed area (Recast convention) of a mesh polygon.
    fn poly_area2(poly: &[u16], verts: &[[i32; 3]]) -> i64 {
        let n = poly.len();
        (0..n)
            .map(|i| {
                let a = verts[usize::from(poly[i])];
                let b = verts[usize::from(poly[prev(i, n)])];
                i64::from(a[0]) * i64::from(b[2]) - i64::from(b[0]) * i64::from(a[2])
            })
            .sum()
    }

    fn is_convex(poly: &[u16], verts: &[[i32; 3]]) -> bool {
        let n = poly.len();
        (0..n).all(|i| {
            let a = verts[usize::from(poly[prev(i, n)])];
            let b = verts[usize::from(poly[i])];
            let c = verts[usize::from(poly[next(i, n)])];
            let (a, b, c) = (
                [i64::from(a[0]), i64::from(a[2])],
                [i64::from(b[0]), i64::from(b[2])],
                [i64::from(c[0]), i64::from(c[2])],
            );
            left_on(a, b, c)
        })
    }

    fn contour_set(contours: Vec<Contour>) -> ContourSet {
        ContourSet {
            contours,
            bmin: [0.0; 3],
            cs: 1.0,
            ch: 1.0,
        }
    }

    fn cfg(nvp: usize) -> RecastConfig {
        RecastConfig {
            max_verts_per_poly: nvp,
            ..Default::default()
        }
    }

    #[test]
    fn test_triangulate_non_convex_l() {
        let verts = l_outline();
        let (tris, complete) = triangulate(&verts);
        assert!(complete);
        assert_eq!(tris.len(), verts.len() - 2);
        let mut total = 0i64;
        for t in &tris {
            let (a, b, c) = (cpt(&verts[t[0]]), cpt(&verts[t[1]]), cpt(&verts[t[2]]));
            let a2 = -crate::recast::grid::area2(a, b, c);
            assert!(a2 > 0, "triangle {t:?} has wrong winding or is degenerate");
            total += a2;
            let cx = (a[0] + b[0] + c[0]) as f64 / 3.0;
            let cz = (a[1] + b[1] + c[1]) as f64 / 3.0;
            assert!(in_l(cx, cz), "centroid ({cx}, {cz}) lies outside the L");
        }
        // The L covers 12*4 + 4*8 = 80 cells.
        assert_eq!(total, 160);
    }

    #[test]
    fn test_triangle_only_mesh_when_nvp_is_three() {
        let cs = contour_set(vec![Contour {
            verts: l_outline(),
            raw_verts: Vec::new(),
            region: 1,
        }]);
        let pm = build_poly_mesh(&cs, &cfg(3)).expect("mesh");
        assert_eq!(pm.nvp, 3);
        assert_eq!(pm.npolys, 4);
        assert!(pm.polys.iter().all(|p| p.len() == 3));
        assert_eq!(pm.nverts, 6);
    }

    #[test]
    fn test_convex_merge_respects_nvp_and_area() {
        let cs = contour_set(vec![Contour {
            verts: l_outline(),
            raw_verts: Vec::new(),
            region: 1,
        }]);
        for nvp in 4..=MAX_VERTS_PER_POLY {
            let pm = build_poly_mesh(&cs, &cfg(nvp)).expect("mesh");
            assert_eq!(pm.npolys, pm.polys.len());
            // An L is not convex: it needs at least two convex pieces.
            assert!(pm.npolys >= 2, "nvp {nvp}: {} polys", pm.npolys);
            assert!(pm.npolys < 4, "nvp {nvp}: triangles were not merged");
            let mut total = 0i64;
            for p in &pm.polys {
                assert!(p.len() <= nvp);
                assert!(is_convex(p, &pm.verts), "polygon {p:?} is not convex");
                let a = poly_area2(p, &pm.verts);
                assert!(a > 0);
                total += a;
            }
            assert_eq!(total, 160, "nvp {nvp}");
            assert!(pm.regs.iter().all(|&r| r == 1));
        }
    }

    #[test]
    fn test_adjacent_contours_share_welded_vertices() {
        let square = |x0: i32, x1: i32| -> Vec<[i32; 4]> {
            [(x0, 0), (x0, 4), (x1, 4), (x1, 0)]
                .iter()
                .map(|&(x, z)| [x, 3, z, 0])
                .collect()
        };
        let cs = contour_set(vec![
            Contour {
                verts: square(0, 4),
                raw_verts: Vec::new(),
                region: 1,
            },
            Contour {
                verts: square(4, 8),
                raw_verts: Vec::new(),
                region: 2,
            },
        ]);
        let pm = build_poly_mesh(&cs, &cfg(6)).expect("mesh");
        assert_eq!(pm.nverts, 6, "the two shared portal corners are welded");
        assert_eq!(pm.npolys, 2);
        assert_eq!(pm.polys[0].len(), 4);
        assert_eq!(pm.polys[1].len(), 4);
        let (verts, tris) = poly_mesh_to_nav_mesh_primitives(&pm);
        assert_eq!(verts.len(), 6);
        assert_eq!(tris.len(), 4);
    }

    #[test]
    fn test_invalid_nvp_is_rejected() {
        let cs = contour_set(Vec::new());
        assert!(build_poly_mesh(&cs, &cfg(2)).is_err());
        assert!(build_poly_mesh(&cs, &cfg(7)).is_err());
        assert!(build_poly_mesh(&cs, &cfg(6)).is_ok());
    }

    #[test]
    fn test_triangulate_bridged_hole_polygon() {
        // 10x10 square with a 2x2 hole, bridged as the contour stage does.
        let outline: Vec<[i32; 4]> = [(0, 0), (0, 10), (10, 10), (10, 0)]
            .iter()
            .map(|&(x, z)| [x, 0, z, 0])
            .collect();
        let hole: Vec<[i32; 4]> = [(4, 4), (6, 4), (6, 6), (4, 6)]
            .iter()
            .map(|&(x, z)| [x, 0, z, 0])
            .collect();
        // Bridge outline vertex 0 <-> hole vertex 0.
        let mut ring = Vec::new();
        for i in 0..=outline.len() {
            ring.push(outline[i % outline.len()]);
        }
        for i in 0..=hole.len() {
            ring.push(hole[i % hole.len()]);
        }
        let cs = contour_set(vec![Contour {
            verts: ring,
            raw_verts: Vec::new(),
            region: 1,
        }]);
        let pm = build_poly_mesh(&cs, &cfg(6)).expect("mesh");
        let total: i64 = pm.polys.iter().map(|p| poly_area2(p, &pm.verts)).sum();
        assert_eq!(total, 2 * (100 - 4));
        for p in &pm.polys {
            assert!(is_convex(p, &pm.verts));
            // No polygon may cover the hole.
            let n = p.len() as f64;
            let cx = p
                .iter()
                .map(|&i| f64::from(pm.verts[usize::from(i)][0]))
                .sum::<f64>()
                / n;
            let cz = p
                .iter()
                .map(|&i| f64::from(pm.verts[usize::from(i)][2]))
                .sum::<f64>()
                / n;
            assert!(
                !((4.0..=6.0).contains(&cx) && (4.0..=6.0).contains(&cz)),
                "polygon {p:?} centroid ({cx}, {cz}) lies inside the hole"
            );
        }
    }
}
