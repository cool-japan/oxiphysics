// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contour tracing and simplification (port of Recast's `rcBuildContours`).
//!
//! For every region, the boundary between the region and anything else
//! (walls, ledges, other regions) is walked span edge by span edge, producing a
//! raw ring of cell-corner vertices tagged with the region on the other side of
//! each edge.  The raw ring is then simplified:
//!
//! * vertices where the neighbouring region changes are always kept, so that
//!   portals between regions are identical when seen from both sides;
//! * wall edges are refined Douglas–Peucker style until every raw vertex lies
//!   within `edge_max_error` **cells** of the simplified outline;
//! * wall edges longer than `edge_max_len` (world units) are split.
//!
//! A region that wraps around an obstacle (a pillar, an enclosed hole) yields an
//! outline ring plus one hole ring per obstacle.  Holes are bridged into their
//! outline (Recast's `mergeRegionHoles`) so that every region ends up as a single
//! weakly-simple polygon that the polygon-mesh stage can triangulate.
//!
//! Contours are emitted in ascending region-id order (rings of one region in
//! row-major scan order), independent of any hashing.

use std::collections::BTreeMap;

use super::RecastConfig;
use super::compact::CompactHeightfield;
use super::grid::{Pt, cell_span_range, intersect, left, left_on, next, prev, rc_neighbour};
use super::region::RegionHeightfield;

/// A contour is a closed polygon bounding a region.
///
/// Vertices are stored as `[x, y, z, neighbour_region]`:
/// `x`/`z` are cell-corner coordinates (world position `bmin + x * cs`),
/// `y` is the floor height of the corner in voxels (world `bmin.y + y * ch`)
/// and `neighbour_region` is the region across the edge that starts at this
/// vertex (0 = wall / unwalkable).
///
/// Outlines are wound so that [`Contour::signed_area2`] is positive; any holes
/// of the region have already been bridged into `verts`.
#[derive(Debug, Clone)]
pub struct Contour {
    /// Simplified polygon vertices (holes bridged in).
    pub verts: Vec<[i32; 4]>,
    /// Raw, unsimplified boundary ring of the outline (cell-corner resolution).
    pub raw_verts: Vec<[i32; 4]>,
    /// Region this contour bounds.
    pub region: u16,
}

impl Contour {
    /// Twice the signed area of the simplified polygon in cell units
    /// (positive for outlines, negative for unmerged hole rings).
    pub fn signed_area2(&self) -> i64 {
        signed_area2(&self.verts)
    }
}

/// All contours of a region heightfield.
pub struct ContourSet {
    /// One contour per region outline, sorted by region id.
    pub contours: Vec<Contour>,
    /// World-space minimum corner of the grid.
    pub bmin: [f64; 3],
    /// Cell size (XZ).
    pub cs: f64,
    /// Cell height (Y).
    pub ch: f64,
}

/// Build simplified contours for every region of `rhf`.
///
/// `chf` must describe the same grid as `rhf.compact` (it is used to verify
/// that the inputs belong together); `cfg.edge_max_error` (cells) and
/// `cfg.edge_max_len` (world units, `0` disables splitting) drive the
/// simplification.
pub fn build_contours(
    rhf: &RegionHeightfield,
    chf: &CompactHeightfield,
    cfg: &RecastConfig,
) -> Result<ContourSet, String> {
    let field = &rhf.compact;
    if field.width != chf.width
        || field.height != chf.height
        || field.spans.len() != chf.spans.len()
    {
        return Err(format!(
            "region heightfield ({}x{}, {} spans) does not match compact heightfield \
             ({}x{}, {} spans)",
            field.width,
            field.height,
            field.spans.len(),
            chf.width,
            chf.height,
            chf.spans.len()
        ));
    }
    if !(cfg.edge_max_error.is_finite() && cfg.edge_max_error >= 0.0) {
        return Err(format!(
            "edge_max_error must be finite and >= 0, got {}",
            cfg.edge_max_error
        ));
    }
    if !(cfg.edge_max_len.is_finite() && cfg.edge_max_len >= 0.0) {
        return Err(format!(
            "edge_max_len must be finite and >= 0, got {}",
            cfg.edge_max_len
        ));
    }

    let mut set = ContourSet {
        contours: Vec::new(),
        bmin: field.bmin,
        cs: field.cs,
        ch: field.ch,
    };
    let n = field.spans.len();
    if n == 0 {
        return Ok(set);
    }

    let max_error = cfg.edge_max_error;
    let max_edge_len = if field.cs > 0.0 {
        let cells = (cfg.edge_max_len / field.cs).floor();
        if cells >= i64::MAX as f64 {
            i64::MAX
        } else {
            cells as i64
        }
    } else {
        0
    };

    // Mark the edges of every span that face something other than its region.
    let mut flags = vec![0u8; n];
    for z in 0..field.height {
        for x in 0..field.width {
            let (start, end) = cell_span_range(field, x, z);
            for (i, flag) in flags.iter_mut().enumerate().take(end).skip(start) {
                let reg = field.spans[i].reg;
                if reg == 0 {
                    continue;
                }
                let mut same = 0u8;
                for rdir in 0..4 {
                    let r = rc_neighbour(field, x, z, i, rdir)
                        .map_or(0, |(_, _, ni)| field.spans[ni].reg);
                    if r == reg {
                        same |= 1 << rdir;
                    }
                }
                *flag = same ^ 0xF;
            }
        }
    }

    // Walk every boundary loop (a span may border several loops).
    let mut contours: Vec<Contour> = Vec::new();
    for z in 0..field.height {
        for x in 0..field.width {
            let (start, end) = cell_span_range(field, x, z);
            for i in start..end {
                while flags[i] != 0 {
                    let reg = field.spans[i].reg;
                    let Some(raw) = walk_contour(field, (x, z, i), &mut flags) else {
                        continue;
                    };
                    let mut simplified = simplify_contour(&raw, max_error, max_edge_len);
                    remove_degenerate_segments(&mut simplified);
                    if (simplified.len() < 3 || signed_area2(&simplified) == 0) && max_error > 0.0 {
                        // The error bound collapsed a thin (e.g. one cell wide)
                        // region into a line; keep its exact corners instead of
                        // silently dropping walkable area.
                        simplified = simplify_contour(&raw, 0.0, max_edge_len);
                        remove_degenerate_segments(&mut simplified);
                    }
                    if simplified.len() >= 3 {
                        contours.push(Contour {
                            verts: simplified,
                            raw_verts: raw,
                            region: reg,
                        });
                    }
                }
            }
        }
    }

    // Order contours by region id; the stable sort keeps the scan order of
    // the rings within one region, so the output is fully deterministic.
    contours.sort_by_key(|c| c.region);
    merge_holes(&mut contours)?;
    set.contours = contours;
    Ok(set)
}

/// Floor height of the corner that ends edge `rdir` of span `i`: the highest
/// floor among the (up to) four spans sharing that corner.
fn corner_height(chf: &CompactHeightfield, x: usize, z: usize, i: usize, rdir: usize) -> i32 {
    let mut ch = i32::from(chf.spans[i].y);
    let dirp = (rdir + 1) & 3;
    if let Some((ax, az, ai)) = rc_neighbour(chf, x, z, i, rdir) {
        ch = ch.max(i32::from(chf.spans[ai].y));
        if let Some((_, _, ai2)) = rc_neighbour(chf, ax, az, ai, dirp) {
            ch = ch.max(i32::from(chf.spans[ai2].y));
        }
    }
    if let Some((ax, az, ai)) = rc_neighbour(chf, x, z, i, dirp) {
        ch = ch.max(i32::from(chf.spans[ai].y));
        if let Some((_, _, ai2)) = rc_neighbour(chf, ax, az, ai, rdir) {
            ch = ch.max(i32::from(chf.spans[ai2].y));
        }
    }
    ch
}

/// Walk one boundary loop starting at the first flagged edge of `start`.
///
/// Returns the raw ring `[x, y, z, neighbour_region]` (one vertex per boundary
/// edge, at the corner where the edge ends), clearing the visited edge flags.
/// Returns `None` if the loop cannot be closed (inconsistent connectivity).
fn walk_contour(
    chf: &CompactHeightfield,
    start: (usize, usize, usize),
    flags: &mut [u8],
) -> Option<Vec<[i32; 4]>> {
    let (mut x, mut z, mut i) = start;
    let start_i = i;
    let reg = chf.spans.get(i)?.reg;
    let mut dir = (0..4).find(|&d| flags[i] & (1 << d) != 0)?;
    let start_dir = dir;
    let mut points: Vec<[i32; 4]> = Vec::new();

    // Every flagged edge is consumed once and every move enters a new
    // (span, direction) state, so 8 * spans iterations always suffice.
    let max_iter = flags.len().saturating_mul(8).saturating_add(16);
    for _ in 0..max_iter {
        if flags[i] & (1 << dir) != 0 {
            let py = corner_height(chf, x, z, i, dir);
            let (mut px, mut pz) = (x as i32, z as i32);
            match dir {
                0 => pz += 1,
                1 => {
                    px += 1;
                    pz += 1;
                }
                2 => px += 1,
                _ => {}
            }
            let r = rc_neighbour(chf, x, z, i, dir).map_or(0, |(_, _, ni)| chf.spans[ni].reg);
            points.push([px, py, pz, i32::from(r)]);
            flags[i] &= !(1u8 << dir);
            dir = (dir + 1) & 3; // rotate clockwise
        } else {
            let (nx, nz, ni) = rc_neighbour(chf, x, z, i, dir)?;
            if chf.spans[ni].reg != reg {
                return None;
            }
            x = nx;
            z = nz;
            i = ni;
            dir = (dir + 3) & 3; // rotate counter-clockwise
        }
        if i == start_i && dir == start_dir {
            return Some(points);
        }
    }
    None
}

/// Squared distance from point `(x, z)` to segment `p -> q`.
fn distance_pt_seg_sq(x: i32, z: i32, px: i32, pz: i32, qx: i32, qz: i32) -> f64 {
    let pqx = f64::from(qx - px);
    let pqz = f64::from(qz - pz);
    let mut dx = f64::from(x - px);
    let mut dz = f64::from(z - pz);
    let d = pqx * pqx + pqz * pqz;
    let mut t = pqx * dx + pqz * dz;
    if d > 0.0 {
        t /= d;
    }
    t = t.clamp(0.0, 1.0);
    dx = f64::from(px) + t * pqx - f64::from(x);
    dz = f64::from(pz) + t * pqz - f64::from(z);
    dx * dx + dz * dz
}

/// Simplify a raw contour ring (port of Recast's `simplifyContour`).
///
/// Returns vertices `[x, y, z, neighbour_region]`, where the region tag is the
/// region across the edge leaving that vertex.
fn simplify_contour(points: &[[i32; 4]], max_error: f64, max_edge_len: i64) -> Vec<[i32; 4]> {
    let pn = points.len();
    if pn == 0 {
        return Vec::new();
    }
    // Working representation: [x, y, z, raw_index].
    let mut simplified: Vec<[i32; 4]> = Vec::new();

    let has_connections = points.iter().any(|p| p[3] != 0);
    if has_connections {
        // Keep every vertex where the neighbouring region changes.
        for i in 0..pn {
            let ii = (i + 1) % pn;
            if points[i][3] != points[ii][3] {
                simplified.push([points[i][0], points[i][1], points[i][2], i as i32]);
            }
        }
    }
    if simplified.is_empty() {
        // No portals: seed with the lower-left and upper-right vertices.
        let (mut ll, mut ur) = (0usize, 0usize);
        for (i, p) in points.iter().enumerate() {
            let (llp, urp) = (points[ll], points[ur]);
            if p[0] < llp[0] || (p[0] == llp[0] && p[2] < llp[2]) {
                ll = i;
            }
            if p[0] > urp[0] || (p[0] == urp[0] && p[2] > urp[2]) {
                ur = i;
            }
        }
        simplified.push([points[ll][0], points[ll][1], points[ll][2], ll as i32]);
        simplified.push([points[ur][0], points[ur][1], points[ur][2], ur as i32]);
    }

    // Insert raw vertices until every wall edge is within the error bound.
    let max_error_sq = max_error * max_error;
    let mut i = 0usize;
    while i < simplified.len() {
        let ii = (i + 1) % simplified.len();
        let (mut ax, mut az, ai) = (
            simplified[i][0],
            simplified[i][2],
            simplified[i][3] as usize,
        );
        let (mut bx, mut bz, bi) = (
            simplified[ii][0],
            simplified[ii][2],
            simplified[ii][3] as usize,
        );

        // Traverse the segment in lexicographic order so that both sides of a
        // shared edge are simplified identically.
        let (mut ci, cinc, endi) = if bx > ax || (bx == ax && bz > az) {
            ((ai + 1) % pn, 1usize, bi)
        } else {
            std::mem::swap(&mut ax, &mut bx);
            std::mem::swap(&mut az, &mut bz);
            ((bi + pn - 1) % pn, pn - 1, ai)
        };

        let mut maxd = 0.0f64;
        let mut maxi: Option<usize> = None;
        // Tessellate only wall edges.
        if points[ci][3] == 0 {
            while ci != endi {
                let d = distance_pt_seg_sq(points[ci][0], points[ci][2], ax, az, bx, bz);
                if d > maxd {
                    maxd = d;
                    maxi = Some(ci);
                }
                ci = (ci + cinc) % pn;
            }
        }

        match maxi {
            Some(mi) if maxd > max_error_sq => {
                let p = points[mi];
                simplified.insert(i + 1, [p[0], p[1], p[2], mi as i32]);
            }
            _ => i += 1,
        }
    }

    // Split wall edges that are too long.
    if max_edge_len > 0 {
        let max_len_sq = max_edge_len.saturating_mul(max_edge_len);
        let mut i = 0usize;
        while i < simplified.len() {
            let ii = (i + 1) % simplified.len();
            let (ax, az, ai) = (
                i64::from(simplified[i][0]),
                i64::from(simplified[i][2]),
                simplified[i][3] as usize,
            );
            let (bx, bz, bi) = (
                i64::from(simplified[ii][0]),
                i64::from(simplified[ii][2]),
                simplified[ii][3] as usize,
            );
            let ci = (ai + 1) % pn;
            let mut maxi: Option<usize> = None;
            if points[ci][3] == 0 {
                let dx = bx - ax;
                let dz = bz - az;
                if dx * dx + dz * dz > max_len_sq {
                    // Round based on lexicographic order so that the split is
                    // consistent regardless of traversal direction.
                    let n = if bi < ai { bi + pn - ai } else { bi - ai };
                    if n > 1 {
                        maxi = Some(if bx > ax || (bx == ax && bz > az) {
                            (ai + n / 2) % pn
                        } else {
                            (ai + n.div_ceil(2)) % pn
                        });
                    }
                }
            }
            match maxi {
                Some(mi) => {
                    let p = points[mi];
                    simplified.insert(i + 1, [p[0], p[1], p[2], mi as i32]);
                }
                None => i += 1,
            }
        }
    }

    // Replace the raw index with the region of the edge leaving the vertex.
    for v in &mut simplified {
        let ai = (v[3] as usize + 1) % pn;
        v[3] = points[ai][3];
    }
    simplified
}

/// Remove consecutive vertices that coincide on the XZ plane.
fn remove_degenerate_segments(verts: &mut Vec<[i32; 4]>) {
    let mut i = 0usize;
    while i < verts.len() && verts.len() > 1 {
        let ni = (i + 1) % verts.len();
        if verts[i][0] == verts[ni][0] && verts[i][2] == verts[ni][2] {
            verts.remove(i);
        } else {
            i += 1;
        }
    }
}

/// Twice the signed XZ area of a vertex ring (positive = outline winding).
fn signed_area2(verts: &[[i32; 4]]) -> i64 {
    let n = verts.len();
    let mut area = 0i64;
    for i in 0..n {
        let vi = verts[i];
        let vj = verts[prev(i, n)];
        area += i64::from(vi[0]) * i64::from(vj[2]) - i64::from(vj[0]) * i64::from(vi[2]);
    }
    area
}

#[inline]
fn pt(v: &[i32; 4]) -> Pt {
    [i64::from(v[0]), i64::from(v[2])]
}

/// Even-odd point-in-polygon test on the XZ plane.
fn point_in_ring(p: Pt, verts: &[[i32; 4]]) -> bool {
    let n = verts.len();
    let (px, pz) = (p[0] as f64 + 0.25, p[1] as f64 + 0.125);
    let mut inside = false;
    for i in 0..n {
        let a = verts[i];
        let b = verts[prev(i, n)];
        let (ax, az) = (f64::from(a[0]), f64::from(a[2]));
        let (bx, bz) = (f64::from(b[0]), f64::from(b[2]));
        if (az > pz) != (bz > pz) && px < (bx - ax) * (pz - az) / (bz - az) + ax {
            inside = !inside;
        }
    }
    inside
}

/// Index of the leftmost (min x, then min z) vertex of a ring.
fn leftmost_vertex(verts: &[[i32; 4]]) -> usize {
    let mut best = 0usize;
    for (i, v) in verts.iter().enumerate() {
        let b = verts[best];
        if v[0] < b[0] || (v[0] == b[0] && v[2] < b[2]) {
            best = i;
        }
    }
    best
}

/// Diagonal from outline vertex `i` to `pj` lies inside the outline's cone at `i`.
fn in_cone(i: usize, verts: &[[i32; 4]], pj: Pt) -> bool {
    let n = verts.len();
    let pi = pt(&verts[i]);
    let pi1 = pt(&verts[next(i, n)]);
    let pin1 = pt(&verts[prev(i, n)]);
    if left_on(pin1, pi, pi1) {
        return left(pi, pj, pin1) && left(pj, pi, pi1);
    }
    !(left_on(pi, pj, pi1) && left_on(pj, pi, pin1))
}

/// Segment `d0 -> d1` crosses an edge of `verts` (edges incident to vertex
/// `skip` and edges sharing an endpoint with the segment are ignored).
fn intersect_seg_contour(d0: Pt, d1: Pt, skip: Option<usize>, verts: &[[i32; 4]]) -> bool {
    let n = verts.len();
    for k in 0..n {
        let k1 = next(k, n);
        if skip.is_some_and(|s| s == k || s == k1) {
            continue;
        }
        let p0 = pt(&verts[k]);
        let p1 = pt(&verts[k1]);
        if d0 == p0 || d1 == p0 || d0 == p1 || d1 == p1 {
            continue;
        }
        if intersect(d0, d1, p0, p1) {
            return true;
        }
    }
    false
}

/// Splice ring `cb` into ring `ca` through the bridge `ca[ia] <-> cb[ib]`.
fn merge_rings(ca: &[[i32; 4]], cb: &[[i32; 4]], ia: usize, ib: usize) -> Vec<[i32; 4]> {
    let (na, nb) = (ca.len(), cb.len());
    let mut out = Vec::with_capacity(na + nb + 2);
    for i in 0..=na {
        out.push(ca[(ia + i) % na]);
    }
    for i in 0..=nb {
        out.push(cb[(ib + i) % nb]);
    }
    out
}

/// Bridge the hole rings of `holes` (sorted left to right) into `outline`.
fn merge_holes_into_outline(
    outline: &mut Vec<[i32; 4]>,
    holes: &[Vec<[i32; 4]>],
    region: u16,
) -> Result<(), String> {
    for (h, hole) in holes.iter().enumerate() {
        let mut best_vertex = leftmost_vertex(hole);
        let mut bridge: Option<usize> = None;
        for _ in 0..hole.len() {
            let corner = pt(&hole[best_vertex]);
            // Candidate outline vertices whose cone contains the hole vertex,
            // shortest first.
            let mut diags: Vec<(i64, usize)> = (0..outline.len())
                .filter(|&j| in_cone(j, outline, corner))
                .map(|j| {
                    let p = pt(&outline[j]);
                    let (dx, dz) = (p[0] - corner[0], p[1] - corner[1]);
                    (dx * dx + dz * dz, j)
                })
                .collect();
            diags.sort_unstable();
            bridge = diags.into_iter().map(|(_, j)| j).find(|&j| {
                let p = pt(&outline[j]);
                !intersect_seg_contour(p, corner, Some(j), outline)
                    && !holes[h..]
                        .iter()
                        .any(|other| intersect_seg_contour(p, corner, None, other))
            });
            if bridge.is_some() {
                break;
            }
            best_vertex = (best_vertex + 1) % hole.len();
        }
        let Some(ia) = bridge else {
            return Err(format!(
                "failed to bridge a hole into the outline of region {region}"
            ));
        };
        *outline = merge_rings(outline, hole, ia, best_vertex);
    }
    Ok(())
}

/// Bridge every hole ring into the outline of its region and drop the hole
/// contours from the set.
fn merge_holes(contours: &mut Vec<Contour>) -> Result<(), String> {
    let is_hole: Vec<bool> = contours.iter().map(|c| c.signed_area2() < 0).collect();
    if !is_hole.iter().any(|&h| h) {
        return Ok(());
    }

    // region -> (outline indices, hole indices), iterated in region order.
    let mut by_region: BTreeMap<u16, (Vec<usize>, Vec<usize>)> = BTreeMap::new();
    for (ci, c) in contours.iter().enumerate() {
        let entry = by_region.entry(c.region).or_default();
        if is_hole[ci] {
            entry.1.push(ci);
        } else {
            entry.0.push(ci);
        }
    }

    for (region, (outlines, holes)) in &by_region {
        if holes.is_empty() || outlines.is_empty() {
            continue;
        }
        // Assign each hole to the outline that contains it.
        let mut assigned: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &hi in holes {
            let probe = pt(&contours[hi].verts[leftmost_vertex(&contours[hi].verts)]);
            let owner = outlines
                .iter()
                .copied()
                .find(|&oi| point_in_ring(probe, &contours[oi].verts))
                .unwrap_or(outlines[0]);
            assigned.entry(owner).or_default().push(hi);
        }
        for (oi, mut hole_ids) in assigned {
            // Merge holes left to right.
            hole_ids.sort_by_key(|&hi| {
                let v = contours[hi].verts[leftmost_vertex(&contours[hi].verts)];
                (v[0], v[2], hi)
            });
            let hole_rings: Vec<Vec<[i32; 4]>> = hole_ids
                .iter()
                .map(|&hi| contours[hi].verts.clone())
                .collect();
            let mut outline = std::mem::take(&mut contours[oi].verts);
            merge_holes_into_outline(&mut outline, &hole_rings, *region)?;
            contours[oi].verts = outline;
        }
    }

    // Holes are now part of their outlines; holes of a region without any
    // outline cannot be triangulated and carry no walkable area.
    let mut idx = 0usize;
    contours.retain(|_| {
        let keep = !is_hole[idx];
        idx += 1;
        keep
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recast::compact::build_compact_heightfield;
    use crate::recast::rasterize::Heightfield;

    fn cfg() -> RecastConfig {
        RecastConfig {
            cell_size: 1.0,
            cell_height: 0.1,
            agent_height: 1.0,
            agent_max_climb: 0.2,
            edge_max_error: 1.3,
            edge_max_len: 0.0,
            ..Default::default()
        }
    }

    /// Region heightfield on a `w x h` grid; `reg_of(x, z)` gives the region
    /// of each cell (`None` = no walkable span).
    fn region_field(
        w: usize,
        h: usize,
        reg_of: impl Fn(usize, usize) -> Option<u16>,
    ) -> (RegionHeightfield, CompactHeightfield) {
        let mut hf = Heightfield::new(w, h, [0.0; 3], [w as f64, 1.0, h as f64], 1.0, 0.1)
            .expect("small heightfield");
        for z in 0..h {
            for x in 0..w {
                if reg_of(x, z).is_some() {
                    hf.add_span(x, z, 0, 1, true);
                }
            }
        }
        let chf = build_compact_heightfield(&hf, &cfg());
        // Start from the real region builder's output, then impose the
        // desired partition so the contour stage is tested in isolation.
        let mut rhf = crate::recast::region::build_regions(&chf, &cfg()).expect("regions");
        for z in 0..h {
            for x in 0..w {
                let (s, e) = cell_span_range(&rhf.compact, x, z);
                for i in s..e {
                    rhf.compact.spans[i].reg = reg_of(x, z).unwrap_or(0);
                }
            }
        }
        (rhf, chf)
    }

    fn corners(c: &Contour) -> std::collections::BTreeSet<(i32, i32)> {
        c.verts.iter().map(|v| (v[0], v[2])).collect()
    }

    /// Reflex vertex count of an outline-wound ring.
    fn reflex_count(verts: &[[i32; 4]]) -> usize {
        let n = verts.len();
        (0..n)
            .filter(|&i| {
                let a = pt(&verts[prev(i, n)]);
                let b = pt(&verts[i]);
                let c = pt(&verts[next(i, n)]);
                !left_on(a, b, c)
            })
            .count()
    }

    fn in_l(x: usize, z: usize) -> bool {
        // Cells 1..=12 x 1..=4 plus 1..=4 x 5..=12 form an L.
        ((1..=12).contains(&x) && (1..=4).contains(&z))
            || ((1..=4).contains(&x) && (5..=12).contains(&z))
    }

    #[test]
    fn test_l_shaped_region_gives_non_convex_contour() {
        let (rhf, chf) = region_field(14, 14, |x, z| in_l(x, z).then_some(1));
        let set = build_contours(&rhf, &chf, &cfg()).expect("contours");
        assert_eq!(set.contours.len(), 1);
        let c = &set.contours[0];
        assert_eq!(c.region, 1);
        let expected: std::collections::BTreeSet<(i32, i32)> =
            [(1, 1), (13, 1), (13, 5), (5, 5), (5, 13), (1, 13)]
                .into_iter()
                .collect();
        assert_eq!(corners(c), expected);
        assert_eq!(c.verts.len(), 6);
        assert_eq!(
            reflex_count(&c.verts),
            1,
            "an L has exactly one reflex corner"
        );
        // Area of the L: 12*4 + 4*8 = 80 cells.
        assert_eq!(c.signed_area2(), 160);
        // Raw ring has one vertex per boundary edge: perimeter = 2*(12+12) = 48.
        assert_eq!(c.raw_verts.len(), 48);
        assert!(c.verts.iter().all(|v| v[3] == 0), "all edges are walls");
    }

    #[test]
    fn test_single_cell_wide_region_is_not_degenerate() {
        let (rhf, chf) = region_field(10, 3, |x, z| (z == 1 && (1..=8).contains(&x)).then_some(4));
        let expected: std::collections::BTreeSet<(i32, i32)> =
            [(1, 1), (9, 1), (9, 2), (1, 2)].into_iter().collect();
        // Both a tight and the default (1.3 cell) error bound keep the strip:
        // with 1.3 the plain Douglas-Peucker pass would collapse it to a line.
        for err in [0.5, 1.3] {
            let c = RecastConfig {
                edge_max_error: err,
                ..cfg()
            };
            let set = build_contours(&rhf, &chf, &c).expect("contours");
            assert_eq!(set.contours.len(), 1, "edge_max_error {err}");
            assert_eq!(corners(&set.contours[0]), expected, "edge_max_error {err}");
            assert_eq!(set.contours[0].signed_area2(), 16, "edge_max_error {err}");
        }
    }

    #[test]
    fn test_region_around_pillar_has_bridged_hole() {
        // 10x10 floor with a 2x2 pillar in the middle.
        let (rhf, chf) = region_field(12, 12, |x, z| {
            let floor = (1..=10).contains(&x) && (1..=10).contains(&z);
            let pillar = (5..=6).contains(&x) && (5..=6).contains(&z);
            (floor && !pillar).then_some(1)
        });
        let set = build_contours(&rhf, &chf, &cfg()).expect("contours");
        assert_eq!(
            set.contours.len(),
            1,
            "hole must be merged into the outline"
        );
        let c = &set.contours[0];
        // outline (4) + hole (4) + two duplicated bridge endpoints.
        assert_eq!(c.verts.len(), 10);
        // Net area = 100 - 4 cells.
        assert_eq!(c.signed_area2(), 192);
        let pillar: std::collections::BTreeSet<(i32, i32)> =
            [(5, 5), (7, 5), (7, 7), (5, 7)].into_iter().collect();
        assert!(pillar.is_subset(&corners(c)));
    }

    #[test]
    fn test_adjacent_regions_share_portal_vertices() {
        let (rhf, chf) = region_field(14, 8, |x, z| {
            ((1..=12).contains(&x) && (1..=6).contains(&z)).then_some(if x <= 6 { 1 } else { 2 })
        });
        let set = build_contours(&rhf, &chf, &cfg()).expect("contours");
        assert_eq!(set.contours.len(), 2);
        let a = corners(&set.contours[0]);
        let b = corners(&set.contours[1]);
        let shared: Vec<_> = a.intersection(&b).copied().collect();
        assert_eq!(shared, vec![(7, 1), (7, 7)]);
        // Each region is a 6x6 rectangle.
        assert_eq!(set.contours[0].signed_area2(), 72);
        assert_eq!(set.contours[1].signed_area2(), 72);
        // The portal edge is tagged with the neighbouring region.
        assert!(set.contours[0].verts.iter().any(|v| v[3] == 2));
        assert!(set.contours[1].verts.iter().any(|v| v[3] == 1));
    }

    #[test]
    fn test_contours_are_deterministic_and_region_ordered() {
        let make = || {
            let (rhf, chf) = region_field(20, 6, |x, z| {
                ((1..=18).contains(&x) && (1..=4).contains(&z)).then_some(match x {
                    1..=6 => 9,
                    7..=12 => 2,
                    _ => 5,
                })
            });
            build_contours(&rhf, &chf, &cfg()).expect("contours")
        };
        let a = make();
        let b = make();
        let regs: Vec<u16> = a.contours.iter().map(|c| c.region).collect();
        assert_eq!(
            regs,
            vec![2, 5, 9],
            "region-id order, independent of hashing"
        );
        for (ca, cb) in a.contours.iter().zip(b.contours.iter()) {
            assert_eq!(ca.verts, cb.verts);
        }
    }

    #[test]
    fn test_edge_max_len_splits_long_wall_edges() {
        let (rhf, chf) = region_field(26, 6, |x, z| {
            ((1..=24).contains(&x) && (1..=4).contains(&z)).then_some(1)
        });
        let unsplit = build_contours(&rhf, &chf, &cfg()).expect("contours");
        assert_eq!(unsplit.contours[0].verts.len(), 4);
        let split_cfg = RecastConfig {
            edge_max_len: 6.0,
            ..cfg()
        };
        let split = build_contours(&rhf, &chf, &split_cfg).expect("contours");
        let v = &split.contours[0].verts;
        assert!(v.len() > 4, "long edges must be split, got {}", v.len());
        let n = v.len();
        for i in 0..n {
            let a = v[i];
            let b = v[next(i, n)];
            let (dx, dz) = (i64::from(b[0] - a[0]), i64::from(b[2] - a[2]));
            assert!(
                dx * dx + dz * dz <= 36,
                "edge {a:?}->{b:?} longer than 6 cells"
            );
        }
        assert_eq!(split.contours[0].signed_area2(), 2 * 24 * 4);
    }

    #[test]
    fn test_edge_max_error_controls_staircase_simplification() {
        // A diagonal staircase boundary: row z has cells 1..=(z+1).
        let (rhf, chf) = region_field(12, 12, |x, z| {
            ((1..=10).contains(&z) && x >= 1 && x <= z).then_some(1)
        });
        let coarse = build_contours(&rhf, &chf, &cfg()).expect("contours");
        let fine_cfg = RecastConfig {
            edge_max_error: 0.0,
            ..cfg()
        };
        let fine = build_contours(&rhf, &chf, &fine_cfg).expect("contours");
        let coarse_n = coarse.contours[0].verts.len();
        let fine_n = fine.contours[0].verts.len();
        assert!(coarse_n < fine_n, "coarse {coarse_n} vs fine {fine_n}");
        // With zero error every staircase corner is kept.
        assert_eq!(
            fine_n,
            fine.contours[0].raw_verts.len() - count_collinear(&fine.contours[0].raw_verts)
        );
        // Every raw vertex is within the error bound of the coarse polygon.
        let poly = &coarse.contours[0].verts;
        for r in &coarse.contours[0].raw_verts {
            let n = poly.len();
            let d = (0..n)
                .map(|i| {
                    let a = poly[i];
                    let b = poly[next(i, n)];
                    distance_pt_seg_sq(r[0], r[2], a[0], a[2], b[0], b[2])
                })
                .fold(f64::INFINITY, f64::min);
            assert!(
                d <= 1.3 * 1.3 + 1e-9,
                "raw vertex {r:?} deviates by {}",
                d.sqrt()
            );
        }
    }

    fn count_collinear(ring: &[[i32; 4]]) -> usize {
        let n = ring.len();
        (0..n)
            .filter(|&i| {
                let a = pt(&ring[prev(i, n)]);
                let b = pt(&ring[i]);
                let c = pt(&ring[next(i, n)]);
                crate::recast::grid::collinear(a, b, c)
            })
            .count()
    }

    #[test]
    fn test_mismatched_inputs_are_rejected() {
        let (rhf, _) = region_field(6, 6, |x, z| {
            ((1..=4).contains(&x) && (1..=4).contains(&z)).then_some(1)
        });
        let (_, other) = region_field(8, 8, |x, z| {
            ((1..=4).contains(&x) && (1..=4).contains(&z)).then_some(1)
        });
        assert!(build_contours(&rhf, &other, &cfg()).is_err());
    }
}
