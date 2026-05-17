// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contour tracing: extract closed boundary polygons for each region.

use super::RecastConfig;
use super::compact::CompactHeightfield;
use super::region::RegionHeightfield;

/// A contour is a closed polygon bounding a region.
/// Vertices stored as `[x, y, z, reg_neighbour]` in integer voxel coordinates.
#[derive(Debug, Clone)]
pub struct Contour {
    pub verts: Vec<[i32; 4]>,
    pub region: u16,
}

pub struct ContourSet {
    pub contours: Vec<Contour>,
    pub bmin: [f64; 3],
    pub cs: f64,
    pub ch: f64,
}

pub fn build_contours(
    rhf: &RegionHeightfield,
    chf: &CompactHeightfield,
    _cfg: &RecastConfig,
) -> Result<ContourSet, String> {
    let n = chf.spans.len();
    if n == 0 {
        return Ok(ContourSet {
            contours: vec![],
            bmin: chf.bmin,
            cs: chf.cs,
            ch: chf.ch,
        });
    }

    // Find all unique non-zero region IDs
    let mut region_ids: std::collections::HashSet<u16> = std::collections::HashSet::new();
    for s in &rhf.compact.spans {
        if s.reg != 0 {
            region_ids.insert(s.reg);
        }
    }

    let mut contours = Vec::new();
    for &reg in &region_ids {
        if let Some(c) = trace_contour_for_region(reg, rhf, chf) {
            contours.push(c);
        }
    }

    Ok(ContourSet {
        contours,
        bmin: chf.bmin,
        cs: chf.cs,
        ch: chf.ch,
    })
}

fn trace_contour_for_region(
    region: u16,
    rhf: &RegionHeightfield,
    chf: &CompactHeightfield,
) -> Option<Contour> {
    // Collect all spans belonging to this region with their cell coordinates
    let mut pts: Vec<[i32; 4]> = Vec::new();
    for iz in 0..chf.height {
        for ix in 0..chf.width {
            let cell_idx = iz * chf.width + ix;
            let ci = chf.cells[cell_idx].index as usize;
            let cc = chf.cells[cell_idx].count as usize;
            for i in 0..cc {
                let si = ci + i;
                if rhf.compact.spans[si].reg == region {
                    let y = rhf.compact.spans[si].y as i32;
                    pts.push([ix as i32, y, iz as i32, 0]);
                }
            }
        }
    }

    if pts.is_empty() {
        return None;
    }

    // Sort by (x, z) and deduplicate for stable hull
    pts.sort_by(|a, b| a[0].cmp(&b[0]).then(a[2].cmp(&b[2])));
    pts.dedup_by(|a, b| a[0] == b[0] && a[2] == b[2]);

    if pts.len() < 3 {
        return Some(Contour { verts: pts, region });
    }

    let hull = convex_hull_2d(&pts);
    Some(Contour {
        verts: hull,
        region,
    })
}

/// 2D signed cross product using x (index 0) and z (index 2) components.
fn cross2d(o: &[i32; 4], a: &[i32; 4], b: &[i32; 4]) -> i64 {
    (a[0] - o[0]) as i64 * (b[2] - o[2]) as i64 - (a[2] - o[2]) as i64 * (b[0] - o[0]) as i64
}

/// Andrew's monotone chain convex hull on the (x, z) plane.
fn convex_hull_2d(pts: &[[i32; 4]]) -> Vec<[i32; 4]> {
    let n = pts.len();
    if n < 3 {
        return pts.to_vec();
    }
    let mut hull: Vec<[i32; 4]> = Vec::with_capacity(2 * n);

    // Lower hull
    for p in pts {
        while hull.len() >= 2 && cross2d(&hull[hull.len() - 2], &hull[hull.len() - 1], p) <= 0 {
            hull.pop();
        }
        hull.push(*p);
    }

    // Upper hull
    let lower_len = hull.len() + 1;
    for p in pts.iter().rev() {
        while hull.len() >= lower_len
            && cross2d(&hull[hull.len() - 2], &hull[hull.len() - 1], p) <= 0
        {
            hull.pop();
        }
        hull.push(*p);
    }

    hull.pop(); // last point is duplicate of first
    hull
}
