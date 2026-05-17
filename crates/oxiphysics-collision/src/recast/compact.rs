// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compact heightfield: dense array of walkable spans with precomputed connectivity.

use super::RecastConfig;
use super::rasterize::Heightfield;

/// Shared direction table: 4 cardinal neighbours (dx, dz).
/// Dir 0 = +X, Dir 1 = +Z, Dir 2 = -X, Dir 3 = -Z.
pub(crate) const DIRS: [(i32, i32); 4] = [(1, 0), (0, 1), (-1, 0), (0, -1)];

#[derive(Debug, Clone, Default)]
pub struct CompactCell {
    /// Index into `CompactHeightfield::spans` for the first span in this cell.
    pub index: u32,
    /// Number of walkable spans in this cell.
    pub count: u16,
}

#[derive(Debug, Clone, Default)]
pub struct CompactSpan {
    /// Floor height (top of the walkable span) in voxel units.
    pub y: u16,
    /// Region ID (0 = unassigned).
    pub reg: u16,
    /// Connectivity info: 4 directions × 8 bits each.
    /// For direction `d`: bits `d*8 .. d*8+7`.
    /// Value 0xFF in those bits = no connection (border or no matching span).
    pub con: u32,
    /// Clearance above this span (in voxels), capped at walkable_height.
    pub h: u16,
}

pub struct CompactHeightfield {
    pub width: usize,
    pub height: usize, // Z dimension
    pub span_count: usize,
    pub walkable_height: u16,
    pub walkable_climb: u16,
    pub bmin: [f64; 3],
    pub bmax: [f64; 3],
    pub cs: f64,
    pub ch: f64,
    pub cells: Vec<CompactCell>, // width * height cells
    pub spans: Vec<CompactSpan>, // all walkable spans
    pub dist: Vec<u16>,          // distance-to-border field per span
    pub max_dist: u16,
}

pub fn build_compact_heightfield(hf: &Heightfield, cfg: &RecastConfig) -> CompactHeightfield {
    let walkable_height = (cfg.agent_height / cfg.cell_height).ceil() as u16;
    let walkable_climb = (cfg.agent_max_climb / cfg.cell_height) as u16;
    let w = hf.width;
    let h = hf.height;
    let n = w * h;

    let mut cells = vec![CompactCell::default(); n];
    let mut spans: Vec<CompactSpan> = Vec::new();

    // Build compact spans from walkable open-heightfield spans
    for iz in 0..h {
        for ix in 0..w {
            let col = &hf.spans[hf.col_idx(ix, iz)];
            let cell_idx = iz * w + ix;
            cells[cell_idx].index = spans.len() as u32;
            let start = spans.len();
            for (si, s) in col.iter().enumerate() {
                if !s.walkable {
                    continue;
                }
                let y = s.smax;
                let clearance = if si + 1 < col.len() {
                    col[si + 1].smin.saturating_sub(s.smax)
                } else {
                    u16::MAX
                };
                spans.push(CompactSpan {
                    y,
                    reg: 0,
                    con: 0xFFFFFFFF, // all directions = no connection initially
                    h: clearance.min(walkable_height),
                });
            }
            cells[cell_idx].count = (spans.len() - start) as u16;
        }
    }

    let span_count = spans.len();

    // Build connectivity using the shared DIRS table
    for iz in 0..h {
        for ix in 0..w {
            let cell_idx = iz * w + ix;
            let ci = cells[cell_idx].index as usize;
            let cc = cells[cell_idx].count as usize;
            for i in 0..cc {
                let si = ci + i;
                let s_y = spans[si].y as i32;
                let mut con = 0xFFFFFFFFu32;
                for (dir, &(dx, dz)) in DIRS.iter().enumerate() {
                    let nx = ix as i32 + dx;
                    let nz = iz as i32 + dz;
                    if nx < 0 || nz < 0 || nx >= w as i32 || nz >= h as i32 {
                        // Border: keep 0xFF in this direction slot (already set)
                        continue;
                    }
                    let ncell_idx = nz as usize * w + nx as usize;
                    let nci = cells[ncell_idx].index as usize;
                    let ncc = cells[ncell_idx].count as usize;
                    // Find the best neighbour span (closest floor height within climb)
                    let mut best = 0xFFu32;
                    for ni in 0..ncc {
                        let n_y = spans[nci + ni].y as i32;
                        if (n_y - s_y).abs() <= walkable_climb as i32 {
                            best = ni as u32;
                            break;
                        }
                    }
                    // Clear the 8-bit slot for this direction, then set new value
                    con &= !(0xFFu32 << (dir * 8));
                    con |= (best & 0xFF) << (dir * 8);
                }
                spans[si].con = con;
            }
        }
    }

    // Distance transform (two-pass, Chebyshev approximation)
    let mut dist = vec![u16::MAX; span_count];

    // Mark border spans: any span with at least one direction having no connection (0xFF)
    for si in 0..span_count {
        let con = spans[si].con;
        let is_border = (0..4).any(|dir| ((con >> (dir * 8)) & 0xFF) == 0xFF);
        if is_border {
            dist[si] = 0;
        }
    }

    // Forward pass (+X, +Z directions)
    for iz in 0..h {
        for ix in 0..w {
            let cell_idx = iz * w + ix;
            let ci = cells[cell_idx].index as usize;
            let cc = cells[cell_idx].count as usize;
            for i in 0..cc {
                let si = ci + i;
                for (dir, &(dx, dz)) in DIRS.iter().enumerate() {
                    let nc_raw = (spans[si].con >> (dir * 8)) & 0xFF;
                    if nc_raw == 0xFF {
                        continue;
                    }
                    let nx = ix as i32 + dx;
                    let nz = iz as i32 + dz;
                    if nx < 0 || nz < 0 || nx >= w as i32 || nz >= h as i32 {
                        continue;
                    }
                    let ncell_idx = nz as usize * w + nx as usize;
                    let nci = cells[ncell_idx].index as usize;
                    let nsi = nci + nc_raw as usize;
                    if dist[nsi] != u16::MAX {
                        let new_dist = dist[nsi].saturating_add(2);
                        if new_dist < dist[si] {
                            dist[si] = new_dist;
                        }
                    }
                }
            }
        }
    }

    // Backward pass (-X, -Z directions)
    for iz in (0..h).rev() {
        for ix in (0..w).rev() {
            let cell_idx = iz * w + ix;
            let ci = cells[cell_idx].index as usize;
            let cc = cells[cell_idx].count as usize;
            for i in 0..cc {
                let si = ci + i;
                for (dir, &(dx, dz)) in DIRS.iter().enumerate() {
                    let nc_raw = (spans[si].con >> (dir * 8)) & 0xFF;
                    if nc_raw == 0xFF {
                        continue;
                    }
                    let nx = ix as i32 + dx;
                    let nz = iz as i32 + dz;
                    if nx < 0 || nz < 0 || nx >= w as i32 || nz >= h as i32 {
                        continue;
                    }
                    let ncell_idx = nz as usize * w + nx as usize;
                    let nci = cells[ncell_idx].index as usize;
                    let nsi = nci + nc_raw as usize;
                    if dist[nsi] != u16::MAX {
                        let new_dist = dist[nsi].saturating_add(2);
                        if new_dist < dist[si] {
                            dist[si] = new_dist;
                        }
                    }
                }
            }
        }
    }

    let max_dist = dist
        .iter()
        .cloned()
        .filter(|&d| d != u16::MAX)
        .max()
        .unwrap_or(0);

    CompactHeightfield {
        width: w,
        height: h,
        span_count,
        walkable_height,
        walkable_climb,
        bmin: hf.bmin,
        bmax: hf.bmax,
        cs: hf.cs,
        ch: hf.ch,
        cells,
        spans,
        dist,
        max_dist,
    }
}
