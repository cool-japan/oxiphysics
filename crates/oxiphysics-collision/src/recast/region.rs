// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Watershed region partitioning of the compact heightfield.

use super::RecastConfig;
use super::compact::{CompactHeightfield, DIRS};

/// Compact heightfield with region IDs assigned to each span.
pub struct RegionHeightfield {
    pub compact: CompactHeightfield,
}

pub fn build_regions(
    chf: &CompactHeightfield,
    cfg: &RecastConfig,
) -> Result<RegionHeightfield, String> {
    let n = chf.spans.len();
    let mut regions = vec![0u16; n];
    let mut reg_id: u16 = 1;

    let min_region_area = cfg.region_min_size * cfg.region_min_size;

    // No walkable area at all — return empty regions
    if chf.max_dist == 0 || n == 0 {
        let compact_copy = clone_compact_hf(chf);
        return Ok(RegionHeightfield {
            compact: compact_copy,
        });
    }

    // Watershed: iterate from max distance to 0, flood-filling from seeds
    let level_step = 2u16;
    let max_level = (chf.max_dist / level_step) * level_step;
    let mut level = max_level;

    loop {
        // Find all unvisited spans at the current distance band and flood-fill from them
        for seed_si in 0..n {
            if regions[seed_si] != 0 {
                continue;
            }
            if chf.dist[seed_si] == u16::MAX {
                continue;
            }
            if chf.dist[seed_si] < level || chf.dist[seed_si] >= level + level_step {
                continue;
            }

            // Flood fill from this seed
            let mut stack = vec![seed_si];
            regions[seed_si] = reg_id;
            let mut count = 0usize;

            while let Some(cur) = stack.pop() {
                count += 1;
                let cy = chf.spans[cur].y as i32;
                let (ix, iz) = span_to_cell(cur, chf);

                for (dir, &(dx, dz)) in DIRS.iter().enumerate() {
                    let nc_raw = (chf.spans[cur].con >> (dir * 8)) & 0xFF;
                    if nc_raw == 0xFF {
                        continue;
                    }
                    let nx = ix as i32 + dx;
                    let nz = iz as i32 + dz;
                    if nx < 0 || nz < 0 || nx >= chf.width as i32 || nz >= chf.height as i32 {
                        continue;
                    }
                    let ncell_idx = nz as usize * chf.width + nx as usize;
                    let nci = chf.cells[ncell_idx].index as usize;
                    let nsi = nci + nc_raw as usize;
                    if regions[nsi] != 0 {
                        continue;
                    }
                    let ny = chf.spans[nsi].y as i32;
                    if (ny - cy).abs() <= chf.walkable_climb as i32 {
                        regions[nsi] = reg_id;
                        stack.push(nsi);
                    }
                }
            }

            // If region too small, unmark it so it can be absorbed later
            if count < min_region_area {
                for r in regions.iter_mut().take(n) {
                    if *r == reg_id {
                        *r = 0;
                    }
                }
            } else {
                reg_id += 1;
            }
        }

        if level == 0 {
            break;
        }
        level = level.saturating_sub(level_step);
    }

    // Apply region IDs to compact spans
    let mut compact_copy = clone_compact_hf(chf);
    for (si, r) in regions.iter().enumerate() {
        compact_copy.spans[si].reg = *r;
    }

    Ok(RegionHeightfield {
        compact: compact_copy,
    })
}

/// Reverse-lookup: given a span index, return its (ix, iz) cell coordinates.
fn span_to_cell(si: usize, chf: &CompactHeightfield) -> (usize, usize) {
    for iz in 0..chf.height {
        for ix in 0..chf.width {
            let cell_idx = iz * chf.width + ix;
            let ci = chf.cells[cell_idx].index as usize;
            let cc = chf.cells[cell_idx].count as usize;
            if si >= ci && si < ci + cc {
                return (ix, iz);
            }
        }
    }
    (0, 0)
}

fn clone_compact_hf(chf: &CompactHeightfield) -> CompactHeightfield {
    CompactHeightfield {
        width: chf.width,
        height: chf.height,
        span_count: chf.span_count,
        walkable_height: chf.walkable_height,
        walkable_climb: chf.walkable_climb,
        bmin: chf.bmin,
        bmax: chf.bmax,
        cs: chf.cs,
        ch: chf.ch,
        cells: chf.cells.clone(),
        spans: chf.spans.clone(),
        dist: chf.dist.clone(),
        max_dist: chf.max_dist,
    }
}
