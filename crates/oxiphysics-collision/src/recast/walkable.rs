// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Walkable-surface filtering passes on the open heightfield.

use super::RecastConfig;
use super::rasterize::Heightfield;

/// Propagate walkability upward to spans that can be stepped up to from below
/// (low-hanging obstacles: a walkable span just below a gap smaller than climb).
pub fn filter_low_hanging_obstacles(hf: &mut Heightfield, cfg: &RecastConfig) {
    let walkable_climb = (cfg.agent_max_climb / cfg.cell_height) as u16;
    for col in hf.spans.iter_mut() {
        let n = col.len();
        for i in 0..n {
            if !col[i].walkable {
                continue;
            }
            if i + 1 < n {
                let cur_smax = col[i].smax;
                let next_smin = col[i + 1].smin;
                // If gap between this span's top and the next span's bottom is within
                // climb range, the next span is reachable — mark it walkable too.
                if next_smin.saturating_sub(cur_smax) <= walkable_climb {
                    col[i + 1].walkable = true;
                }
            }
        }
    }
}

/// Remove walkable flags from spans whose neighbour has a ledge drop > agent_max_climb.
pub fn filter_ledge_spans(hf: &mut Heightfield, cfg: &RecastConfig) {
    let walkable_climb = (cfg.agent_max_climb / cfg.cell_height) as i32;
    let w = hf.width;
    let h = hf.height;

    // Collect indices to un-mark (avoid borrow issues with simultaneous read/write)
    let mut to_unmark: Vec<(usize, usize)> = Vec::new();

    for iz in 0..h {
        for ix in 0..w {
            let col_idx = hf.col_idx(ix, iz);
            for si in 0..hf.spans[col_idx].len() {
                if !hf.spans[col_idx][si].walkable {
                    continue;
                }
                let s_smax = hf.spans[col_idx][si].smax as i32;

                // Check 4 cardinal neighbours
                let mut is_ledge = false;
                for (dx, dz) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                    let nx = ix as i32 + dx;
                    let nz = iz as i32 + dz;
                    if nx < 0 || nz < 0 || nx >= w as i32 || nz >= h as i32 {
                        is_ledge = true;
                        break;
                    }
                    let ncol = hf.col_idx(nx as usize, nz as usize);
                    if hf.spans[ncol].is_empty() {
                        is_ledge = true;
                        break;
                    }
                    // Find the minimum floor of reachable neighbour spans
                    let min_floor = hf.spans[ncol]
                        .iter()
                        .map(|ns| ns.smin as i32)
                        .min()
                        .unwrap_or(0);
                    if (s_smax - min_floor).abs() > walkable_climb {
                        is_ledge = true;
                        break;
                    }
                }
                if is_ledge {
                    to_unmark.push((col_idx, si));
                }
            }
        }
    }

    for (col_idx, si) in to_unmark {
        hf.spans[col_idx][si].walkable = false;
    }
}

/// Remove walkable flags from spans that don't have enough vertical clearance
/// above them to fit an agent.
pub fn filter_walkable_low_height(hf: &mut Heightfield, cfg: &RecastConfig) {
    let walkable_height = (cfg.agent_height / cfg.cell_height).ceil() as u16;
    for col in hf.spans.iter_mut() {
        let n = col.len();
        for i in 0..n {
            if !col[i].walkable {
                continue;
            }
            let clearance = if i + 1 < n {
                col[i + 1].smin.saturating_sub(col[i].smax)
            } else {
                u16::MAX
            };
            if clearance < walkable_height {
                col[i].walkable = false;
            }
        }
    }
}
