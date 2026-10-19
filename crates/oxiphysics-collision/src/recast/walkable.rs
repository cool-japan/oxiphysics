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

/// Height used for "no ceiling" and "no neighbour span above" in the ledge
/// filter (the maximum voxel height, as in Recast's `MAX_HEIGHT`).
const MAX_HEIGHT: i32 = 0xffff;

/// Offsets `(dx, dz)` of the four cardinal neighbour columns.
const NEIGHBOUR_OFFSETS: [(isize, isize); 4] = [(-1, 0), (0, 1), (1, 0), (0, -1)];

/// Remove walkable flags from spans that lie on a ledge or on a steep slope.
///
/// Port of Recast's `rcFilterLedgeSpans`. A walkable span has its floor at
/// `bot = smax` and its ceiling at `top`, the `smin` of the next span of the
/// column (no ceiling for the highest span). For each of the four neighbour
/// columns every *traversable gap* is examined: the open space below the
/// first neighbour span (whose floor is taken as `-walkable_climb`) and the
/// space above each neighbour span (floor `nbot = smax`, ceiling `ntop` = the
/// next neighbour span's `smin`). A gap is traversable from the current span
/// when the overlap of the two open intervals,
/// `min(top, ntop) - max(bot, nbot)`, exceeds the agent height. From the
/// traversable gaps the filter tracks
///
/// * `minh`, the lowest neighbour floor relative to `bot` (the deepest drop);
/// * the lowest and highest neighbour floor reachable within
///   `walkable_climb` of `bot` (including `bot` itself).
///
/// The span is marked unwalkable when `minh < -walkable_climb` (a drop larger
/// than the agent can step down: a ledge, also towards an empty column or the
/// edge of the grid, which counts as a floor at `-walkable_climb`) or when the
/// reachable neighbour floors differ by more than `walkable_climb` (a slope
/// that is too steep).
pub fn filter_ledge_spans(hf: &mut Heightfield, cfg: &RecastConfig) {
    let walkable_climb = (cfg.agent_max_climb / cfg.cell_height) as i32;
    let walkable_height = (cfg.agent_height / cfg.cell_height).ceil() as i32;
    let w = hf.width;
    let h = hf.height;

    // Collect the spans to un-mark first: every decision must see the
    // original flags and span layout.
    let mut ledges: Vec<(usize, usize)> = Vec::new();

    for iz in 0..h {
        for ix in 0..w {
            let col_idx = hf.col_idx(ix, iz);
            let column = &hf.spans[col_idx];
            for (si, span) in column.iter().enumerate() {
                if !span.walkable {
                    continue;
                }
                let bot = i32::from(span.smax);
                let top = column
                    .get(si + 1)
                    .map_or(MAX_HEIGHT, |next| i32::from(next.smin));

                // Lowest neighbour floor relative to this span, and the range
                // of neighbour floors that can be stepped onto.
                let mut min_neighbour_height = MAX_HEIGHT;
                let mut accessible_min = bot;
                let mut accessible_max = bot;

                for (dx, dz) in NEIGHBOUR_OFFSETS {
                    let neighbour_col = ix
                        .checked_add_signed(dx)
                        .filter(|&nx| nx < w)
                        .zip(iz.checked_add_signed(dz).filter(|&nz| nz < h))
                        .map(|(nx, nz)| hf.col_idx(nx, nz));
                    let Some(neighbour_col) = neighbour_col else {
                        // Outside the grid: treated as a floor `walkable_climb` below 0.
                        min_neighbour_height = min_neighbour_height.min(-walkable_climb - bot);
                        continue;
                    };
                    let neighbours = &hf.spans[neighbour_col];

                    // The open space from minus infinity up to the first
                    // neighbour span (the whole column when it is empty).
                    let first_floor = -walkable_climb;
                    let first_ceiling = neighbours
                        .first()
                        .map_or(MAX_HEIGHT, |first| i32::from(first.smin));
                    if top.min(first_ceiling) - bot.max(first_floor) > walkable_height {
                        min_neighbour_height = min_neighbour_height.min(first_floor - bot);
                    }

                    // The open space above every neighbour span.
                    for (ni, neighbour) in neighbours.iter().enumerate() {
                        let floor = i32::from(neighbour.smax);
                        let ceiling = neighbours
                            .get(ni + 1)
                            .map_or(MAX_HEIGHT, |next| i32::from(next.smin));
                        if top.min(ceiling) - bot.max(floor) <= walkable_height {
                            continue;
                        }
                        min_neighbour_height = min_neighbour_height.min(floor - bot);
                        if (floor - bot).abs() <= walkable_climb {
                            accessible_min = accessible_min.min(floor);
                            accessible_max = accessible_max.max(floor);
                        }
                    }
                }

                let drops_too_far = min_neighbour_height < -walkable_climb;
                let too_steep = accessible_max - accessible_min > walkable_climb;
                if drops_too_far || too_steep {
                    ledges.push((col_idx, si));
                }
            }
        }
    }

    for (col_idx, si) in ledges {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Configuration with one world unit per voxel so that the agent
    /// parameters are voxel counts.
    fn voxel_config(agent_height: f64, agent_max_climb: f64) -> RecastConfig {
        RecastConfig {
            cell_size: 1.0,
            cell_height: 1.0,
            agent_height,
            agent_max_climb,
            ..RecastConfig::default()
        }
    }

    fn empty_field(width: usize, height: usize) -> Heightfield {
        Heightfield::new(
            width,
            height,
            [0.0; 3],
            [width as f64, 64.0, height as f64],
            1.0,
            1.0,
        )
        .expect("small grid")
    }

    fn walkable_at(hf: &Heightfield, ix: usize, iz: usize, span: usize) -> bool {
        hf.spans[hf.col_idx(ix, iz)][span].walkable
    }

    /// A closed box: the columns under its side walls hold one merged span
    /// from the bottom to the top face, the interior columns hold a separate
    /// bottom span and top span (the rasteriser only covers surfaces).
    fn closed_box() -> Heightfield {
        let mut hf = empty_field(12, 12);
        for iz in 2..10 {
            for ix in 2..10 {
                let on_wall = ix == 2 || ix == 9 || iz == 2 || iz == 9;
                if on_wall {
                    hf.add_span(ix, iz, 0, 11, true);
                } else {
                    hf.add_span(ix, iz, 0, 1, false);
                    hf.add_span(ix, iz, 10, 11, true);
                }
            }
        }
        hf
    }

    #[test]
    fn closed_box_top_stays_walkable_and_its_rim_is_a_ledge() {
        let mut hf = closed_box();
        filter_ledge_spans(&mut hf, &voxel_config(3.0, 4.0));
        for iz in 3..9 {
            for ix in 3..9 {
                assert!(
                    walkable_at(&hf, ix, iz, 1),
                    "top of the box at ({ix}, {iz}) must stay walkable"
                );
                assert!(!walkable_at(&hf, ix, iz, 0), "bottom was never walkable");
            }
        }
        // The rim of the top face drops 11 voxels into the empty neighbours.
        for iz in 2..10 {
            for ix in 2..10 {
                if ix == 2 || ix == 9 || iz == 2 || iz == 9 {
                    assert!(
                        !walkable_at(&hf, ix, iz, 0),
                        "rim span at ({ix}, {iz}) is a ledge"
                    );
                }
            }
        }
    }

    #[test]
    fn two_storey_slab_keeps_its_upper_floor() {
        let mut hf = empty_field(10, 10);
        for iz in 1..9 {
            for ix in 1..9 {
                hf.add_span(ix, iz, 0, 1, true);
                hf.add_span(ix, iz, 20, 21, true);
            }
        }
        filter_ledge_spans(&mut hf, &voxel_config(3.0, 4.0));
        for iz in 2..8 {
            for ix in 2..8 {
                assert!(walkable_at(&hf, ix, iz, 0), "ground floor at ({ix}, {iz})");
                assert!(walkable_at(&hf, ix, iz, 1), "upper floor at ({ix}, {iz})");
            }
        }
        for ix in 1..9 {
            assert!(!walkable_at(&hf, ix, 1, 0), "ground floor edge");
            assert!(!walkable_at(&hf, ix, 1, 1), "upper floor edge");
        }
    }

    #[test]
    fn a_cliff_marks_its_lip_but_not_its_foot() {
        // Plateau of height 20 for ix < 4, low floor of height 1 for ix >= 4.
        let mut hf = empty_field(8, 4);
        for iz in 0..4 {
            for ix in 0..8 {
                if ix < 4 {
                    hf.add_span(ix, iz, 0, 20, true);
                } else {
                    hf.add_span(ix, iz, 0, 1, true);
                }
            }
        }
        filter_ledge_spans(&mut hf, &voxel_config(3.0, 4.0));
        for iz in 1..3 {
            assert!(
                !walkable_at(&hf, 3, iz, 0),
                "lip of the plateau drops 19 voxels"
            );
            assert!(
                walkable_at(&hf, 4, iz, 0),
                "foot of the cliff only sees a wall, not a drop"
            );
            assert!(walkable_at(&hf, 2, iz, 0), "plateau interior");
            assert!(walkable_at(&hf, 5, iz, 0), "low floor interior");
        }
    }

    /// A row of columns whose floors rise by `rise` voxels per column.
    fn stairs(rise: u16) -> Heightfield {
        let mut hf = empty_field(7, 3);
        for iz in 0..3 {
            for ix in 0..7usize {
                let floor = 5 + rise * ix as u16;
                hf.add_span(ix, iz, floor - 1, floor, true);
            }
        }
        hf
    }

    #[test]
    fn steps_within_the_climb_are_walkable_steeper_ones_are_not() {
        let config = voxel_config(3.0, 4.0);
        // Neighbouring floors 2 apart: the reachable range is exactly the climb.
        let mut hf = stairs(2);
        filter_ledge_spans(&mut hf, &config);
        for ix in 1..6 {
            assert!(walkable_at(&hf, ix, 1, 0), "step of 2 at column {ix}");
        }
        // Neighbouring floors 3 apart: the reachable range 6 exceeds the climb.
        let mut hf = stairs(3);
        filter_ledge_spans(&mut hf, &config);
        for ix in 1..6 {
            assert!(!walkable_at(&hf, ix, 1, 0), "slope of 3 at column {ix}");
        }
        // Neighbouring floors 5 apart: the drop exceeds the climb.
        let mut hf = stairs(5);
        filter_ledge_spans(&mut hf, &config);
        for ix in 1..6 {
            assert!(!walkable_at(&hf, ix, 1, 0), "drop of 5 at column {ix}");
        }
    }

    #[test]
    fn a_flat_floor_loses_only_its_border() {
        let mut hf = empty_field(9, 9);
        for iz in 0..9 {
            for ix in 0..9 {
                hf.add_span(ix, iz, 4, 5, true);
            }
        }
        filter_ledge_spans(&mut hf, &voxel_config(3.0, 4.0));
        for iz in 0..9 {
            for ix in 0..9 {
                let border = ix == 0 || iz == 0 || ix == 8 || iz == 8;
                assert_eq!(
                    walkable_at(&hf, ix, iz, 0),
                    !border,
                    "span ({ix}, {iz}) of a flat floor"
                );
            }
        }
    }

    #[test]
    fn a_floor_next_to_a_deep_empty_column_is_a_ledge_next_to_a_shallow_one_is_not() {
        let config = voxel_config(3.0, 4.0);
        // Floor at height 3 beside an empty column: drop to -climb is 7 > climb.
        let mut high = empty_field(3, 3);
        for iz in 0..3 {
            for ix in 0..2 {
                high.add_span(ix, iz, 2, 3, true);
            }
        }
        filter_ledge_spans(&mut high, &config);
        assert!(
            !walkable_at(&high, 1, 1, 0),
            "the empty column is 7 voxels below the floor"
        );
        // Floor at height 0 beside an empty column: the virtual floor of the
        // empty column is exactly `climb` below, which an agent can step down.
        let mut low = empty_field(3, 3);
        for iz in 0..3 {
            for ix in 0..2 {
                low.add_span(ix, iz, 0, 0, true);
            }
        }
        filter_ledge_spans(&mut low, &config);
        assert!(
            walkable_at(&low, 1, 1, 0),
            "a drop of exactly the climb is still allowed"
        );
    }
}
