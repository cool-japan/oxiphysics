// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Region filtering and merging (port of Recast's `mergeAndFilterRegions`).
//!
//! After watershed partitioning, the floor is usually split into many small
//! regions.  This pass
//!
//! 1. removes connected groups of regions whose total span count is below
//!    `min_region_area` (isolated specks of walkable area);
//! 2. repeatedly merges every region that is small (`<= merge_region_area`
//!    spans) or not touching any wall into its smallest neighbour, provided the
//!    two regions share exactly one contiguous border (so merging never
//!    creates a hole) and do not overlap vertically;
//! 3. compacts the surviving region ids to `1..=n` in ascending order.

use super::compact::CompactHeightfield;
use super::grid::{cell_span_range, rc_neighbour};
use super::region::RegionHeightfield;

/// Per-region bookkeeping used while merging.
#[derive(Debug, Clone)]
struct RegionInfo {
    /// Current id (changes when merged into another region; 0 = removed).
    id: u16,
    /// Number of spans owned by this region.
    span_count: usize,
    /// Neighbouring region ids in boundary-walk order (0 = wall/unwalkable).
    connections: Vec<u16>,
    /// Regions stacked above or below this one in the same columns.
    floors: Vec<u16>,
    /// The region overlaps itself vertically (e.g. a spiral ramp).
    overlap: bool,
    /// Visited flag for the small-group removal pass.
    visited: bool,
    /// Pending id remap flag for the compaction pass.
    remap: bool,
}

impl RegionInfo {
    fn new(id: u16) -> Self {
        Self {
            id,
            span_count: 0,
            connections: Vec::new(),
            floors: Vec::new(),
            overlap: false,
            visited: false,
            remap: false,
        }
    }

    fn add_unique_floor(&mut self, floor: u16) {
        if !self.floors.contains(&floor) {
            self.floors.push(floor);
        }
    }

    fn is_connected_to_border(&self) -> bool {
        self.connections.contains(&0)
    }

    fn remove_adjacent_duplicate_connections(&mut self) {
        let mut i = 0;
        while i < self.connections.len() && self.connections.len() > 1 {
            let ni = (i + 1) % self.connections.len();
            if self.connections[i] == self.connections[ni] {
                self.connections.remove(i);
            } else {
                i += 1;
            }
        }
    }

    fn replace_neighbour(&mut self, old_id: u16, new_id: u16) {
        let mut changed = false;
        for c in &mut self.connections {
            if *c == old_id {
                *c = new_id;
                changed = true;
            }
        }
        for f in &mut self.floors {
            if *f == old_id {
                *f = new_id;
            }
        }
        if changed {
            self.remove_adjacent_duplicate_connections();
        }
    }
}

/// Region `a` may absorb (or be absorbed by) region `b`.
fn can_merge_with_region(a: &RegionInfo, b: &RegionInfo) -> bool {
    // More than one shared border segment would enclose a hole.
    let shared = a.connections.iter().filter(|&&c| c == b.id).count();
    if shared > 1 {
        return false;
    }
    !a.floors.contains(&b.id)
}

/// Merge `b` into `a`, splicing their ordered neighbour lists at the shared
/// border.  Returns `false` when the two regions are not adjacent.
fn merge_regions(a: &mut RegionInfo, b: &mut RegionInfo) -> bool {
    let a_id = a.id;
    let b_id = b.id;
    let acon = a.connections.clone();
    let Some(ins_a) = acon.iter().position(|&c| c == b_id) else {
        // `b` sits entirely inside a hole of `a` (its only neighbour is `a`,
        // while `a`'s outer boundary never touches `b`).  Absorbing it fills
        // the hole and leaves `a`'s outer neighbour list unchanged.
        if b.connections.as_slice() == [a_id] {
            absorb_counts(a, b);
            return true;
        }
        return false;
    };
    let Some(ins_b) = b.connections.iter().position(|&c| c == a_id) else {
        return false;
    };

    let mut merged = Vec::with_capacity(acon.len() + b.connections.len());
    let na = acon.len();
    for i in 0..na.saturating_sub(1) {
        merged.push(acon[(ins_a + 1 + i) % na]);
    }
    let nb = b.connections.len();
    for i in 0..nb.saturating_sub(1) {
        merged.push(b.connections[(ins_b + 1 + i) % nb]);
    }
    a.connections = merged;
    a.remove_adjacent_duplicate_connections();
    absorb_counts(a, b);
    true
}

/// Move `b`'s span count and floor list into `a`.
fn absorb_counts(a: &mut RegionInfo, b: &mut RegionInfo) {
    for &f in &b.floors {
        a.add_unique_floor(f);
    }
    a.span_count += b.span_count;
    b.span_count = 0;
    b.connections.clear();
}

/// Region id of the neighbour of span `i` across Recast direction `rdir`
/// (0 when there is no walkable neighbour).
fn neighbour_region(
    chf: &CompactHeightfield,
    src_reg: &[u16],
    x: usize,
    z: usize,
    i: usize,
    rdir: usize,
) -> u16 {
    rc_neighbour(chf, x, z, i, rdir)
        .and_then(|(_, _, ni)| src_reg.get(ni).copied())
        .unwrap_or(0)
}

/// Walk the boundary of the region owning span `i` (starting at edge `dir`)
/// and record the ordered list of neighbouring region ids.
fn walk_region_contour(
    chf: &CompactHeightfield,
    src_reg: &[u16],
    start: (usize, usize, usize),
    start_dir: usize,
) -> Vec<u16> {
    let (mut x, mut z, mut i) = start;
    let own = src_reg.get(i).copied().unwrap_or(0);
    let mut dir = start_dir;
    let mut cont = Vec::new();
    let mut cur_reg = neighbour_region(chf, src_reg, x, z, i, dir);
    cont.push(cur_reg);

    // Each boundary edge is visited once and each interior step moves to a
    // new (span, direction) state, so the walk is bounded by 8 * spans.
    let max_iter = chf.spans.len().saturating_mul(8).saturating_add(16);
    for _ in 0..max_iter {
        let r = neighbour_region(chf, src_reg, x, z, i, dir);
        if r != own {
            // Solid edge: record a region change and rotate clockwise.
            if r != cur_reg {
                cur_reg = r;
                cont.push(cur_reg);
            }
            dir = (dir + 1) & 3;
        } else {
            match rc_neighbour(chf, x, z, i, dir) {
                Some((nx, nz, ni)) => {
                    x = nx;
                    z = nz;
                    i = ni;
                }
                None => break,
            }
            dir = (dir + 3) & 3;
        }
        if i == start.2 && dir == start_dir {
            break;
        }
    }

    // Remove adjacent duplicates (cyclically).
    if cont.len() > 1 {
        let mut j = 0;
        while j < cont.len() && cont.len() > 1 {
            let nj = (j + 1) % cont.len();
            if cont[j] == cont[nj] {
                cont.remove(j);
            } else {
                j += 1;
            }
        }
    }
    cont
}

/// Filter out tiny region groups and merge small regions into neighbours.
///
/// * `min_region_area` — connected groups of regions with fewer spans than
///   this are removed (their spans get region 0).
/// * `merge_region_area` — regions with at most this many spans (and regions
///   that do not touch any wall) are merged into their smallest compatible
///   neighbour.
///
/// Region ids are compacted to `1..=n`; returns `n`.
pub fn merge_and_filter_regions(
    rhf: &mut RegionHeightfield,
    min_region_area: usize,
    merge_region_area: usize,
) -> Result<u16, String> {
    let chf = &rhf.compact;
    let src_reg: Vec<u16> = chf.spans.iter().map(|s| s.reg).collect();
    let max_id = src_reg.iter().copied().max().unwrap_or(0);
    if max_id == 0 {
        return Ok(0);
    }
    let nreg = usize::from(max_id) + 1;
    let mut regions: Vec<RegionInfo> = (0..nreg).map(|i| RegionInfo::new(i as u16)).collect();

    // Gather span counts, vertical overlaps and ordered neighbour lists.
    for z in 0..chf.height {
        for x in 0..chf.width {
            let (start, end) = cell_span_range(chf, x, z);
            for i in start..end {
                let r = src_reg[i];
                if r == 0 {
                    continue;
                }
                let ri = usize::from(r);
                regions[ri].span_count += 1;

                for (j, &floor_id) in src_reg.iter().enumerate().take(end).skip(start) {
                    if i == j || floor_id == 0 {
                        continue;
                    }
                    if floor_id == r {
                        regions[ri].overlap = true;
                    }
                    regions[ri].add_unique_floor(floor_id);
                }

                if !regions[ri].connections.is_empty() {
                    continue;
                }
                let solid_dir =
                    (0..4).find(|&rdir| neighbour_region(chf, &src_reg, x, z, i, rdir) != r);
                if let Some(dir) = solid_dir {
                    regions[ri].connections = walk_region_contour(chf, &src_reg, (x, z, i), dir);
                }
            }
        }
    }

    // Ids that own no span (gaps in a sparse id range) are not regions.
    for reg in regions.iter_mut() {
        if reg.span_count == 0 {
            reg.id = 0;
        }
    }

    // Remove connected groups that are too small.
    let mut stack: Vec<usize> = Vec::new();
    let mut trace: Vec<usize> = Vec::new();
    for i in 1..nreg {
        if regions[i].id == 0 || regions[i].span_count == 0 || regions[i].visited {
            continue;
        }
        let mut span_count = 0usize;
        stack.clear();
        trace.clear();
        regions[i].visited = true;
        stack.push(i);
        while let Some(ri) = stack.pop() {
            span_count += regions[ri].span_count;
            trace.push(ri);
            let connections = regions[ri].connections.clone();
            for c in connections {
                let ci = usize::from(c);
                if c == 0 || ci >= nreg {
                    continue;
                }
                if regions[ci].visited || regions[ci].id == 0 {
                    continue;
                }
                regions[ci].visited = true;
                stack.push(ci);
            }
        }
        if span_count < min_region_area {
            for &ri in &trace {
                regions[ri].span_count = 0;
                regions[ri].id = 0;
            }
        }
    }

    // Merge small regions into their smallest compatible neighbour.
    loop {
        let mut merge_count = 0usize;
        for i in 1..nreg {
            {
                let reg = &regions[i];
                if reg.id == 0 || reg.overlap || reg.span_count == 0 {
                    continue;
                }
                if reg.span_count > merge_region_area && reg.is_connected_to_border() {
                    continue;
                }
            }
            let mut smallest = usize::MAX;
            let mut merge_id = regions[i].id;
            for &c in &regions[i].connections {
                let ci = usize::from(c);
                if c == 0 || ci >= nreg {
                    continue;
                }
                let mreg = &regions[ci];
                if mreg.id == 0 || mreg.overlap {
                    continue;
                }
                if mreg.span_count < smallest
                    && can_merge_with_region(&regions[i], mreg)
                    && can_merge_with_region(mreg, &regions[i])
                {
                    smallest = mreg.span_count;
                    merge_id = mreg.id;
                }
            }
            let old_id = regions[i].id;
            if merge_id == old_id {
                continue;
            }
            let ti = usize::from(merge_id);
            if ti == i || ti >= nreg {
                continue;
            }
            let merged = {
                let (target, source) = if ti < i {
                    let (lo, hi) = regions.split_at_mut(i);
                    (&mut lo[ti], &mut hi[0])
                } else {
                    let (lo, hi) = regions.split_at_mut(ti);
                    (&mut hi[0], &mut lo[i])
                };
                merge_regions(target, source)
            };
            if merged {
                for reg in regions.iter_mut() {
                    if reg.id == 0 {
                        continue;
                    }
                    if reg.id == old_id {
                        reg.id = merge_id;
                    }
                    reg.replace_neighbour(old_id, merge_id);
                }
                merge_count += 1;
            }
        }
        if merge_count == 0 {
            break;
        }
    }

    // Compact region ids in ascending order of their slot.
    for reg in regions.iter_mut() {
        reg.remap = reg.id != 0;
    }
    let mut next_id: u16 = 0;
    for i in 0..nreg {
        if !regions[i].remap {
            continue;
        }
        let old_id = regions[i].id;
        next_id = next_id
            .checked_add(1)
            .ok_or_else(|| "region id overflow while compacting ids".to_string())?;
        for reg in regions.iter_mut().skip(i) {
            if reg.remap && reg.id == old_id {
                reg.id = next_id;
                reg.remap = false;
            }
        }
    }

    for (span, &r) in rhf.compact.spans.iter_mut().zip(src_reg.iter()) {
        span.reg = regions.get(usize::from(r)).map_or(0, |reg| reg.id);
    }
    Ok(next_id)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recast::RecastConfig;
    use crate::recast::compact::build_compact_heightfield;
    use crate::recast::rasterize::Heightfield;

    fn cfg() -> RecastConfig {
        RecastConfig {
            cell_size: 1.0,
            cell_height: 0.1,
            agent_height: 1.0,
            agent_max_climb: 0.2,
            ..Default::default()
        }
    }

    /// Region heightfield over a `w x h` grid; `reg_of(x, z)` gives the region
    /// of each cell (`None` = no walkable span).
    fn region_field(
        w: usize,
        h: usize,
        reg_of: impl Fn(usize, usize) -> Option<u16>,
    ) -> RegionHeightfield {
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
        // desired partition so the merge pass is tested in isolation.
        let mut rhf = crate::recast::region::build_regions(&chf, &cfg()).expect("regions");
        for z in 0..h {
            for x in 0..w {
                let (s, e) = cell_span_range(&rhf.compact, x, z);
                for i in s..e {
                    rhf.compact.spans[i].reg = reg_of(x, z).unwrap_or(0);
                }
            }
        }
        rhf
    }

    fn region_ids(rhf: &RegionHeightfield) -> std::collections::BTreeSet<u16> {
        rhf.compact.spans.iter().map(|s| s.reg).collect()
    }

    #[test]
    fn test_small_neighbour_regions_are_merged() {
        // 12 x 6 floor split into a left (4 wide) and right (8 wide) region.
        let mut rhf = region_field(14, 8, |x, z| {
            if (1..=12).contains(&x) && (1..=6).contains(&z) {
                Some(if x <= 4 { 1 } else { 2 })
            } else {
                None
            }
        });
        let n = merge_and_filter_regions(&mut rhf, 1, 400).expect("merge");
        assert_eq!(n, 1);
        assert_eq!(region_ids(&rhf), [1u16].into_iter().collect());
    }

    #[test]
    fn test_large_wall_regions_are_kept_apart() {
        let mut rhf = region_field(14, 8, |x, z| {
            if (1..=12).contains(&x) && (1..=6).contains(&z) {
                Some(if x <= 6 { 1 } else { 2 })
            } else {
                None
            }
        });
        // Each region has 36 spans; merge threshold 10 keeps them separate.
        let n = merge_and_filter_regions(&mut rhf, 1, 10).expect("merge");
        assert_eq!(n, 2);
        assert_eq!(region_ids(&rhf), [1u16, 2].into_iter().collect());
    }

    #[test]
    fn test_enclosed_region_is_merged_even_if_large() {
        // Region 2 is fully enclosed by region 1, so it never touches a wall.
        let mut rhf = region_field(12, 12, |x, z| {
            if (1..=10).contains(&x) && (1..=10).contains(&z) {
                Some(if (3..=8).contains(&x) && (3..=8).contains(&z) {
                    2
                } else {
                    1
                })
            } else {
                None
            }
        });
        let n = merge_and_filter_regions(&mut rhf, 1, 0).expect("merge");
        assert_eq!(n, 1);
        assert_eq!(region_ids(&rhf), [1u16].into_iter().collect());
    }

    #[test]
    fn test_tiny_isolated_region_is_removed() {
        // A 2x2 speck far away from a 6x6 floor.
        let mut rhf = region_field(16, 8, |x, z| {
            if (1..=6).contains(&x) && (1..=6).contains(&z) {
                Some(3)
            } else if (12..=13).contains(&x) && (2..=3).contains(&z) {
                Some(7)
            } else {
                None
            }
        });
        let n = merge_and_filter_regions(&mut rhf, 8, 0).expect("merge");
        assert_eq!(n, 1);
        assert_eq!(region_ids(&rhf), [0u16, 1].into_iter().collect());
        let speck = cell_span_range(&rhf.compact, 12, 2).0;
        assert_eq!(rhf.compact.spans[speck].reg, 0);
        let floor = cell_span_range(&rhf.compact, 3, 3).0;
        assert_eq!(rhf.compact.spans[floor].reg, 1);
    }

    #[test]
    fn test_region_walk_records_ordered_neighbours() {
        // Region 1 in the middle with region 2 to the west and 3 to the east.
        let rhf = region_field(11, 5, |x, z| {
            if (1..=9).contains(&x) && (1..=3).contains(&z) {
                Some(if x <= 3 {
                    2
                } else if x >= 7 {
                    3
                } else {
                    1
                })
            } else {
                None
            }
        });
        let chf = &rhf.compact;
        let src: Vec<u16> = chf.spans.iter().map(|s| s.reg).collect();
        let start = cell_span_range(chf, 4, 1).0;
        let cont = walk_region_contour(chf, &src, (4, 1, start), 0);
        // Walls (0) separate the west and east neighbours on both sides.
        assert_eq!(cont.len(), 4);
        assert_eq!(cont.iter().filter(|&&c| c == 0).count(), 2);
        assert!(cont.contains(&2));
        assert!(cont.contains(&3));
    }
}
