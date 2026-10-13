// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Watershed partitioning of a compact heightfield.
//!
//! This is the region-growing core of Recast's `rcBuildRegions`.  The distance
//! field of the compact heightfield counts two units per cell, so the spans are
//! processed from the highest level down to level 0 in steps of two:
//!
//! 1. At every level the regions that already exist first grow into the still
//!    unassigned spans whose distance is at least the level
//!    ([`Partition::expand`]).  A region therefore claims the cells next to it
//!    before a new seed may appear there, and the growth is bounded to a few
//!    rounds per level so a basin never floods past a narrow neck.
//! 2. The spans that are still unassigned become seeds of new regions
//!    ([`Partition::flood_region`]).  A flood spreads over spans with
//!    `dist >= level - 2` only and never touches a span next to another region
//!    (diagonal contacts included), which keeps basins apart.
//! 3. A final unbounded expansion at level 0 hands the leftover spans to the
//!    neighbouring regions.
//!
//! The result is the classic watershed: one region per basin of the distance
//! field, separated along the narrow necks between basins, so doorways,
//! corridors and pillars split a room into several regions.  Merging small
//! regions and dropping specks is left to
//! [`merge_and_filter_regions`](crate::recast::region_merge::merge_and_filter_regions).
//!
//! Every stage is linear in the number of spans.  The neighbours are looked up
//! in a table that is built once (bounds-checked through
//! [`rc_neighbour`](crate::recast::grid::rc_neighbour)), the spans are bucketed
//! by level once, and the expansion revisits only the spans that link to a span
//! assigned in the previous round.

use crate::recast::compact::CompactHeightfield;
use crate::recast::grid::{cell_span_range, rc_neighbour};

/// Marker for "no span" in the neighbour table.
const NO_SPAN: u32 = u32::MAX;

/// Distance-field units between two watershed levels (the field counts 2 per cell).
const LEVEL_STEP: u32 = 2;

/// Expansion rounds per level before new seeds may appear (Recast's `expandIters`).
const EXPAND_ROUNDS: u32 = 8;

/// Generation-stamped set membership: `begin` empties the set in O(1).
struct Marks {
    stamps: Vec<u32>,
    generation: u32,
}

impl Marks {
    fn new(len: usize) -> Self {
        Self {
            stamps: vec![0; len],
            generation: 0,
        }
    }

    /// Start a new, empty generation.
    fn begin(&mut self) {
        match self.generation.checked_add(1) {
            Some(next) => self.generation = next,
            None => {
                self.stamps.fill(0);
                self.generation = 1;
            }
        }
    }

    fn mark(&mut self, index: usize) {
        self.stamps[index] = self.generation;
    }

    fn is_marked(&self, index: usize) -> bool {
        self.stamps[index] == self.generation
    }
}

/// Span connectivity of a compact heightfield in Recast direction order
/// (`0 = -X, 1 = +Z, 2 = +X, 3 = -Z`).
struct SpanGraph {
    /// `neighbours[i][rdir]`: the span reached from span `i` across `rdir`.
    neighbours: Vec<[u32; 4]>,
    /// `(target, source)` pairs of links `source -> target` that have no link
    /// back from `target`, sorted.  Connections are chosen per span (first span
    /// within the climb height), so with several floors per cell a link need
    /// not be mutual; the expansion needs the reverse direction of such links.
    one_way_in: Vec<(u32, u32)>,
}

impl SpanGraph {
    /// Build the neighbour table; fails when the span count does not fit the
    /// 32-bit span indices.
    fn build(chf: &CompactHeightfield) -> Result<Self, String> {
        let count = chf.spans.len();
        if count >= NO_SPAN as usize {
            return Err(format!(
                "{count} spans exceed the supported maximum of {}",
                NO_SPAN - 1
            ));
        }
        let mut neighbours = vec![[NO_SPAN; 4]; count];
        for z in 0..chf.height {
            for x in 0..chf.width {
                let (start, end) = cell_span_range(chf, x, z);
                for (offset, slots) in neighbours[start..end].iter_mut().enumerate() {
                    for (rdir, slot) in slots.iter_mut().enumerate() {
                        if let Some((_, _, neighbour)) =
                            rc_neighbour(chf, x, z, start + offset, rdir)
                        {
                            // `neighbour < count < u32::MAX`, checked above.
                            *slot = neighbour as u32;
                        }
                    }
                }
            }
        }
        let mut one_way_in = Vec::new();
        for (source, slots) in neighbours.iter().enumerate() {
            for &target in slots {
                if target != NO_SPAN && !neighbours[target as usize].contains(&(source as u32)) {
                    one_way_in.push((target, source as u32));
                }
            }
        }
        one_way_in.sort_unstable();
        Ok(Self {
            neighbours,
            one_way_in,
        })
    }

    /// Spans with a one-way link to `target` (a link `source -> target` that
    /// `target` does not return).
    fn one_way_sources(&self, target: u32) -> impl Iterator<Item = u32> + '_ {
        let first = self.one_way_in.partition_point(|&(t, _)| t < target);
        self.one_way_in[first..]
            .iter()
            .take_while(move |&&(t, _)| t == target)
            .map(|&(_, source)| source)
    }

    /// Flag the spans whose connected component holds at least
    /// `min_region_area` spans.  Smaller components can never survive the
    /// region filter, so they are left out of the watershed (this also keeps
    /// specks from using up the 16-bit region ids).
    fn active_spans(&self, min_region_area: usize) -> Vec<bool> {
        let count = self.neighbours.len();
        if min_region_area <= 1 {
            return vec![true; count];
        }
        let mut parent: Vec<u32> = (0..count as u32).collect();
        for (source, slots) in self.neighbours.iter().enumerate() {
            for &target in slots {
                if target != NO_SPAN {
                    let a = find_root(&mut parent, source as u32);
                    let b = find_root(&mut parent, target);
                    if a != b {
                        parent[a.max(b) as usize] = a.min(b);
                    }
                }
            }
        }
        let mut size = vec![0usize; count];
        for i in 0..count as u32 {
            size[find_root(&mut parent, i) as usize] += 1;
        }
        (0..count as u32)
            .map(|i| size[find_root(&mut parent, i) as usize] >= min_region_area)
            .collect()
    }
}

/// Union-find root of `start` with path halving.
fn find_root(parent: &mut [u32], start: u32) -> u32 {
    let mut node = start;
    while parent[node as usize] != node {
        let grandparent = parent[parent[node as usize] as usize];
        parent[node as usize] = grandparent;
        node = grandparent;
    }
    node
}

/// Spans bucketed by `level / 2`, each bucket in ascending span order.
struct Bands {
    offsets: Vec<usize>,
    items: Vec<u32>,
}

impl Bands {
    fn build(levels: &[u32], active: &[bool], max_level: u32) -> Self {
        let band_count = (max_level / LEVEL_STEP) as usize + 1;
        let mut offsets = vec![0usize; band_count + 1];
        for (&level, &on) in levels.iter().zip(active) {
            if on {
                offsets[(level / LEVEL_STEP) as usize + 1] += 1;
            }
        }
        let mut running = 0usize;
        for slot in offsets.iter_mut() {
            running += *slot;
            *slot = running;
        }
        let mut cursor = offsets.clone();
        let mut items = vec![0u32; running];
        for (i, (&level, &on)) in levels.iter().zip(active).enumerate() {
            if on {
                let band = (level / LEVEL_STEP) as usize;
                items[cursor[band]] = i as u32;
                cursor[band] += 1;
            }
        }
        Self { offsets, items }
    }

    /// The spans whose level lies in `[2 * band, 2 * band + 1]`.
    fn band(&self, band: u32) -> &[u32] {
        let band = band as usize;
        match (self.offsets.get(band), self.offsets.get(band + 1)) {
            (Some(&start), Some(&end)) => &self.items[start..end],
            _ => &[],
        }
    }
}

/// Working state of the watershed.
struct Partition<'a> {
    graph: &'a SpanGraph,
    /// Distance-field level of every span.
    levels: &'a [u32],
    /// Region of every span (`0` = unassigned).
    region: Vec<u16>,
    /// Distance of an expanded span from the seed flood of its region; the
    /// expansion hands a span to the neighbour with the smallest value.
    seed_distance: Vec<u32>,
    /// Next region id to hand out.
    next_region: u32,
    /// Candidate set of the running [`Partition::expand`] call.
    candidates: Marks,
    /// Spans already queued for the next expansion round.
    queued: Marks,
    flood_stack: Vec<u32>,
}

impl<'a> Partition<'a> {
    fn new(graph: &'a SpanGraph, levels: &'a [u32]) -> Self {
        let count = levels.len();
        Self {
            graph,
            levels,
            region: vec![0; count],
            seed_distance: vec![0; count],
            next_region: 1,
            candidates: Marks::new(count),
            queued: Marks::new(count),
            flood_stack: Vec::new(),
        }
    }

    /// First region other than `own` found next to `span`: across one of its
    /// links, or across the link of that neighbour one turn further (the
    /// diagonal).  `0` when no foreign region touches the span.
    fn conflicting_region(&self, span: u32, own: u16) -> u16 {
        for (dir, &neighbour) in self.graph.neighbours[span as usize].iter().enumerate() {
            if neighbour == NO_SPAN {
                continue;
            }
            let adjacent = self.region[neighbour as usize];
            if adjacent != 0 && adjacent != own {
                return adjacent;
            }
            let diagonal = self.graph.neighbours[neighbour as usize][(dir + 1) & 3];
            if diagonal != NO_SPAN {
                let across = self.region[diagonal as usize];
                if across != 0 && across != own {
                    return across;
                }
            }
        }
        0
    }

    /// Hand out the next region id, or fail once all 65535 ids are in use.
    fn allocate_region(&mut self) -> Result<u16, String> {
        let id = u16::try_from(self.next_region).map_err(|_| {
            format!(
                "region id overflow: the heightfield needs more than {} regions",
                u16::MAX
            )
        })?;
        self.next_region += 1;
        Ok(id)
    }

    /// Flood a new region `region` from `seed` over unassigned spans with
    /// `dist >= level - 2`.  A span that touches another region is released
    /// again and left to the expansion.
    fn flood_region(&mut self, seed: u32, level: u32, region: u16) {
        let lowest = level.saturating_sub(LEVEL_STEP);
        let mut stack = std::mem::take(&mut self.flood_stack);
        stack.clear();
        stack.push(seed);
        self.region[seed as usize] = region;
        self.seed_distance[seed as usize] = 0;
        while let Some(current) = stack.pop() {
            if self.conflicting_region(current, region) != 0 {
                self.region[current as usize] = 0;
                continue;
            }
            for &neighbour in &self.graph.neighbours[current as usize] {
                if neighbour == NO_SPAN {
                    continue;
                }
                let n = neighbour as usize;
                if self.levels[n] >= lowest && self.region[n] == 0 {
                    self.region[n] = region;
                    self.seed_distance[n] = 0;
                    stack.push(neighbour);
                }
            }
        }
        self.flood_stack = stack;
    }

    /// The region (and expansion distance) `span` would take from its linked
    /// neighbours: the neighbouring region whose span is nearest to its seed.
    fn nearest_region(&self, span: u32) -> Option<(u16, u32)> {
        let mut best: Option<(u16, u32)> = None;
        for &neighbour in &self.graph.neighbours[span as usize] {
            if neighbour == NO_SPAN {
                continue;
            }
            let region = self.region[neighbour as usize];
            if region == 0 {
                continue;
            }
            let distance = self.seed_distance[neighbour as usize].saturating_add(LEVEL_STEP);
            if best.is_none_or(|(_, d)| distance < d) {
                best = Some((region, distance));
            }
        }
        best
    }

    /// Queue `span` for the next expansion round when it is an unassigned
    /// candidate that has not been queued yet.
    fn queue_candidate(&mut self, span: u32, next: &mut Vec<u32>) {
        if span == NO_SPAN {
            return;
        }
        let s = span as usize;
        if self.region[s] == 0 && self.candidates.is_marked(s) && !self.queued.is_marked(s) {
            self.queued.mark(s);
            next.push(span);
        }
    }

    /// Grow the existing regions into the unassigned spans of `candidates`.
    ///
    /// Each round every candidate that links to an assigned span takes the
    /// region of its nearest such neighbour; the assignments of one round are
    /// applied together, so the result does not depend on the candidate order.
    /// Only the spans linking to a span assigned in the previous round are
    /// examined again.  For `level > 0` the growth stops after
    /// [`EXPAND_ROUNDS`] rounds, at level 0 it runs until nothing changes.
    fn expand(&mut self, candidates: &[u32], level: u32) {
        self.candidates.begin();
        let mut frontier: Vec<u32> = Vec::with_capacity(candidates.len());
        for &span in candidates {
            if self.region[span as usize] == 0 {
                self.candidates.mark(span as usize);
                frontier.push(span);
            }
        }
        let mut assigned: Vec<(u32, u16, u32)> = Vec::new();
        let mut next: Vec<u32> = Vec::new();
        let mut rounds = 0;
        while !frontier.is_empty() {
            assigned.clear();
            for &span in &frontier {
                if let Some((region, distance)) = self.nearest_region(span) {
                    assigned.push((span, region, distance));
                }
            }
            if assigned.is_empty() {
                break;
            }
            for &(span, region, distance) in &assigned {
                self.region[span as usize] = region;
                self.seed_distance[span as usize] = distance;
            }
            self.queued.begin();
            next.clear();
            for &(span, _, _) in &assigned {
                let graph = self.graph;
                for &neighbour in &graph.neighbours[span as usize] {
                    self.queue_candidate(neighbour, &mut next);
                }
                for source in graph.one_way_sources(span) {
                    self.queue_candidate(source, &mut next);
                }
            }
            std::mem::swap(&mut frontier, &mut next);
            if level > 0 {
                rounds += 1;
                if rounds >= EXPAND_ROUNDS {
                    break;
                }
            }
        }
    }

    /// Run the level loop over the active spans.
    fn grow(&mut self, active: &[bool], bands: &Bands, max_level: u32) -> Result<(), String> {
        // Like Recast, start at the smallest even level that is not below the
        // maximum and step down by two before every pass.  The loop always runs
        // the level-0 pass, which also partitions fields whose distance is 0
        // everywhere (one-cell-wide corridors).
        let mut level = (max_level + 1) & !1;
        let mut first = true;
        let mut stack: Vec<u32> = Vec::new();
        let mut leftover: Vec<u32> = Vec::new();
        loop {
            level = level.saturating_sub(LEVEL_STEP);
            let band = level / LEVEL_STEP;
            stack.clear();
            if first {
                // The first pass takes every span at or above its level.
                for (i, &on) in active.iter().enumerate() {
                    if on && self.region[i] == 0 && self.levels[i] / LEVEL_STEP >= band {
                        stack.push(i as u32);
                    }
                }
                first = false;
            } else {
                stack.extend(
                    bands
                        .band(band)
                        .iter()
                        .copied()
                        .filter(|&s| self.region[s as usize] == 0),
                );
            }
            stack.extend(
                leftover
                    .iter()
                    .copied()
                    .filter(|&s| self.region[s as usize] == 0),
            );
            self.expand(&stack, level);
            for &seed in &stack {
                if self.region[seed as usize] != 0 || self.conflicting_region(seed, 0) != 0 {
                    continue;
                }
                let region = self.allocate_region()?;
                self.flood_region(seed, level, region);
            }
            leftover.clear();
            leftover.extend(
                stack
                    .iter()
                    .copied()
                    .filter(|&s| self.region[s as usize] == 0),
            );
            if level == 0 {
                break;
            }
        }
        self.expand(&leftover, 0);
        Ok(())
    }
}

/// Partition the walkable spans of `chf` into watershed regions.
///
/// Returns one region id per span (`0` = not part of any region).  Region ids
/// are `1..=n` in order of creation.  Connected groups of fewer than
/// `min_region_area` spans are skipped.  Fails when the field needs more than
/// 65535 regions or the compact heightfield is inconsistent.
pub(super) fn partition(
    chf: &CompactHeightfield,
    min_region_area: usize,
) -> Result<Vec<u16>, String> {
    let count = chf.spans.len();
    if count == 0 {
        return Ok(Vec::new());
    }
    if chf.dist.len() != count {
        return Err(format!(
            "compact heightfield has {count} spans but {} distance values",
            chf.dist.len()
        ));
    }
    let graph = SpanGraph::build(chf)?;
    let active = graph.active_spans(min_region_area);
    // A span without a distance (`u16::MAX`) is treated as a border span.
    let levels: Vec<u32> = chf
        .dist
        .iter()
        .map(|&d| if d == u16::MAX { 0 } else { u32::from(d) })
        .collect();
    let max_level = levels
        .iter()
        .zip(&active)
        .filter(|&(_, &on)| on)
        .map(|(&level, _)| level)
        .max();
    let Some(max_level) = max_level else {
        return Ok(vec![0; count]);
    };
    let bands = Bands::build(&levels, &active, max_level);
    let mut state = Partition::new(&graph, &levels);
    state.grow(&active, &bands, max_level)?;
    Ok(state.region)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recast::RecastConfig;
    use crate::recast::compact::{CompactCell, CompactSpan, build_compact_heightfield};
    use crate::recast::rasterize::Heightfield;
    use std::collections::BTreeSet;

    fn test_config() -> RecastConfig {
        RecastConfig {
            cell_size: 1.0,
            cell_height: 0.1,
            agent_height: 1.0,
            agent_max_climb: 0.2,
            ..Default::default()
        }
    }

    /// Single-layer floor on a `width x depth` grid; `walkable(x, z)` selects the
    /// cells that carry a floor span.
    fn floor(
        width: usize,
        depth: usize,
        walkable: impl Fn(usize, usize) -> bool,
    ) -> CompactHeightfield {
        let mut hf = Heightfield::new(
            width,
            depth,
            [0.0; 3],
            [width as f64, 1.0, depth as f64],
            1.0,
            0.1,
        )
        .expect("heightfield");
        for z in 0..depth {
            for x in 0..width {
                if walkable(x, z) {
                    hf.add_span(x, z, 0, 1, true);
                }
            }
        }
        build_compact_heightfield(&hf, &test_config())
    }

    /// Region of the floor span in cell `(x, z)`.
    fn region_at(chf: &CompactHeightfield, regions: &[u16], x: usize, z: usize) -> u16 {
        let (start, end) = cell_span_range(chf, x, z);
        assert_eq!(end - start, 1, "cell ({x}, {z}) must hold exactly one span");
        regions[start]
    }

    /// Distinct regions found in the walkable cells selected by `cells`.
    fn regions_in(
        chf: &CompactHeightfield,
        regions: &[u16],
        cells: impl Fn(usize, usize) -> bool,
    ) -> BTreeSet<u16> {
        let mut found = BTreeSet::new();
        for z in 0..chf.height {
            for x in 0..chf.width {
                if cells(x, z) {
                    found.insert(region_at(chf, regions, x, z));
                }
            }
        }
        found
    }

    /// Deterministic xorshift generator.
    fn next_random(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    /// Isolated floor cells on a 512 x 256 checkerboard (`count` of them).
    fn checkerboard(count: usize) -> CompactHeightfield {
        let placed = std::cell::Cell::new(0usize);
        floor(512, 256, |x, z| {
            if (x + z).is_multiple_of(2) && placed.get() < count {
                placed.set(placed.get() + 1);
                true
            } else {
                false
            }
        })
    }

    /// A row of single-span cells.  `east[i]` / `west[i]` say whether span `i`
    /// links to its right / left neighbour; `dist` is the distance field.
    fn line(east: &[bool], west: &[bool], dist: &[u16]) -> CompactHeightfield {
        let n = east.len();
        let slot = |linked: bool| if linked { 0u32 } else { 0xFF };
        let spans: Vec<CompactSpan> = (0..n)
            .map(|i| CompactSpan {
                y: 0,
                reg: 0,
                con: slot(east[i]) | (0xFF << 8) | (slot(west[i]) << 16) | (0xFF << 24),
                h: 10,
            })
            .collect();
        CompactHeightfield {
            width: n,
            height: 1,
            span_count: n,
            walkable_height: 10,
            walkable_climb: 2,
            bmin: [0.0; 3],
            bmax: [n as f64, 1.0, 1.0],
            cs: 1.0,
            ch: 0.1,
            cells: (0..n)
                .map(|i| CompactCell {
                    index: i as u32,
                    count: 1,
                })
                .collect(),
            spans,
            dist: dist.to_vec(),
            max_dist: dist.iter().copied().max().unwrap_or(0),
        }
    }

    #[test]
    fn empty_heightfield_has_no_regions() {
        let chf = floor(4, 4, |_, _| false);
        assert!(chf.spans.is_empty());
        assert_eq!(partition(&chf, 0).expect("partition"), Vec::<u16>::new());
    }

    #[test]
    fn flat_room_is_a_single_region() {
        let chf = floor(20, 20, |_, _| true);
        let regions = partition(&chf, 0).expect("partition");
        assert_eq!(regions.len(), chf.spans.len());
        let ids: BTreeSet<u16> = regions.iter().copied().collect();
        assert_eq!(ids, BTreeSet::from([1]), "one basin gives one region");
    }

    #[test]
    fn long_rectangular_room_is_a_single_region() {
        // The distance ridge is a plateau 30 cells long: a single flood covers it.
        let chf = floor(42, 12, |x, z| {
            (1..=40).contains(&x) && (1..=10).contains(&z)
        });
        let regions = partition(&chf, 0).expect("partition");
        let ids: BTreeSet<u16> = regions.iter().copied().collect();
        assert_eq!(ids, BTreeSet::from([1]));
    }

    #[test]
    fn rooms_joined_by_a_doorway_split_into_two_regions() {
        let left = |x: usize, z: usize| (1..=12).contains(&x) && (1..=12).contains(&z);
        let right = |x: usize, z: usize| (16..=27).contains(&x) && (1..=12).contains(&z);
        let door = |x: usize, z: usize| (13..=15).contains(&x) && z == 6;
        let chf = floor(29, 14, |x, z| left(x, z) || right(x, z) || door(x, z));
        // The two rooms and the doorway form a single connected component, so
        // the old flood fill produced exactly one region for all of it.
        assert_eq!(chf.spans.len(), 2 * 144 + 3);
        let regions = partition(&chf, 0).expect("partition");
        let in_left = regions_in(&chf, &regions, left);
        let in_right = regions_in(&chf, &regions, right);
        assert_eq!(in_left.len(), 1, "left room: {in_left:?}");
        assert_eq!(in_right.len(), 1, "right room: {in_right:?}");
        assert!(in_left.is_disjoint(&in_right), "{in_left:?} / {in_right:?}");
        assert!(!in_left.contains(&0) && !in_right.contains(&0));
        // Every door cell joins one of the two rooms.
        let rooms: BTreeSet<u16> = in_left.union(&in_right).copied().collect();
        assert!(regions_in(&chf, &regions, door).is_subset(&rooms));
    }

    #[test]
    fn regions_split_at_the_neck_of_a_dumbbell() {
        // Two 9 x 9 chambers joined by a 3-wide, 8-long neck.
        let chamber = |x: usize, z: usize| {
            (1..=9).contains(&z) && ((1..=9).contains(&x) || (18..=26).contains(&x))
        };
        let neck = |x: usize, z: usize| (10..=17).contains(&x) && (4..=6).contains(&z);
        let chf = floor(28, 11, |x, z| chamber(x, z) || neck(x, z));
        let regions = partition(&chf, 0).expect("partition");
        let west = regions_in(&chf, &regions, |x, z| chamber(x, z) && x <= 9);
        let east = regions_in(&chf, &regions, |x, z| chamber(x, z) && x >= 18);
        assert_eq!(west.len(), 1, "{west:?}");
        assert_eq!(east.len(), 1, "{east:?}");
        assert!(west.is_disjoint(&east));
    }

    #[test]
    fn obstacle_field_regions_are_complete_contiguous_and_connected() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let (width, depth) = (60, 48);
        let blocked: Vec<bool> = (0..width * depth)
            .map(|_| next_random(&mut state).is_multiple_of(7))
            .collect();
        let chf = floor(width, depth, |x, z| !blocked[z * width + x]);
        let regions = partition(&chf, 0).expect("partition");
        assert!(
            regions.iter().all(|&r| r != 0),
            "every walkable span belongs to a region"
        );
        let ids: BTreeSet<u16> = regions.iter().copied().collect();
        assert_eq!(ids.first().copied(), Some(1));
        assert_eq!(ids.last().copied(), u16::try_from(ids.len()).ok());
        assert!(ids.len() > 1, "obstacles must produce several regions");
        // Each region is one connected piece of the span graph.
        let graph = SpanGraph::build(&chf).expect("graph");
        for &id in &ids {
            let members: Vec<usize> = (0..regions.len()).filter(|&i| regions[i] == id).collect();
            let mut seen = vec![false; regions.len()];
            let mut stack = vec![members[0]];
            seen[members[0]] = true;
            let mut reached = 0;
            while let Some(current) = stack.pop() {
                reached += 1;
                for &next in &graph.neighbours[current] {
                    if next != NO_SPAN && regions[next as usize] == id && !seen[next as usize] {
                        seen[next as usize] = true;
                        stack.push(next as usize);
                    }
                }
            }
            assert_eq!(reached, members.len(), "region {id} is split in pieces");
        }
    }

    #[test]
    fn partition_is_deterministic() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let blocked: Vec<bool> = (0..40 * 40)
            .map(|_| next_random(&mut state).is_multiple_of(5))
            .collect();
        let chf = floor(40, 40, |x, z| !blocked[z * 40 + x]);
        let first = partition(&chf, 4).expect("partition");
        let second = partition(&chf, 4).expect("partition");
        assert_eq!(first, second);
    }

    #[test]
    fn one_cell_wide_corridor_is_still_partitioned() {
        let chf = floor(24, 3, |_, z| z == 1);
        assert_eq!(
            chf.max_dist, 0,
            "every span of a 1-wide corridor is a border"
        );
        let regions = partition(&chf, 0).expect("partition");
        let ids: BTreeSet<u16> = regions.iter().copied().collect();
        assert_eq!(ids, BTreeSet::from([1]));
    }

    #[test]
    fn components_below_the_minimum_area_are_skipped() {
        // An 8 x 8 floor and a 2 x 2 speck.
        let floor_cell = |x: usize, z: usize| (1..=8).contains(&x) && (1..=8).contains(&z);
        let speck = |x: usize, z: usize| (12..=13).contains(&x) && (4..=5).contains(&z);
        let chf = floor(16, 10, |x, z| floor_cell(x, z) || speck(x, z));
        let filtered = partition(&chf, 16).expect("partition");
        assert_eq!(regions_in(&chf, &filtered, speck), BTreeSet::from([0]));
        assert_eq!(regions_in(&chf, &filtered, floor_cell), BTreeSet::from([1]));
        let unfiltered = partition(&chf, 0).expect("partition");
        let kept = regions_in(&chf, &unfiltered, speck);
        assert_eq!(kept.len(), 1);
        assert!(!kept.contains(&0));
        assert!(kept.is_disjoint(&regions_in(&chf, &unfiltered, floor_cell)));
    }

    #[test]
    fn region_ids_are_exhausted_cleanly() {
        let full = checkerboard(65_536);
        assert_eq!(full.spans.len(), 65_536);
        let err = partition(&full, 0).expect_err("65536 regions do not fit 16-bit ids");
        assert!(err.contains("region id overflow"), "{err}");
        // One span less needs exactly the 65535 ids that exist.
        let fits = checkerboard(65_535);
        let regions = partition(&fits, 0).expect("65535 regions fit");
        let ids: BTreeSet<u16> = regions.iter().copied().collect();
        assert_eq!(ids.len(), 65_535);
        assert_eq!(ids.first().copied(), Some(1));
        assert_eq!(ids.last().copied(), Some(u16::MAX));
    }

    #[test]
    fn expansion_follows_links_that_are_not_mutual() {
        // s0 <-> s1 are linked both ways; s2 -> s1, s3 -> s2 and s4 -> s3 are
        // one-way links.  s2..s4 all sit at level 0, so they are assigned in one
        // expansion call that has to chain through the one-way links.
        let east = [true, false, false, false, false];
        let west = [false, true, true, true, true];
        let chf = line(&east, &west, &[8, 6, 0, 0, 0]);
        let graph = SpanGraph::build(&chf).expect("graph");
        assert_eq!(graph.one_way_in, vec![(1, 2), (2, 3), (3, 4)]);
        assert_eq!(graph.one_way_sources(2).collect::<Vec<_>>(), vec![3]);
        assert_eq!(graph.one_way_sources(0).count(), 0);
        let regions = partition(&chf, 0).expect("partition");
        assert_eq!(regions, vec![1; 5]);
    }

    #[test]
    fn stacked_floors_with_a_one_way_link_stay_in_one_region() {
        // Cell 0 holds two floors (tops at 1 and 4), cell 1 one floor (top 2).
        // Both floors of cell 0 can step onto cell 1, but cell 1 links back to
        // the lower floor only, so the upper floor has a one-way link.
        let mut hf =
            Heightfield::new(2, 1, [0.0; 3], [2.0, 1.0, 1.0], 1.0, 0.1).expect("heightfield");
        hf.add_span(0, 0, 0, 1, true);
        hf.add_span(0, 0, 3, 4, true);
        hf.add_span(1, 0, 0, 2, true);
        let chf = build_compact_heightfield(&hf, &test_config());
        assert_eq!(chf.spans.len(), 3);
        let graph = SpanGraph::build(&chf).expect("graph");
        assert_eq!(graph.one_way_in, vec![(2, 1)]);
        assert_eq!(partition(&chf, 3).expect("partition"), vec![1, 1, 1]);
        // The three spans form one component, so a larger minimum removes it.
        assert_eq!(partition(&chf, 4).expect("partition"), vec![0, 0, 0]);
    }

    #[test]
    fn marks_survive_generation_wraparound() {
        let mut marks = Marks::new(3);
        marks.generation = u32::MAX - 1;
        marks.begin();
        assert_eq!(marks.generation, u32::MAX);
        marks.mark(1);
        assert!(marks.is_marked(1) && !marks.is_marked(0));
        marks.begin();
        assert_eq!(marks.generation, 1, "the counter restarts after clearing");
        assert!(!marks.is_marked(1), "old marks are gone after the wrap");
        marks.mark(0);
        assert!(marks.is_marked(0) && !marks.is_marked(1));
    }

    #[test]
    fn inconsistent_distance_field_is_rejected() {
        let mut chf = floor(6, 6, |_, _| true);
        chf.dist.pop();
        let err = partition(&chf, 0).expect_err("length mismatch");
        assert!(err.contains("distance values"), "{err}");
    }

    #[test]
    fn dangling_connections_are_ignored() {
        // Every span claims a slot-0 neighbour in all four directions, which
        // points out of the grid along its border.
        let mut chf = floor(6, 6, |_, _| true);
        for span in &mut chf.spans {
            span.con = 0;
        }
        let regions = partition(&chf, 0).expect("partition");
        assert!(regions.iter().all(|&r| r != 0));
    }

    #[test]
    fn spans_without_a_distance_are_treated_as_borders() {
        let mut chf = floor(8, 8, |_, _| true);
        chf.dist.fill(u16::MAX);
        let regions = partition(&chf, 0).expect("partition");
        let ids: BTreeSet<u16> = regions.iter().copied().collect();
        assert_eq!(ids, BTreeSet::from([1]));
    }
}
