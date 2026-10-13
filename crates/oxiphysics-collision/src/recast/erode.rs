// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Agent-radius erosion of the walkable area (port of Recast's
//! `rcErodeWalkableArea`).
//!
//! Every walkable span closer than `agent_radius` to an obstacle or ledge is
//! made unwalkable, so that an agent of that radius following the resulting
//! navmesh never clips into walls.  The distance is a two-pass chamfer
//! transform (orthogonal step = 2, diagonal step = 3) measured from spans that
//! lack at least one walkable neighbour; spans with distance `< 2 * radius`
//! are removed.

use super::compact::CompactHeightfield;
use super::grid::{cell_span_range, rc_neighbour};
use super::rasterize::Heightfield;

/// Walkable radius in cells for an agent of radius `agent_radius`
/// (`ceil(agent_radius / cell_size)`, clamped to the `u16` range).
pub fn walkable_radius_cells(agent_radius: f64, cell_size: f64) -> u16 {
    if !(agent_radius.is_finite() && cell_size.is_finite()) || cell_size <= 0.0 {
        return 0;
    }
    let cells = (agent_radius / cell_size).ceil();
    if cells <= 0.0 {
        0
    } else if cells >= f64::from(u16::MAX) {
        u16::MAX
    } else {
        cells as u16
    }
}

/// Chamfer distance (orthogonal = 2, diagonal = 3) of every compact span to
/// the nearest span that is missing one of its four walkable neighbours.
pub fn erosion_distance_field(chf: &CompactHeightfield) -> Vec<u16> {
    let w = chf.width;
    let h = chf.height;
    let n = chf.spans.len();
    let mut dist = vec![u16::MAX; n];

    // Mark boundary spans.
    for z in 0..h {
        for x in 0..w {
            let (start, end) = cell_span_range(chf, x, z);
            for (i, d) in dist.iter_mut().enumerate().take(end).skip(start) {
                let connected = (0..4)
                    .filter(|&rdir| rc_neighbour(chf, x, z, i, rdir).is_some())
                    .count();
                if connected != 4 {
                    *d = 0;
                }
            }
        }
    }

    // Relax `dist[i]` from neighbour `ni` with step cost `cost`.
    let relax = |dist: &mut [u16], i: usize, ni: usize, cost: u16| {
        let nd = dist[ni].saturating_add(cost);
        if nd < dist[i] {
            dist[i] = nd;
        }
    };

    // Pass 1: neighbours (-1, 0), (-1, -1), (0, -1), (1, -1).
    for z in 0..h {
        for x in 0..w {
            let (start, end) = cell_span_range(chf, x, z);
            for i in start..end {
                if let Some((ax, az, ai)) = rc_neighbour(chf, x, z, i, 0) {
                    relax(&mut dist, i, ai, 2);
                    if let Some((_, _, aai)) = rc_neighbour(chf, ax, az, ai, 3) {
                        relax(&mut dist, i, aai, 3);
                    }
                }
                if let Some((ax, az, ai)) = rc_neighbour(chf, x, z, i, 3) {
                    relax(&mut dist, i, ai, 2);
                    if let Some((_, _, aai)) = rc_neighbour(chf, ax, az, ai, 2) {
                        relax(&mut dist, i, aai, 3);
                    }
                }
            }
        }
    }

    // Pass 2: neighbours (1, 0), (1, 1), (0, 1), (-1, 1).
    for z in (0..h).rev() {
        for x in (0..w).rev() {
            let (start, end) = cell_span_range(chf, x, z);
            for i in start..end {
                if let Some((ax, az, ai)) = rc_neighbour(chf, x, z, i, 2) {
                    relax(&mut dist, i, ai, 2);
                    if let Some((_, _, aai)) = rc_neighbour(chf, ax, az, ai, 1) {
                        relax(&mut dist, i, aai, 3);
                    }
                }
                if let Some((ax, az, ai)) = rc_neighbour(chf, x, z, i, 1) {
                    relax(&mut dist, i, ai, 2);
                    if let Some((_, _, aai)) = rc_neighbour(chf, ax, az, ai, 0) {
                        relax(&mut dist, i, aai, 3);
                    }
                }
            }
        }
    }

    dist
}

/// Erode the walkable area of `hf` by `radius` cells.
///
/// `chf` must be the compact heightfield built from `hf` (it supplies the
/// span connectivity).  Eroded spans are marked non-walkable in `hf`; the
/// caller rebuilds the compact heightfield afterwards.  Returns the number of
/// spans that were eroded.
pub fn erode_walkable_area(hf: &mut Heightfield, chf: &CompactHeightfield, radius: u16) -> usize {
    if radius == 0 || chf.spans.is_empty() {
        return 0;
    }
    let dist = erosion_distance_field(chf);
    let threshold = u32::from(radius) * 2;

    let mut eroded = 0usize;
    for z in 0..chf.height.min(hf.height) {
        for x in 0..chf.width.min(hf.width) {
            let (start, end) = cell_span_range(chf, x, z);
            if start == end {
                continue;
            }
            let col_idx = hf.col_idx(x, z);
            let Some(column) = hf.spans.get_mut(col_idx) else {
                continue;
            };
            // The k-th compact span of a cell is the k-th walkable span of
            // the corresponding open-heightfield column.
            let mut walkable_spans = column.iter_mut().filter(|s| s.walkable);
            for &d in &dist[start..end] {
                let Some(span) = walkable_spans.next() else {
                    break;
                };
                if u32::from(d) < threshold {
                    span.walkable = false;
                    eroded += 1;
                }
            }
        }
    }
    eroded
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recast::RecastConfig;
    use crate::recast::compact::build_compact_heightfield;

    /// An `n x n` block of walkable cells surrounded by one empty cell.
    fn square_heightfield(n: usize) -> Heightfield {
        let size = n + 2;
        let mut hf = Heightfield::new(
            size,
            size,
            [0.0; 3],
            [size as f64, 1.0, size as f64],
            1.0,
            0.1,
        )
        .expect("small heightfield");
        for z in 1..=n {
            for x in 1..=n {
                hf.add_span(x, z, 0, 1, true);
            }
        }
        hf
    }

    fn cfg() -> RecastConfig {
        RecastConfig {
            cell_size: 1.0,
            cell_height: 0.1,
            agent_height: 1.0,
            agent_max_climb: 0.2,
            ..Default::default()
        }
    }

    #[test]
    fn test_walkable_radius_cells_rounds_up() {
        assert_eq!(walkable_radius_cells(0.6, 0.3), 2);
        assert_eq!(walkable_radius_cells(0.61, 0.3), 3);
        assert_eq!(walkable_radius_cells(0.0, 0.3), 0);
        assert_eq!(walkable_radius_cells(f64::NAN, 0.3), 0);
    }

    #[test]
    fn test_distance_field_rings() {
        let hf = square_heightfield(7);
        let chf = build_compact_heightfield(&hf, &cfg());
        let dist = erosion_distance_field(&chf);
        assert_eq!(chf.spans.len(), 49);
        // Centre cell (4, 4) is three rings in: distance 6.
        let (start, _) = cell_span_range(&chf, 4, 4);
        assert_eq!(dist[start], 6);
        // Corner cell is on the boundary.
        let (corner, _) = cell_span_range(&chf, 1, 1);
        assert_eq!(dist[corner], 0);
        // Second ring.
        let (ring1, _) = cell_span_range(&chf, 2, 3);
        assert_eq!(dist[ring1], 2);
    }

    #[test]
    fn test_erosion_shrinks_square_by_radius() {
        let mut hf = square_heightfield(10);
        let chf = build_compact_heightfield(&hf, &cfg());
        assert_eq!(chf.spans.len(), 100);
        let eroded = erode_walkable_area(&mut hf, &chf, 2);
        assert_eq!(eroded, 100 - 36);
        let chf2 = build_compact_heightfield(&hf, &cfg());
        assert_eq!(chf2.spans.len(), 36);
        // The surviving block is the inner 6x6 square: cells 3..=8.
        for z in 0..hf.height {
            for x in 0..hf.width {
                let (s, e) = cell_span_range(&chf2, x, z);
                let inside = (3..=8).contains(&x) && (3..=8).contains(&z);
                assert_eq!(e - s, usize::from(inside), "cell ({x},{z})");
            }
        }
    }

    #[test]
    fn test_erosion_removes_narrow_corridor() {
        // A 3-cell-wide corridor cannot fit an agent of radius 2 cells.
        let mut hf = Heightfield::new(22, 5, [0.0; 3], [22.0, 1.0, 5.0], 1.0, 0.1)
            .expect("small heightfield");
        for z in 1..=3 {
            for x in 1..=20 {
                hf.add_span(x, z, 0, 1, true);
            }
        }
        let chf = build_compact_heightfield(&hf, &cfg());
        let eroded = erode_walkable_area(&mut hf, &chf, 2);
        assert_eq!(eroded, 60);
        let chf2 = build_compact_heightfield(&hf, &cfg());
        assert!(chf2.spans.is_empty());
    }

    #[test]
    fn test_zero_radius_is_noop() {
        let mut hf = square_heightfield(4);
        let chf = build_compact_heightfield(&hf, &cfg());
        assert_eq!(erode_walkable_area(&mut hf, &chf, 0), 0);
        assert_eq!(build_compact_heightfield(&hf, &cfg()).spans.len(), 16);
    }
}
