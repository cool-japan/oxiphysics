// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared helpers for the Recast stages that walk region boundaries.
//!
//! The compact heightfield stores its connectivity using the local direction
//! table [`DIRS`](super::compact::DIRS) (`0 = +X, 1 = +Z, 2 = -X, 3 = -Z`).
//! The boundary-walking algorithms (region merging, contour tracing) and the
//! integer polygon predicates (hole merging, ear clipping) are ports of Recast,
//! whose orientation conventions depend on Recast's own direction order
//! (`0 = -X, 1 = +Z, 2 = +X, 3 = -Z`, i.e. clockwise rotation when the walk
//! rotates `dir -> dir + 1`).  This module exposes the connectivity in Recast's
//! order so that the ported algorithms keep their exact geometric meaning.

use super::compact::{CompactHeightfield, CompactSpan};

/// X offset of each Recast direction (`0 = -X, 1 = +Z, 2 = +X, 3 = -Z`).
pub(crate) const RC_DIR_OFFSET_X: [i64; 4] = [-1, 0, 1, 0];
/// Z offset of each Recast direction (`0 = -X, 1 = +Z, 2 = +X, 3 = -Z`).
pub(crate) const RC_DIR_OFFSET_Z: [i64; 4] = [0, 1, 0, -1];

/// Maps a Recast direction to the compact heightfield's local direction slot.
const RC_TO_COMPACT_DIR: [usize; 4] = [2, 1, 0, 3];

/// Sentinel stored in a connection slot when there is no neighbour.
const NOT_CONNECTED: u32 = 0xFF;

/// Connection slot (index inside the neighbour cell) of `span` towards Recast
/// direction `rdir`, or `None` when the span has no neighbour that way.
pub(crate) fn rc_con(span: &CompactSpan, rdir: usize) -> Option<usize> {
    let local = RC_TO_COMPACT_DIR[rdir & 3];
    let raw = (span.con >> (local * 8)) & 0xFF;
    if raw == NOT_CONNECTED {
        None
    } else {
        Some(raw as usize)
    }
}

/// Neighbour of span `span_idx` (located in cell `(x, z)`) across Recast
/// direction `rdir`, returned as `(nx, nz, neighbour_span_index)`.
///
/// Every index is bounds-checked, so a corrupted connection value yields
/// `None` instead of an out-of-bounds access.
pub(crate) fn rc_neighbour(
    chf: &CompactHeightfield,
    x: usize,
    z: usize,
    span_idx: usize,
    rdir: usize,
) -> Option<(usize, usize, usize)> {
    let span = chf.spans.get(span_idx)?;
    let slot = rc_con(span, rdir)?;
    let nx = x as i64 + RC_DIR_OFFSET_X[rdir & 3];
    let nz = z as i64 + RC_DIR_OFFSET_Z[rdir & 3];
    if nx < 0 || nz < 0 || nx >= chf.width as i64 || nz >= chf.height as i64 {
        return None;
    }
    let (nx, nz) = (nx as usize, nz as usize);
    let cell = chf.cells.get(nz * chf.width + nx)?;
    if slot >= cell.count as usize {
        return None;
    }
    let ni = cell.index as usize + slot;
    if ni >= chf.spans.len() {
        return None;
    }
    Some((nx, nz, ni))
}

/// Span index range `[start, end)` of cell `(x, z)`, clamped to the span array.
pub(crate) fn cell_span_range(chf: &CompactHeightfield, x: usize, z: usize) -> (usize, usize) {
    match chf.cells.get(z * chf.width + x) {
        Some(cell) => {
            let start = (cell.index as usize).min(chf.spans.len());
            let end = (start + cell.count as usize).min(chf.spans.len());
            (start, end)
        }
        None => (0, 0),
    }
}

/// A 2-D integer point `(x, z)` used by the polygon predicates.
pub(crate) type Pt = [i64; 2];

/// Twice the signed area of triangle `(a, b, c)` in Recast's convention.
#[inline]
pub(crate) fn area2(a: Pt, b: Pt, c: Pt) -> i64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
}

/// `c` lies strictly to the left of the directed line `a -> b` (Recast convention).
#[inline]
pub(crate) fn left(a: Pt, b: Pt, c: Pt) -> bool {
    area2(a, b, c) < 0
}

/// `c` lies to the left of or on the directed line `a -> b` (Recast convention).
#[inline]
pub(crate) fn left_on(a: Pt, b: Pt, c: Pt) -> bool {
    area2(a, b, c) <= 0
}

/// `a`, `b` and `c` are collinear.
#[inline]
pub(crate) fn collinear(a: Pt, b: Pt, c: Pt) -> bool {
    area2(a, b, c) == 0
}

/// Segments `ab` and `cd` properly intersect (share an interior point and are
/// not collinear at any endpoint).
pub(crate) fn intersect_prop(a: Pt, b: Pt, c: Pt, d: Pt) -> bool {
    if collinear(a, b, c) || collinear(a, b, d) || collinear(c, d, a) || collinear(c, d, b) {
        return false;
    }
    (left(a, b, c) ^ left(a, b, d)) && (left(c, d, a) ^ left(c, d, b))
}

/// `c` lies on the closed segment `ab` (requires collinearity).
pub(crate) fn between(a: Pt, b: Pt, c: Pt) -> bool {
    if !collinear(a, b, c) {
        return false;
    }
    if a[0] != b[0] {
        (a[0] <= c[0] && c[0] <= b[0]) || (a[0] >= c[0] && c[0] >= b[0])
    } else {
        (a[1] <= c[1] && c[1] <= b[1]) || (a[1] >= c[1] && c[1] >= b[1])
    }
}

/// Segments `ab` and `cd` intersect, properly or improperly.
pub(crate) fn intersect(a: Pt, b: Pt, c: Pt, d: Pt) -> bool {
    intersect_prop(a, b, c, d)
        || between(a, b, c)
        || between(a, b, d)
        || between(c, d, a)
        || between(c, d, b)
}

/// Previous index in a cyclic sequence of length `n` (`n > 0`).
#[inline]
pub(crate) fn prev(i: usize, n: usize) -> usize {
    if i == 0 { n - 1 } else { i - 1 }
}

/// Next index in a cyclic sequence of length `n` (`n > 0`).
#[inline]
pub(crate) fn next(i: usize, n: usize) -> usize {
    if i + 1 < n { i + 1 } else { 0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_area2_orientation_matches_recast() {
        // Recast's convention: a clockwise walk in (x, z) has negative area2 turns.
        assert!(left([0, 0], [0, 1], [1, 1]));
        assert!(!left([0, 0], [0, 1], [-1, 1]));
        assert!(left_on([0, 0], [1, 0], [2, 0]));
        assert!(collinear([0, 0], [1, 1], [3, 3]));
    }

    #[test]
    fn test_segment_intersection_predicates() {
        assert!(intersect_prop([0, 0], [2, 2], [0, 2], [2, 0]));
        assert!(!intersect_prop([0, 0], [1, 0], [2, 0], [3, 0]));
        assert!(between([0, 0], [4, 0], [2, 0]));
        assert!(!between([0, 0], [4, 0], [5, 0]));
        assert!(intersect([0, 0], [4, 0], [2, 0], [2, 3]));
        assert!(!intersect([0, 0], [1, 0], [0, 1], [1, 1]));
    }

    #[test]
    fn test_cyclic_index_helpers() {
        assert_eq!(prev(0, 5), 4);
        assert_eq!(next(4, 5), 0);
        assert_eq!(next(1, 5), 2);
        assert_eq!(prev(3, 5), 2);
    }
}
