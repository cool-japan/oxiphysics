// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Marching-cubes lookup tables and cube topology shared by every isosurface
//! extractor of this crate.
//!
//! The tables are the classic Lorensen & Cline / Bourke tables for the corner
//! and edge numbering below:
//!
//! ```text
//!        7 ---------- 6          corner  offset (di, dj, dk)
//!       /|           /|            0      (0, 0, 0)
//!      4 ---------- 5 |            1      (1, 0, 0)
//!      | |          | |            2      (1, 1, 0)
//!      | 3 ---------|-2            3      (0, 1, 0)
//!      |/           |/             4      (0, 0, 1)
//!      0 ---------- 1              5      (1, 0, 1)
//!                                  6      (1, 1, 1)
//!                                  7      (0, 1, 1)
//! ```
//!
//! A corner is *inside* when its field value is below the iso-value, and bit `c`
//! of the 8-bit case index is set for an inside corner `c` (see [`case_index`]).
//! Edge `e` joins the corners listed in [`MC_EDGE_CORNERS`]; bit `e` of
//! [`MC_EDGE_TABLE`]`[case]` is set when that edge is crossed by the surface, and
//! [`MC_TRI_TABLE`]`[case]` lists the triangles of the case as consecutive
//! triples of edge indices, terminated by `-1` (see [`case_triangles`]).
//!
//! As stored, the right-hand-rule normal of every triangle points towards
//! *decreasing* field values (into the solid).  Extractors that must emit
//! outward-facing triangles reverse the vertex order of each triple.
//!
//! The consistency of the tables with the cube topology is checked for all 256
//! cases by the unit tests at the bottom of this file: edge sets, closed oriented
//! patches bounded by cube faces, the winding relative to the field gradient and
//! face-by-face agreement of neighbouring cells (every checkerboard face
//! separates its inside corners, so cells never leave cracks between them).

/// Grid offsets `(di, dj, dk)` of the eight cube corners.
pub(crate) const MC_CORNER_OFFSETS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
];

/// The two corners joined by each of the twelve cube edges.
pub(crate) const MC_EDGE_CORNERS: [[usize; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];

/// Lattice description of each cube edge: the grid offset `(di, dj, dk)` of the
/// edge's lower end point and the axis (0 = x, 1 = y, 2 = z) the edge runs along.
///
/// Two neighbouring cells that share a physical edge see the same
/// `(lower end point, axis)` pair, so `3 * node_index(lower end point) + axis`
/// is a dense global edge id that welds the vertices of adjacent cells.
pub(crate) const MC_EDGE_ORIGIN_AXIS: [([usize; 3], usize); 12] = [
    ([0, 0, 0], 0),
    ([1, 0, 0], 1),
    ([0, 1, 0], 0),
    ([0, 0, 0], 1),
    ([0, 0, 1], 0),
    ([1, 0, 1], 1),
    ([0, 1, 1], 0),
    ([0, 0, 1], 1),
    ([0, 0, 0], 2),
    ([1, 0, 0], 2),
    ([1, 1, 0], 2),
    ([0, 1, 0], 2),
];

/// The corners of each cube edge ordered by increasing lattice coordinate along
/// the edge: `[lower end, upper end]`.
///
/// Neighbouring cells see a shared edge from opposite ends (edge 1 of one cell is
/// edge 3 of its +x neighbour, with the corners in the opposite order).  Always
/// interpolating from the lower end to the upper end makes the crossing point of
/// a shared edge bit-identical in every cell that owns it.
pub(crate) const MC_EDGE_ASCENDING: [[usize; 2]; 12] = ascending_edge_corners();

/// Builds [`MC_EDGE_ASCENDING`] from [`MC_EDGE_CORNERS`] and [`MC_EDGE_ORIGIN_AXIS`].
const fn ascending_edge_corners() -> [[usize; 2]; 12] {
    let mut ascending = [[0usize; 2]; 12];
    let mut edge = 0;
    while edge < 12 {
        let [a, b] = MC_EDGE_CORNERS[edge];
        let axis = MC_EDGE_ORIGIN_AXIS[edge].1;
        ascending[edge] = if MC_CORNER_OFFSETS[a][axis] == 0 {
            [a, b]
        } else {
            [b, a]
        };
        edge += 1;
    }
    ascending
}

/// Case index of a cube: bit `c` is set when corner `c` is *inside*, i.e. its
/// field value is strictly below `iso`.
pub(crate) fn case_index(corner_values: &[f64; 8], iso: f64) -> usize {
    corner_values
        .iter()
        .enumerate()
        .fold(0usize, |case, (corner, &value)| {
            if value < iso {
                case | (1 << corner)
            } else {
                case
            }
        })
}

/// Triangles of `case` as triples of cube-edge indices, in table order.
///
/// Yields nothing for the empty cases (0 and 255) and for an out-of-range case.
pub(crate) fn case_triangles(case: usize) -> impl Iterator<Item = [usize; 3]> {
    let row: &[i8] = MC_TRI_TABLE.get(case).map_or(&[][..], |r| &r[..]);
    row.chunks_exact(3).map_while(|triple| {
        Some([
            usize::try_from(triple[0]).ok()?,
            usize::try_from(triple[1]).ok()?,
            usize::try_from(triple[2]).ok()?,
        ])
    })
}

/// Bit mask of the cube edges crossed by the surface, indexed by case.
#[rustfmt::skip]
pub(crate) const MC_EDGE_TABLE: [u16; 256] = [
    0x000,0x109,0x203,0x30a,0x406,0x50f,0x605,0x70c,
    0x80c,0x905,0xa0f,0xb06,0xc0a,0xd03,0xe09,0xf00,
    0x190,0x099,0x393,0x29a,0x596,0x49f,0x795,0x69c,
    0x99c,0x895,0xb9f,0xa96,0xd9a,0xc93,0xf99,0xe90,
    0x230,0x339,0x033,0x13a,0x636,0x73f,0x435,0x53c,
    0xa3c,0xb35,0x83f,0x936,0xe3a,0xf33,0xc39,0xd30,
    0x3a0,0x2a9,0x1a3,0x0aa,0x7a6,0x6af,0x5a5,0x4ac,
    0xbac,0xaa5,0x9af,0x8a6,0xfaa,0xea3,0xda9,0xca0,
    0x460,0x569,0x663,0x76a,0x066,0x16f,0x265,0x36c,
    0xc6c,0xd65,0xe6f,0xf66,0x86a,0x963,0xa69,0xb60,
    0x5f0,0x4f9,0x7f3,0x6fa,0x1f6,0x0ff,0x3f5,0x2fc,
    0xdfc,0xcf5,0xfff,0xef6,0x9fa,0x8f3,0xbf9,0xaf0,
    0x650,0x759,0x453,0x55a,0x256,0x35f,0x055,0x15c,
    0xe5c,0xf55,0xc5f,0xd56,0xa5a,0xb53,0x859,0x950,
    0x7c0,0x6c9,0x5c3,0x4ca,0x3c6,0x2cf,0x1c5,0x0cc,
    0xfcc,0xec5,0xdcf,0xcc6,0xbca,0xac3,0x9c9,0x8c0,
    0x8c0,0x9c9,0xac3,0xbca,0xcc6,0xdcf,0xec5,0xfcc,
    0x0cc,0x1c5,0x2cf,0x3c6,0x4ca,0x5c3,0x6c9,0x7c0,
    0x950,0x859,0xb53,0xa5a,0xd56,0xc5f,0xf55,0xe5c,
    0x15c,0x055,0x35f,0x256,0x55a,0x453,0x759,0x650,
    0xaf0,0xbf9,0x8f3,0x9fa,0xef6,0xfff,0xcf5,0xdfc,
    0x2fc,0x3f5,0x0ff,0x1f6,0x6fa,0x7f3,0x4f9,0x5f0,
    0xb60,0xa69,0x963,0x86a,0xf66,0xe6f,0xd65,0xc6c,
    0x36c,0x265,0x16f,0x066,0x76a,0x663,0x569,0x460,
    0xca0,0xda9,0xea3,0xfaa,0x8a6,0x9af,0xaa5,0xbac,
    0x4ac,0x5a5,0x6af,0x7a6,0x0aa,0x1a3,0x2a9,0x3a0,
    0xd30,0xc39,0xf33,0xe3a,0x936,0x83f,0xb35,0xa3c,
    0x53c,0x435,0x73f,0x636,0x13a,0x033,0x339,0x230,
    0xe90,0xf99,0xc93,0xd9a,0xa96,0xb9f,0x895,0x99c,
    0x69c,0x795,0x49f,0x596,0x29a,0x393,0x099,0x190,
    0xf00,0xe09,0xd03,0xc0a,0xb06,0xa0f,0x905,0x80c,
    0x70c,0x605,0x50f,0x406,0x30a,0x203,0x109,0x000,
];

/// Triangles of every case as triples of edge indices, terminated by `-1`.
#[rustfmt::skip]
pub(crate) const MC_TRI_TABLE: [[i8; 16]; 256] = [
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,8,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,1,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,8,3,9,8,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,8,3,1,2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [9,2,10,0,2,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [2,8,3,2,10,8,10,9,8,-1,-1,-1,-1,-1,-1,-1],
    [3,11,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,11,2,8,11,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,9,0,2,3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,11,2,1,9,11,9,8,11,-1,-1,-1,-1,-1,-1,-1],
    [3,10,1,11,10,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,10,1,0,8,10,8,11,10,-1,-1,-1,-1,-1,-1,-1],
    [3,9,0,3,11,9,11,10,9,-1,-1,-1,-1,-1,-1,-1],
    [9,8,10,10,8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,3,0,7,3,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,1,9,8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,1,9,4,7,1,7,3,1,-1,-1,-1,-1,-1,-1,-1],
    [1,2,10,8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [3,4,7,3,0,4,1,2,10,-1,-1,-1,-1,-1,-1,-1],
    [9,2,10,9,0,2,8,4,7,-1,-1,-1,-1,-1,-1,-1],
    [2,10,9,2,9,7,2,7,3,7,9,4,-1,-1,-1,-1],
    [8,4,7,3,11,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11,4,7,11,2,4,2,0,4,-1,-1,-1,-1,-1,-1,-1],
    [9,0,1,8,4,7,2,3,11,-1,-1,-1,-1,-1,-1,-1],
    [4,7,11,9,4,11,9,11,2,9,2,1,-1,-1,-1,-1],
    [3,10,1,3,11,10,7,8,4,-1,-1,-1,-1,-1,-1,-1],
    [1,11,10,1,4,11,1,0,4,7,11,4,-1,-1,-1,-1],
    [4,7,8,9,0,11,9,11,10,11,0,3,-1,-1,-1,-1],
    [4,7,11,4,11,9,9,11,10,-1,-1,-1,-1,-1,-1,-1],
    [9,5,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [9,5,4,0,8,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,5,4,1,5,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [8,5,4,8,3,5,3,1,5,-1,-1,-1,-1,-1,-1,-1],
    [1,2,10,9,5,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [3,0,8,1,2,10,4,9,5,-1,-1,-1,-1,-1,-1,-1],
    [5,2,10,5,4,2,4,0,2,-1,-1,-1,-1,-1,-1,-1],
    [2,10,5,3,2,5,3,5,4,3,4,8,-1,-1,-1,-1],
    [9,5,4,2,3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,11,2,0,8,11,4,9,5,-1,-1,-1,-1,-1,-1,-1],
    [0,5,4,0,1,5,2,3,11,-1,-1,-1,-1,-1,-1,-1],
    [2,1,5,2,5,8,2,8,11,4,8,5,-1,-1,-1,-1],
    [10,3,11,10,1,3,9,5,4,-1,-1,-1,-1,-1,-1,-1],
    [4,9,5,0,8,1,8,10,1,8,11,10,-1,-1,-1,-1],
    [5,4,0,5,0,11,5,11,10,11,0,3,-1,-1,-1,-1],
    [5,4,8,5,8,10,10,8,11,-1,-1,-1,-1,-1,-1,-1],
    [9,7,8,5,7,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [9,3,0,9,5,3,5,7,3,-1,-1,-1,-1,-1,-1,-1],
    [0,7,8,0,1,7,1,5,7,-1,-1,-1,-1,-1,-1,-1],
    [1,5,3,3,5,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [9,7,8,9,5,7,10,1,2,-1,-1,-1,-1,-1,-1,-1],
    [10,1,2,9,5,0,5,3,0,5,7,3,-1,-1,-1,-1],
    [8,0,2,8,2,5,8,5,7,10,5,2,-1,-1,-1,-1],
    [2,10,5,2,5,3,3,5,7,-1,-1,-1,-1,-1,-1,-1],
    [7,9,5,7,8,9,3,11,2,-1,-1,-1,-1,-1,-1,-1],
    [9,5,7,9,7,2,9,2,0,2,7,11,-1,-1,-1,-1],
    [2,3,11,0,1,8,1,7,8,1,5,7,-1,-1,-1,-1],
    [11,2,1,11,1,7,7,1,5,-1,-1,-1,-1,-1,-1,-1],
    [9,5,8,8,5,7,10,1,3,10,3,11,-1,-1,-1,-1],
    [5,7,0,5,0,9,7,11,0,1,0,10,11,10,0,-1],
    [11,10,0,11,0,3,10,5,0,8,0,7,5,7,0,-1],
    [11,10,5,7,11,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [10,6,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,8,3,5,10,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [9,0,1,5,10,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,8,3,1,9,8,5,10,6,-1,-1,-1,-1,-1,-1,-1],
    [1,6,5,2,6,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,6,5,1,2,6,3,0,8,-1,-1,-1,-1,-1,-1,-1],
    [9,6,5,9,0,6,0,2,6,-1,-1,-1,-1,-1,-1,-1],
    [5,9,8,5,8,2,5,2,6,3,2,8,-1,-1,-1,-1],
    [2,3,11,10,6,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11,0,8,11,2,0,10,6,5,-1,-1,-1,-1,-1,-1,-1],
    [0,1,9,2,3,11,5,10,6,-1,-1,-1,-1,-1,-1,-1],
    [5,10,6,1,9,2,9,11,2,9,8,11,-1,-1,-1,-1],
    [6,3,11,6,5,3,5,1,3,-1,-1,-1,-1,-1,-1,-1],
    [0,8,11,0,11,5,0,5,1,5,11,6,-1,-1,-1,-1],
    [3,11,6,0,3,6,0,6,5,0,5,9,-1,-1,-1,-1],
    [6,5,9,6,9,11,11,9,8,-1,-1,-1,-1,-1,-1,-1],
    [5,10,6,4,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,3,0,4,7,3,6,5,10,-1,-1,-1,-1,-1,-1,-1],
    [1,9,0,5,10,6,8,4,7,-1,-1,-1,-1,-1,-1,-1],
    [10,6,5,1,9,7,1,7,3,7,9,4,-1,-1,-1,-1],
    [6,1,2,6,5,1,4,7,8,-1,-1,-1,-1,-1,-1,-1],
    [1,2,5,5,2,6,3,0,4,3,4,7,-1,-1,-1,-1],
    [8,4,7,9,0,5,0,6,5,0,2,6,-1,-1,-1,-1],
    [7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,-1],
    [3,11,2,7,8,4,10,6,5,-1,-1,-1,-1,-1,-1,-1],
    [5,10,6,4,7,2,4,2,0,2,7,11,-1,-1,-1,-1],
    [0,1,9,4,7,8,2,3,11,5,10,6,-1,-1,-1,-1],
    [9,2,1,9,11,2,9,4,11,7,11,4,5,10,6,-1],
    [8,4,7,3,11,5,3,5,1,5,11,6,-1,-1,-1,-1],
    [5,1,11,5,11,6,1,0,11,7,11,4,0,4,11,-1],
    [0,5,9,0,6,5,0,3,6,11,6,3,8,4,7,-1],
    [6,5,9,6,9,11,4,7,9,7,11,9,-1,-1,-1,-1],
    [10,4,9,6,4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,10,6,4,9,10,0,8,3,-1,-1,-1,-1,-1,-1,-1],
    [10,0,1,10,6,0,6,4,0,-1,-1,-1,-1,-1,-1,-1],
    [8,3,1,8,1,6,8,6,4,6,1,10,-1,-1,-1,-1],
    [1,4,9,1,2,4,2,6,4,-1,-1,-1,-1,-1,-1,-1],
    [3,0,8,1,2,9,2,4,9,2,6,4,-1,-1,-1,-1],
    [0,2,4,4,2,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [8,3,2,8,2,4,4,2,6,-1,-1,-1,-1,-1,-1,-1],
    [10,4,9,10,6,4,11,2,3,-1,-1,-1,-1,-1,-1,-1],
    [0,8,2,2,8,11,4,9,10,4,10,6,-1,-1,-1,-1],
    [3,11,2,0,1,6,0,6,4,6,1,10,-1,-1,-1,-1],
    [6,4,1,6,1,10,4,8,1,2,1,11,8,11,1,-1],
    [9,6,4,9,3,6,9,1,3,11,6,3,-1,-1,-1,-1],
    [8,11,1,8,1,0,11,6,1,9,1,4,6,4,1,-1],
    [3,11,6,3,6,0,0,6,4,-1,-1,-1,-1,-1,-1,-1],
    [6,4,8,11,6,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [7,10,6,7,8,10,8,9,10,-1,-1,-1,-1,-1,-1,-1],
    [0,7,3,0,10,7,0,9,10,6,7,10,-1,-1,-1,-1],
    [10,6,7,1,10,7,1,7,8,1,8,0,-1,-1,-1,-1],
    [10,6,7,10,7,1,1,7,3,-1,-1,-1,-1,-1,-1,-1],
    [1,2,6,1,6,8,1,8,9,8,6,7,-1,-1,-1,-1],
    [2,6,9,2,9,1,6,7,9,0,9,3,7,3,9,-1],
    [7,8,0,7,0,6,6,0,2,-1,-1,-1,-1,-1,-1,-1],
    [7,3,2,6,7,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [2,3,11,10,6,8,10,8,9,8,6,7,-1,-1,-1,-1],
    [2,0,7,2,7,11,0,9,7,6,7,10,9,10,7,-1],
    [1,8,0,1,7,8,1,10,7,6,7,10,2,3,11,-1],
    [11,2,1,11,1,7,10,6,1,6,7,1,-1,-1,-1,-1],
    [8,9,6,8,6,7,9,1,6,11,6,3,1,3,6,-1],
    [0,9,1,11,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [7,8,0,7,0,6,3,11,0,11,6,0,-1,-1,-1,-1],
    [7,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [7,6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [3,0,8,11,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,1,9,11,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [8,1,9,8,3,1,11,7,6,-1,-1,-1,-1,-1,-1,-1],
    [10,1,2,6,11,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,2,10,3,0,8,6,11,7,-1,-1,-1,-1,-1,-1,-1],
    [2,9,0,2,10,9,6,11,7,-1,-1,-1,-1,-1,-1,-1],
    [6,11,7,2,10,3,10,8,3,10,9,8,-1,-1,-1,-1],
    [7,2,3,6,2,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [7,0,8,7,6,0,6,2,0,-1,-1,-1,-1,-1,-1,-1],
    [2,7,6,2,3,7,0,1,9,-1,-1,-1,-1,-1,-1,-1],
    [1,6,2,1,8,6,1,9,8,8,7,6,-1,-1,-1,-1],
    [10,7,6,10,1,7,1,3,7,-1,-1,-1,-1,-1,-1,-1],
    [10,7,6,1,7,10,1,8,7,1,0,8,-1,-1,-1,-1],
    [0,3,7,0,7,10,0,10,9,6,10,7,-1,-1,-1,-1],
    [7,6,10,7,10,8,8,10,9,-1,-1,-1,-1,-1,-1,-1],
    [6,8,4,11,8,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [3,6,11,3,0,6,0,4,6,-1,-1,-1,-1,-1,-1,-1],
    [8,6,11,8,4,6,9,0,1,-1,-1,-1,-1,-1,-1,-1],
    [9,4,6,9,6,3,9,3,1,11,3,6,-1,-1,-1,-1],
    [6,8,4,6,11,8,2,10,1,-1,-1,-1,-1,-1,-1,-1],
    [1,2,10,3,0,11,0,6,11,0,4,6,-1,-1,-1,-1],
    [4,11,8,4,6,11,0,2,9,2,10,9,-1,-1,-1,-1],
    [10,9,3,10,3,2,9,4,3,11,3,6,4,6,3,-1],
    [8,2,3,8,4,2,4,6,2,-1,-1,-1,-1,-1,-1,-1],
    [0,4,2,4,6,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,9,0,2,3,4,2,4,6,4,3,8,-1,-1,-1,-1],
    [1,9,4,1,4,2,2,4,6,-1,-1,-1,-1,-1,-1,-1],
    [8,1,3,8,6,1,8,4,6,6,10,1,-1,-1,-1,-1],
    [10,1,0,10,0,6,6,0,4,-1,-1,-1,-1,-1,-1,-1],
    [4,6,3,4,3,8,6,10,3,0,3,9,10,9,3,-1],
    [10,9,4,6,10,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,9,5,7,6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,8,3,4,9,5,11,7,6,-1,-1,-1,-1,-1,-1,-1],
    [5,0,1,5,4,0,7,6,11,-1,-1,-1,-1,-1,-1,-1],
    [11,7,6,8,3,4,3,5,4,3,1,5,-1,-1,-1,-1],
    [9,5,4,10,1,2,7,6,11,-1,-1,-1,-1,-1,-1,-1],
    [6,11,7,1,2,10,0,8,3,4,9,5,-1,-1,-1,-1],
    [7,6,11,5,4,10,4,2,10,4,0,2,-1,-1,-1,-1],
    [3,4,8,3,5,4,3,2,5,10,5,2,11,7,6,-1],
    [7,2,3,7,6,2,5,4,9,-1,-1,-1,-1,-1,-1,-1],
    [9,5,4,0,8,6,0,6,2,6,8,7,-1,-1,-1,-1],
    [3,6,2,3,7,6,1,5,0,5,4,0,-1,-1,-1,-1],
    [6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,-1],
    [9,5,4,10,1,6,1,7,6,1,3,7,-1,-1,-1,-1],
    [1,6,10,1,7,6,1,0,7,8,7,0,9,5,4,-1],
    [4,0,10,4,10,5,0,3,10,6,10,7,3,7,10,-1],
    [7,6,10,7,10,8,5,4,10,4,8,10,-1,-1,-1,-1],
    [6,9,5,6,11,9,11,8,9,-1,-1,-1,-1,-1,-1,-1],
    [3,6,11,0,6,3,0,5,6,0,9,5,-1,-1,-1,-1],
    [0,11,8,0,5,11,0,1,5,5,6,11,-1,-1,-1,-1],
    [6,11,3,6,3,5,5,3,1,-1,-1,-1,-1,-1,-1,-1],
    [1,2,10,9,5,11,9,11,8,11,5,6,-1,-1,-1,-1],
    [0,11,3,0,6,11,0,9,6,5,6,9,1,2,10,-1],
    [11,8,5,11,5,6,8,0,5,10,5,2,0,2,5,-1],
    [6,11,3,6,3,5,2,10,3,10,5,3,-1,-1,-1,-1],
    [5,8,9,5,2,8,5,6,2,3,8,2,-1,-1,-1,-1],
    [9,5,6,9,6,0,0,6,2,-1,-1,-1,-1,-1,-1,-1],
    [1,5,8,1,8,0,5,6,8,3,8,2,6,2,8,-1],
    [1,5,6,2,1,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,3,6,1,6,10,3,8,6,5,6,9,8,9,6,-1],
    [10,1,0,10,0,6,9,5,0,5,6,0,-1,-1,-1,-1],
    [0,3,8,5,6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [10,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11,5,10,7,5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11,5,10,11,7,5,8,3,0,-1,-1,-1,-1,-1,-1,-1],
    [5,11,7,5,10,11,1,9,0,-1,-1,-1,-1,-1,-1,-1],
    [10,7,5,10,11,7,9,8,1,8,3,1,-1,-1,-1,-1],
    [11,1,2,11,7,1,7,5,1,-1,-1,-1,-1,-1,-1,-1],
    [0,8,3,1,2,7,1,7,5,7,2,11,-1,-1,-1,-1],
    [9,7,5,9,2,7,9,0,2,2,11,7,-1,-1,-1,-1],
    [7,5,2,7,2,11,5,9,2,3,2,8,9,8,2,-1],
    [2,5,10,2,3,5,3,7,5,-1,-1,-1,-1,-1,-1,-1],
    [8,2,0,8,5,2,8,7,5,10,2,5,-1,-1,-1,-1],
    [9,0,1,5,10,3,5,3,7,3,10,2,-1,-1,-1,-1],
    [9,8,2,9,2,1,8,7,2,10,2,5,7,5,2,-1],
    [1,3,5,3,7,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,8,7,0,7,1,1,7,5,-1,-1,-1,-1,-1,-1,-1],
    [9,0,3,9,3,5,5,3,7,-1,-1,-1,-1,-1,-1,-1],
    [9,8,7,5,9,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [5,8,4,5,10,8,10,11,8,-1,-1,-1,-1,-1,-1,-1],
    [5,0,4,5,11,0,5,10,11,11,3,0,-1,-1,-1,-1],
    [0,1,9,8,4,10,8,10,11,10,4,5,-1,-1,-1,-1],
    [10,11,4,10,4,5,11,3,4,9,4,1,3,1,4,-1],
    [2,5,1,2,8,5,2,11,8,4,5,8,-1,-1,-1,-1],
    [0,4,11,0,11,3,4,5,11,2,11,1,5,1,11,-1],
    [0,2,5,0,5,9,2,11,5,4,5,8,11,8,5,-1],
    [9,4,5,2,11,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [2,5,10,3,5,2,3,4,5,3,8,4,-1,-1,-1,-1],
    [5,10,2,5,2,4,4,2,0,-1,-1,-1,-1,-1,-1,-1],
    [3,10,2,3,5,10,3,8,5,4,5,8,0,1,9,-1],
    [5,10,2,5,2,4,1,9,2,9,4,2,-1,-1,-1,-1],
    [8,4,5,8,5,3,3,5,1,-1,-1,-1,-1,-1,-1,-1],
    [0,4,5,1,0,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [8,4,5,8,5,3,9,0,5,0,3,5,-1,-1,-1,-1],
    [9,4,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,11,7,4,9,11,9,10,11,-1,-1,-1,-1,-1,-1,-1],
    [0,8,3,4,9,7,9,11,7,9,10,11,-1,-1,-1,-1],
    [1,10,11,1,11,4,1,4,0,7,4,11,-1,-1,-1,-1],
    [3,1,4,3,4,8,1,10,4,7,4,11,10,11,4,-1],
    [4,11,7,9,11,4,9,2,11,9,1,2,-1,-1,-1,-1],
    [9,7,4,9,11,7,9,1,11,2,11,1,0,8,3,-1],
    [11,7,4,11,4,2,2,4,0,-1,-1,-1,-1,-1,-1,-1],
    [11,7,4,11,4,2,8,3,4,3,2,4,-1,-1,-1,-1],
    [2,9,10,2,7,9,2,3,7,7,4,9,-1,-1,-1,-1],
    [9,10,7,9,7,4,10,2,7,8,7,0,2,0,7,-1],
    [3,7,10,3,10,2,7,4,10,1,10,0,4,0,10,-1],
    [1,10,2,8,7,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,9,1,4,1,7,7,1,3,-1,-1,-1,-1,-1,-1,-1],
    [4,9,1,4,1,7,0,8,1,8,7,1,-1,-1,-1,-1],
    [4,0,3,7,4,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [9,10,8,10,11,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [3,0,9,3,9,11,11,9,10,-1,-1,-1,-1,-1,-1,-1],
    [0,1,10,0,10,8,8,10,11,-1,-1,-1,-1,-1,-1,-1],
    [3,1,10,11,3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,2,11,1,11,9,9,11,8,-1,-1,-1,-1,-1,-1,-1],
    [3,0,9,3,9,11,1,2,9,2,11,9,-1,-1,-1,-1],
    [0,2,11,8,0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [3,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [2,3,8,2,8,10,10,8,9,-1,-1,-1,-1,-1,-1,-1],
    [9,10,2,0,9,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [2,3,8,2,8,10,0,1,8,1,10,8,-1,-1,-1,-1],
    [1,10,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [1,3,8,9,1,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,9,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [0,3,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
];

/// Mesh analysis helpers shared by the isosurface-extraction regression tests.
#[cfg(test)]
pub(crate) mod test_support {
    use std::collections::HashMap;

    /// Welds a triangle soup by the exact bit pattern of its coordinates.
    pub(crate) fn weld_soup(soup: &[[[f64; 3]; 3]]) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let mut ids: HashMap<[u64; 3], usize> = HashMap::new();
        let mut vertices: Vec<[f64; 3]> = Vec::new();
        let mut triangles: Vec<[usize; 3]> = Vec::with_capacity(soup.len());
        for tri in soup {
            let mut idx = [0usize; 3];
            for (k, p) in tri.iter().enumerate() {
                // `+ 0.0` maps -0.0 to +0.0 so both zeros weld together.
                let key = [
                    (p[0] + 0.0).to_bits(),
                    (p[1] + 0.0).to_bits(),
                    (p[2] + 0.0).to_bits(),
                ];
                idx[k] = *ids.entry(key).or_insert_with(|| {
                    vertices.push(*p);
                    vertices.len() - 1
                });
            }
            triangles.push(idx);
        }
        (vertices, triangles)
    }

    /// Number of triangles incident to every undirected edge, as
    /// `(valence, number of edges with that valence)` in ascending valence order.
    pub(crate) fn edge_valence_histogram(triangles: &[[usize; 3]]) -> Vec<(usize, usize)> {
        let mut valence: HashMap<(usize, usize), usize> = HashMap::new();
        for t in triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *valence.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let mut histogram: HashMap<usize, usize> = HashMap::new();
        for &v in valence.values() {
            *histogram.entry(v).or_insert(0) += 1;
        }
        let mut out: Vec<(usize, usize)> = histogram.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// `true` when every directed edge occurs exactly once and so does its reverse,
    /// i.e. the triangles form a closed, consistently oriented surface.
    pub(crate) fn is_closed_and_consistently_oriented(triangles: &[[usize; 3]]) -> bool {
        let mut directed: HashMap<(usize, usize), usize> = HashMap::new();
        for t in triangles {
            for k in 0..3 {
                *directed.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
            }
        }
        directed
            .iter()
            .all(|(&(a, b), &count)| count == 1 && directed.get(&(b, a)) == Some(&1))
    }

    /// Total surface area and signed volume (divergence theorem) of a mesh.  The
    /// volume is positive when the triangles face outwards.
    pub(crate) fn area_and_volume(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> (f64, f64) {
        let mut area = 0.0;
        let mut volume = 0.0;
        for t in triangles {
            let (a, b, c) = (vertices[t[0]], vertices[t[1]], vertices[t[2]]);
            let (u, v) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
            );
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            area += 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let bc = [
                b[1] * c[2] - b[2] * c[1],
                b[2] * c[0] - b[0] * c[2],
                b[0] * c[1] - b[1] * c[0],
            ];
            volume += (a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2]) / 6.0;
        }
        (area, volume)
    }

    /// Number of triangles with (numerically) zero area.
    pub(crate) fn degenerate_triangle_count(
        vertices: &[[f64; 3]],
        triangles: &[[usize; 3]],
    ) -> usize {
        triangles
            .iter()
            .filter(|t| {
                let (a, b, c) = (vertices[t[0]], vertices[t[1]], vertices[t[2]]);
                let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                let n = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                n[0] * n[0] + n[1] * n[1] + n[2] * n[2] < 1e-24
            })
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cube face: `(axis, side)` with `side` in `{0, 1}`.
    type Face = (usize, usize);

    /// The six faces of the cube.
    const ALL_FACES: [Face; 6] = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)];

    /// Edge mask derived from the cube topology alone: edge `e` is crossed when
    /// its two corners lie on different sides of the surface.
    fn derived_edge_mask(case: usize) -> u16 {
        let mut mask = 0u16;
        for (e, [a, b]) in MC_EDGE_CORNERS.iter().enumerate() {
            if ((case >> a) & 1) != ((case >> b) & 1) {
                mask |= 1 << e;
            }
        }
        mask
    }

    /// The triangles of `case` parsed strictly from the raw table row.
    fn triangles_of(case: usize) -> Vec<[usize; 3]> {
        let row = &MC_TRI_TABLE[case];
        let used = row.iter().take_while(|&&e| e != -1).count();
        assert_eq!(used % 3, 0, "case {case}: edge list is not a multiple of 3");
        assert!(
            row[used..].iter().all(|&e| e == -1),
            "case {case}: entries after the terminator must all be -1"
        );
        row[..used]
            .chunks(3)
            .map(|t| {
                let idx = |n: usize| usize::try_from(t[n]).expect("edge index is non-negative");
                [idx(0), idx(1), idx(2)]
            })
            .collect()
    }

    /// The two cube faces that contain edge `e`.
    fn faces_of_edge(e: usize) -> Vec<Face> {
        let [a, b] = MC_EDGE_CORNERS[e];
        ALL_FACES
            .into_iter()
            .filter(|&(axis, side)| {
                MC_CORNER_OFFSETS[a][axis] == side && MC_CORNER_OFFSETS[b][axis] == side
            })
            .collect()
    }

    /// The face containing both edges, if they share one.
    fn common_face(e1: usize, e2: usize) -> Option<Face> {
        let f2 = faces_of_edge(e2);
        faces_of_edge(e1).into_iter().find(|f| f2.contains(f))
    }

    #[test]
    fn every_edge_lies_on_exactly_two_faces() {
        for e in 0..12 {
            assert_eq!(faces_of_edge(e).len(), 2, "edge {e}");
        }
    }

    #[test]
    fn edge_table_matches_cube_topology_for_every_case() {
        for (case, &mask) in MC_EDGE_TABLE.iter().enumerate() {
            assert_eq!(
                mask,
                derived_edge_mask(case),
                "edge mask of case {case} (0x{case:02x}) disagrees with the cube topology"
            );
        }
    }

    #[test]
    fn edge_table_entries_213_to_215_are_the_standard_values() {
        assert_eq!(MC_EDGE_TABLE[213], 0x83f);
        assert_eq!(MC_EDGE_TABLE[214], 0xb35);
        assert_eq!(MC_EDGE_TABLE[215], 0xa3c);
    }

    #[test]
    fn edge_table_is_symmetric_under_complement() {
        for case in 0..256 {
            assert_eq!(
                MC_EDGE_TABLE[case],
                MC_EDGE_TABLE[255 - case],
                "case {case} and its complement cross different edges"
            );
        }
    }

    #[test]
    fn triangle_table_uses_exactly_the_crossed_edges() {
        for (case, &mask) in MC_EDGE_TABLE.iter().enumerate() {
            let mut used = 0u16;
            for tri in triangles_of(case) {
                assert!(
                    tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2],
                    "case {case}: degenerate triangle {tri:?}"
                );
                for e in tri {
                    assert!(e < 12, "case {case}: edge index {e} out of range");
                    used |= 1 << e;
                }
            }
            assert_eq!(
                used, mask,
                "case {case}: triangles reference a different edge set than the edge table"
            );
        }
    }

    #[test]
    fn triangle_counts_follow_the_known_case_structure() {
        assert!(triangles_of(0).is_empty());
        assert!(triangles_of(255).is_empty());
        for corner in 0..8 {
            assert_eq!(
                triangles_of(1 << corner).len(),
                1,
                "one inside corner cuts exactly one triangle"
            );
            assert_eq!(triangles_of(255 ^ (1 << corner)).len(), 1);
        }
        let max = (0..256).map(|c| triangles_of(c).len()).max();
        assert_eq!(max, Some(5), "the largest standard case has five triangles");
    }

    /// Every triangle side is either shared by exactly two triangles with
    /// opposite directions (an interior side of the patch) or lies on one cube
    /// face (a boundary side); on each face the boundary sides pair up exactly
    /// the crossed edges of that face.
    #[test]
    fn triangle_table_is_a_consistently_oriented_patch_bounded_by_cube_faces() {
        for (case, &edge_mask) in MC_EDGE_TABLE.iter().enumerate() {
            let mut sides: Vec<(usize, usize)> = Vec::new();
            for tri in triangles_of(case) {
                for s in 0..3 {
                    sides.push((tri[s], tri[(s + 1) % 3]));
                }
            }
            let count = |a: usize, b: usize| sides.iter().filter(|&&s| s == (a, b)).count();

            let mut boundary: Vec<(usize, usize, Face)> = Vec::new();
            for &(a, b) in &sides {
                let forward = count(a, b);
                let reverse = count(b, a);
                if reverse > 0 {
                    assert_eq!(
                        (forward, reverse),
                        (1, 1),
                        "case {case}: side {a}->{b} is not shared by exactly two opposite triangles"
                    );
                } else {
                    assert_eq!(forward, 1, "case {case}: boundary side {a}->{b} repeated");
                    let face = common_face(a, b).unwrap_or_else(|| {
                        panic!("case {case}: boundary side {a}->{b} does not lie on a cube face")
                    });
                    boundary.push((a, b, face));
                }
            }

            for face in ALL_FACES {
                let (axis, side) = face;
                let crossed: Vec<usize> = (0..12)
                    .filter(|&e| edge_mask & (1 << e) != 0 && faces_of_edge(e).contains(&face))
                    .collect();
                for &e in &crossed {
                    let touching = boundary
                        .iter()
                        .filter(|&&(a, b, f)| f == face && (a == e || b == e))
                        .count();
                    assert_eq!(
                        touching, 1,
                        "case {case}: crossed edge {e} must meet exactly one boundary side on face {face:?}"
                    );
                }
                let on_face = boundary.iter().filter(|&&(_, _, f)| f == face).count();
                assert_eq!(
                    on_face * 2,
                    crossed.len(),
                    "case {case}: boundary sides on face {face:?} must pair up its crossed edges"
                );

                // A checkerboard face (two diagonal inside corners) is ambiguous.
                // Neighbouring cells only agree on it when every case resolves it
                // the same way: each inside corner is cut off on its own.
                let inside: Vec<usize> = (0..8)
                    .filter(|&c| MC_CORNER_OFFSETS[c][axis] == side && (case >> c) & 1 == 1)
                    .collect();
                let joined = MC_EDGE_CORNERS
                    .iter()
                    .any(|&[a, b]| inside.contains(&a) && inside.contains(&b));
                if inside.len() == 2 && !joined {
                    for &corner in &inside {
                        let at_corner: Vec<usize> = crossed
                            .iter()
                            .copied()
                            .filter(|&e| MC_EDGE_CORNERS[e].contains(&corner))
                            .collect();
                        assert_eq!(at_corner.len(), 2, "case {case}: corner {corner}");
                        let cut_off = boundary.iter().any(|&(a, b, f)| {
                            f == face
                                && ((a == at_corner[0] && b == at_corner[1])
                                    || (a == at_corner[1] && b == at_corner[0]))
                        });
                        assert!(
                            cut_off,
                            "case {case}: checkerboard face {face:?} must cut off inside corner {corner}"
                        );
                    }
                }
            }
        }
    }

    /// Unit-cube geometry check: with field value `-1` at inside corners and `+1`
    /// at outside corners every crossing is the edge midpoint, and the sum of the
    /// triangle area vectors (right-hand rule on the stored vertex order) must
    /// point towards decreasing field value, i.e. into the solid.
    #[test]
    fn triangle_table_winds_towards_decreasing_field_for_every_case() {
        let corner_pos = |c: usize| -> [f64; 3] {
            let o = MC_CORNER_OFFSETS[c];
            [o[0] as f64, o[1] as f64, o[2] as f64]
        };
        let midpoint = |e: usize| -> [f64; 3] {
            let [a, b] = MC_EDGE_CORNERS[e];
            let (pa, pb) = (corner_pos(a), corner_pos(b));
            [
                0.5 * (pa[0] + pb[0]),
                0.5 * (pa[1] + pb[1]),
                0.5 * (pa[2] + pb[2]),
            ]
        };
        for case in 1..255 {
            let value = |c: usize| if (case >> c) & 1 == 1 { -1.0 } else { 1.0 };
            let mut weighted = 0.0;
            for tri in triangles_of(case) {
                let (p0, p1, p2) = (midpoint(tri[0]), midpoint(tri[1]), midpoint(tri[2]));
                let (u, v) = (
                    [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]],
                    [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]],
                );
                let normal = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                // Trilinear field gradient at the triangle centroid.
                let g = [
                    (p0[0] + p1[0] + p2[0]) / 3.0,
                    (p0[1] + p1[1] + p2[1]) / 3.0,
                    (p0[2] + p1[2] + p2[2]) / 3.0,
                ];
                let mut grad = [0.0; 3];
                for (c, &o) in MC_CORNER_OFFSETS.iter().enumerate() {
                    let w = |axis: usize| if o[axis] == 1 { g[axis] } else { 1.0 - g[axis] };
                    let dw = |axis: usize| if o[axis] == 1 { 1.0 } else { -1.0 };
                    grad[0] += value(c) * dw(0) * w(1) * w(2);
                    grad[1] += value(c) * w(0) * dw(1) * w(2);
                    grad[2] += value(c) * w(0) * w(1) * dw(2);
                }
                weighted += normal[0] * grad[0] + normal[1] * grad[1] + normal[2] * grad[2];
            }
            assert!(
                weighted < 0.0,
                "case {case}: triangle normals do not point towards decreasing field (sum {weighted})"
            );
        }
    }

    #[test]
    fn edge_origin_axis_matches_edge_corners() {
        for (e, &([oi, oj, ok], axis)) in MC_EDGE_ORIGIN_AXIS.iter().enumerate() {
            let [a, b] = MC_EDGE_CORNERS[e];
            let (pa, pb) = (MC_CORNER_OFFSETS[a], MC_CORNER_OFFSETS[b]);
            let differing: Vec<usize> = (0..3).filter(|&k| pa[k] != pb[k]).collect();
            assert_eq!(differing, vec![axis], "edge {e} must run along axis {axis}");
            let lower = [pa[0].min(pb[0]), pa[1].min(pb[1]), pa[2].min(pb[2])];
            assert_eq!([oi, oj, ok], lower, "edge {e} lower end point");
        }
    }

    #[test]
    fn ascending_corners_run_from_the_lower_to_the_upper_lattice_end() {
        for (e, &[lo, hi]) in MC_EDGE_ASCENDING.iter().enumerate() {
            let axis = MC_EDGE_ORIGIN_AXIS[e].1;
            assert_eq!(MC_CORNER_OFFSETS[lo][axis], 0, "edge {e} lower end");
            assert_eq!(MC_CORNER_OFFSETS[hi][axis], 1, "edge {e} upper end");
            let mut ends = MC_EDGE_CORNERS[e];
            ends.sort_unstable();
            let mut sorted = [lo, hi];
            sorted.sort_unstable();
            assert_eq!(ends, sorted, "edge {e} must join the same two corners");
        }
    }

    #[test]
    fn dense_edge_ids_are_shared_between_neighbouring_cells() {
        // Cell (0,0,0) edge 1 (x = 1 face, running along y) is the same physical
        // edge as cell (1,0,0) edge 3 (x = 0 face, running along y).
        let id = |cell: [usize; 3], edge: usize| {
            let (off, axis) = MC_EDGE_ORIGIN_AXIS[edge];
            let node = [cell[0] + off[0], cell[1] + off[1], cell[2] + off[2]];
            (node, axis)
        };
        assert_eq!(id([0, 0, 0], 1), id([1, 0, 0], 3));
        assert_eq!(id([0, 0, 0], 5), id([1, 0, 0], 7));
        assert_eq!(id([0, 0, 0], 9), id([1, 0, 0], 8));
        assert_eq!(id([0, 0, 0], 2), id([0, 1, 0], 0));
        assert_eq!(id([0, 0, 0], 10), id([1, 1, 0], 8));
    }

    #[test]
    fn case_index_marks_strictly_inside_corners() {
        let values = [-1.0, 1.0, -0.5, 0.0, 2.0, -3.0, 0.25, 0.0];
        // Corners 0, 2 and 5 are strictly below zero; value == iso is outside.
        assert_eq!(case_index(&values, 0.0), 0b0010_0101);
        assert_eq!(case_index(&[1.0; 8], 0.0), 0);
        assert_eq!(case_index(&[-1.0; 8], 0.0), 255);
        assert_eq!(case_index(&values, 0.3), 0b1110_1101);
    }

    #[test]
    fn case_triangles_agrees_with_the_raw_table_rows() {
        for case in 0..256 {
            let from_iter: Vec<[usize; 3]> = case_triangles(case).collect();
            assert_eq!(from_iter, triangles_of(case), "case {case}");
        }
        assert_eq!(case_triangles(256).count(), 0);
    }
}
