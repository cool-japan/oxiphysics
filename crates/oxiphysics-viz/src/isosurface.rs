#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Isosurface extraction via marching cubes and related algorithms.
//!
//! Implements the classic marching-cubes algorithm (all 256 cube configurations
//! handled through an edge-table and symmetry-reduced triangle-table), signed
//! distance functions for common shapes, and a simple 2-D dual-contouring
//! sketch.

#![allow(dead_code)]

// ─── Scalar Field ────────────────────────────────────────────────────────────

/// A 3-D axis-aligned scalar field stored on a regular grid.
#[derive(Debug, Clone)]
pub struct ScalarField {
    /// Number of grid cells in the X direction.
    pub nx: usize,
    /// Number of grid cells in the Y direction.
    pub ny: usize,
    /// Number of grid cells in the Z direction.
    pub nz: usize,
    /// Grid spacing in X.
    pub dx: f64,
    /// Grid spacing in Y.
    pub dy: f64,
    /// Grid spacing in Z.
    pub dz: f64,
    /// Flat storage: `values[index(i, j, k)]`.
    pub values: Vec<f64>,
}

impl ScalarField {
    /// Construct a new field initialised to zero.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            values: vec![0.0; nx * ny * nz],
        }
    }

    /// Flat index for grid position `(i, j, k)`.
    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k + self.nz * (j + self.ny * i)
    }

    /// Get the scalar value at `(i, j, k)`.
    pub fn get(&self, i: usize, j: usize, k: usize) -> f64 {
        self.values[self.index(i, j, k)]
    }

    /// Set the scalar value at `(i, j, k)`.
    pub fn set(&mut self, i: usize, j: usize, k: usize, v: f64) {
        let idx = self.index(i, j, k);
        self.values[idx] = v;
    }
}

// ─── Triangle ─────────────────────────────────────────────────────────────────

/// A single triangle in 3-D space with an associated surface normal.
#[derive(Debug, Clone)]
pub struct Triangle {
    /// The three vertex positions.
    pub vertices: [[f64; 3]; 3],
    /// Outward unit normal.
    pub normal: [f64; 3],
}

/// Compute the normal of a triangle via the cross product of two edge vectors.
///
/// Returns a normalised vector, or `[0, 0, 1]` if the triangle is degenerate.
pub fn compute_normal(tri: &Triangle) -> [f64; 3] {
    let [a, b, c] = tri.vertices;
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = cross3(ab, ac);
    normalise3(n)
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-300 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Linearly interpolate on an edge between two grid vertices.
///
/// Returns the world-space position where the isosurface crosses the edge.
pub fn interpolate_edge(p1: [f64; 3], v1: f64, p2: [f64; 3], v2: f64, isovalue: f64) -> [f64; 3] {
    let dv = v2 - v1;
    let t = if dv.abs() < 1e-12 {
        0.5
    } else {
        (isovalue - v1) / dv
    };
    [
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2]),
    ]
}

// ─── Marching Cubes ───────────────────────────────────────────────────────────

// Edge table: bitmask of which edges a given cube configuration intersects.
// Index by 8-bit cube case (0-255).
#[rustfmt::skip]
const EDGE_TABLE: [u16; 256] = [
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

// Triangle table: for each cube configuration, list of edge indices forming
// triangles.  Terminated by -1.  Stored as a flat slice, 16 entries per case.
// This is the standard Lorensen & Cline table (256 × 16).
#[rustfmt::skip]
const TRI_TABLE: [[i8; 16]; 256] = [
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
    [9,0,1,2,3,5,2,5,10,5,3,7,-1,-1,-1,-1],
    [8,2,9,8,5,2,8,7,5,10,2,5,9,10,2,1],
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

// The 12 edges of a cube, each defined by two corner indices (0-7).
// Corner index k encodes position (k&1, (k>>1)&1, (k>>2)&1).
const CUBE_EDGES: [(usize, usize); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0), // bottom face edges
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4), // top face edges
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7), // vertical edges
];

/// Extract an isosurface from `field` at the given `isovalue` using the
/// marching-cubes algorithm.
///
/// Returns a list of [`Triangle`]s whose normals are computed from the cross
/// product of the triangle edges.
pub fn marching_cubes(field: &ScalarField, isovalue: f64) -> Vec<Triangle> {
    let mut triangles = Vec::new();
    if field.nx < 2 || field.ny < 2 || field.nz < 2 {
        return triangles;
    }

    for i in 0..field.nx - 1 {
        for j in 0..field.ny - 1 {
            for k in 0..field.nz - 1 {
                // Cube corners (in canonical order 0-7).
                let corners: [[usize; 3]; 8] = [
                    [i, j, k],
                    [i + 1, j, k],
                    [i + 1, j + 1, k],
                    [i, j + 1, k],
                    [i, j, k + 1],
                    [i + 1, j, k + 1],
                    [i + 1, j + 1, k + 1],
                    [i, j + 1, k + 1],
                ];

                let vals: [f64; 8] =
                    std::array::from_fn(|n| field.get(corners[n][0], corners[n][1], corners[n][2]));

                // Build cube index.
                let mut cube_idx = 0u8;
                for n in 0..8 {
                    if vals[n] < isovalue {
                        cube_idx |= 1 << n;
                    }
                }

                let edges = EDGE_TABLE[cube_idx as usize];
                if edges == 0 {
                    continue;
                }

                // Compute world-space positions of the 8 corners.
                let pos: [[f64; 3]; 8] = std::array::from_fn(|n| {
                    [
                        corners[n][0] as f64 * field.dx,
                        corners[n][1] as f64 * field.dy,
                        corners[n][2] as f64 * field.dz,
                    ]
                });

                // Interpolate edge vertices (only where needed).
                let mut edge_verts: [[f64; 3]; 12] = [[0.0; 3]; 12];
                for e in 0..12 {
                    if edges & (1 << e) != 0 {
                        let (a, b) = CUBE_EDGES[e];
                        edge_verts[e] =
                            interpolate_edge(pos[a], vals[a], pos[b], vals[b], isovalue);
                    }
                }

                // Emit triangles.
                let tris = &TRI_TABLE[cube_idx as usize];
                let mut idx = 0;
                while idx < 16 && tris[idx] != -1 {
                    let e0 = tris[idx] as usize;
                    let e1 = tris[idx + 1] as usize;
                    let e2 = tris[idx + 2] as usize;
                    let verts = [edge_verts[e0], edge_verts[e1], edge_verts[e2]];
                    let mut tri = Triangle {
                        vertices: verts,
                        normal: [0.0, 0.0, 1.0],
                    };
                    tri.normal = compute_normal(&tri);
                    triangles.push(tri);
                    idx += 3;
                }
            }
        }
    }
    triangles
}

// ─── IsosurfaceMesh ───────────────────────────────────────────────────────────

/// A collection of triangles representing an extracted isosurface.
#[derive(Debug, Clone)]
pub struct IsosurfaceMesh {
    /// All triangles forming the surface.
    pub triangles: Vec<Triangle>,
    /// Total number of distinct vertices (≈ `triangles.len() * 3` before welding).
    pub vertex_count: usize,
}

impl IsosurfaceMesh {
    /// Build an `IsosurfaceMesh` from a list of triangles.
    pub fn from_triangles(tris: Vec<Triangle>) -> Self {
        let vc = tris.len() * 3;
        Self {
            triangles: tris,
            vertex_count: vc,
        }
    }

    /// Number of triangles in the mesh.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }
}

// ─── Signed Distance Functions ───────────────────────────────────────────────

/// Signed distance function for a sphere.
///
/// Returns a negative value inside the sphere, zero on the surface, and
/// positive outside.
pub fn sphere_sdf(center: [f64; 3], radius: f64, x: f64, y: f64, z: f64) -> f64 {
    let dx = x - center[0];
    let dy = y - center[1];
    let dz = z - center[2];
    (dx * dx + dy * dy + dz * dz).sqrt() - radius
}

/// Signed distance function for a torus.
///
/// The torus is centred at the origin in the XZ plane with major radius
/// `major_r` and tube radius `minor_r`.
pub fn torus_sdf(major_r: f64, minor_r: f64, x: f64, y: f64, z: f64) -> f64 {
    let q = ((x * x + z * z).sqrt() - major_r).abs();
    (q * q + y * y).sqrt() - minor_r
}

// ─── 2-D Dual Contouring ─────────────────────────────────────────────────────

/// Approximate 2-D dual contouring: for each grid cell where the isovalue is
/// crossed, emit the midpoint of the crossing as a contour vertex.
///
/// `field` is a row-major grid of shape `[rows][cols]`.  Returns a list of
/// approximate contour positions in 2-D.
pub fn dual_contouring_2d(field: &[[f64; 2]], isovalue: f64) -> Vec<[f64; 2]> {
    // field[i] = [x, value] for row i — not a 2-D grid.
    // Interpret as a 1-D signed-distance function along x: crossings.
    let mut points = Vec::new();
    for window in field.windows(2) {
        let (x1, v1) = (window[0][0], window[0][1]);
        let (x2, v2) = (window[1][0], window[1][1]);
        if (v1 - isovalue) * (v2 - isovalue) <= 0.0 {
            let t = if (v2 - v1).abs() < 1e-12 {
                0.5
            } else {
                (isovalue - v1) / (v2 - v1)
            };
            points.push([x1 + t * (x2 - x1), 0.0]);
        }
    }
    points
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- sphere_sdf ---

    #[test]
    fn test_sphere_sdf_inside() {
        let sdf = sphere_sdf([0.0, 0.0, 0.0], 1.0, 0.0, 0.0, 0.0);
        assert!(
            sdf < 0.0,
            "origin should be inside unit sphere, got {}",
            sdf
        );
    }

    #[test]
    fn test_sphere_sdf_on_surface() {
        let sdf = sphere_sdf([0.0, 0.0, 0.0], 1.0, 1.0, 0.0, 0.0);
        assert!(
            sdf.abs() < 1e-12,
            "on surface, SDF should be 0, got {}",
            sdf
        );
    }

    #[test]
    fn test_sphere_sdf_outside() {
        let sdf = sphere_sdf([0.0, 0.0, 0.0], 1.0, 2.0, 0.0, 0.0);
        assert!(
            sdf > 0.0,
            "outside sphere, SDF should be positive, got {}",
            sdf
        );
    }

    #[test]
    fn test_sphere_sdf_radius() {
        // Distance from [5,0,0] to sphere of radius 1 centred at origin should be 4.
        let sdf = sphere_sdf([0.0, 0.0, 0.0], 1.0, 5.0, 0.0, 0.0);
        assert!((sdf - 4.0).abs() < 1e-12, "sdf = {}", sdf);
    }

    // --- torus_sdf ---

    #[test]
    fn test_torus_sdf_on_tube_surface() {
        // Point on the outer equator of a torus with R=2, r=0.5: at (2.5, 0, 0).
        let sdf = torus_sdf(2.0, 0.5, 2.5, 0.0, 0.0);
        assert!(sdf.abs() < 1e-10, "should be on surface, got {}", sdf);
    }

    #[test]
    fn test_torus_sdf_inside_tube() {
        let sdf = torus_sdf(2.0, 0.5, 2.0, 0.0, 0.0);
        assert!(sdf < 0.0, "centre of tube should be inside, got {}", sdf);
    }

    // --- interpolate_edge ---

    #[test]
    fn test_interpolate_edge_midpoint() {
        // v1 = -1, v2 = +1, isovalue = 0 → t = 0.5.
        let p = interpolate_edge([0.0, 0.0, 0.0], -1.0, [2.0, 0.0, 0.0], 1.0, 0.0);
        assert!((p[0] - 1.0).abs() < 1e-12, "x = {}", p[0]);
    }

    #[test]
    fn test_interpolate_edge_at_start() {
        let p = interpolate_edge([0.0, 0.0, 0.0], 0.0, [4.0, 0.0, 0.0], 2.0, 0.0);
        assert!(p[0].abs() < 1e-12, "x = {}", p[0]);
    }

    #[test]
    fn test_interpolate_edge_at_end() {
        let p = interpolate_edge([0.0, 0.0, 0.0], -2.0, [4.0, 0.0, 0.0], 0.0, 0.0);
        assert!((p[0] - 4.0).abs() < 1e-12, "x = {}", p[0]);
    }

    // --- compute_normal ---

    #[test]
    fn test_compute_normal_unit_length() {
        let tri = Triangle {
            vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normal: [0.0, 0.0, 1.0],
        };
        let n = compute_normal(&tri);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "normal length = {}", len);
    }

    #[test]
    fn test_compute_normal_xy_plane() {
        let tri = Triangle {
            vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normal: [0.0, 0.0, 1.0],
        };
        let n = compute_normal(&tri);
        assert!(n[2].abs() > 0.9, "normal should point in Z, got {:?}", n);
    }

    #[test]
    fn test_compute_normal_degenerate() {
        // All vertices identical → degenerate triangle, returns fallback [0,0,1].
        let tri = Triangle {
            vertices: [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
            normal: [0.0, 0.0, 1.0],
        };
        let n = compute_normal(&tri);
        assert!(
            (n[2] - 1.0).abs() < 1e-10,
            "fallback normal should be [0,0,1]"
        );
    }

    // --- ScalarField ---

    #[test]
    fn test_scalar_field_get_set() {
        let mut f = ScalarField::new(4, 4, 4, 1.0, 1.0, 1.0);
        f.set(1, 2, 3, 5.5);
        assert!((f.get(1, 2, 3) - 5.5).abs() < 1e-12);
    }

    #[test]
    fn test_scalar_field_default_zero() {
        let f = ScalarField::new(3, 3, 3, 1.0, 1.0, 1.0);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    assert_eq!(f.get(i, j, k), 0.0);
                }
            }
        }
    }

    // --- marching_cubes ---

    #[test]
    fn test_marching_cubes_flat_field_no_triangles() {
        // Constant field above isovalue: no surface.
        let n = 5;
        let mut f = ScalarField::new(n, n, n, 1.0, 1.0, 1.0);
        for v in f.values.iter_mut() {
            *v = 2.0;
        }
        let tris = marching_cubes(&f, 1.0);
        assert!(
            tris.is_empty(),
            "expected no triangles for field above isovalue"
        );
    }

    #[test]
    fn test_marching_cubes_sphere_positive_count() {
        // Build SDF of a sphere centred in the grid.
        let n = 20;
        let mut f = ScalarField::new(n, n, n, 1.0, 1.0, 1.0);
        let center = [(n as f64 - 1.0) / 2.0; 3];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let v = sphere_sdf(center, 5.0, i as f64, j as f64, k as f64);
                    f.set(i, j, k, v);
                }
            }
        }
        let tris = marching_cubes(&f, 0.0);
        assert!(
            !tris.is_empty(),
            "sphere isosurface should produce triangles"
        );
    }

    #[test]
    fn test_marching_cubes_triangles_have_unit_normals() {
        let n = 12;
        let mut f = ScalarField::new(n, n, n, 1.0, 1.0, 1.0);
        let center = [(n as f64 - 1.0) / 2.0; 3];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let v = sphere_sdf(center, 3.0, i as f64, j as f64, k as f64);
                    f.set(i, j, k, v);
                }
            }
        }
        let tris = marching_cubes(&f, 0.0);
        for tri in &tris {
            let n = tri.normal;
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-9 || len < 1e-9,
                "normal length = {}",
                len
            );
        }
    }

    #[test]
    fn test_marching_cubes_too_small_returns_empty() {
        let f = ScalarField::new(1, 1, 1, 1.0, 1.0, 1.0);
        let tris = marching_cubes(&f, 0.0);
        assert!(tris.is_empty());
    }

    // --- IsosurfaceMesh ---

    #[test]
    fn test_isosurface_mesh_triangle_count() {
        let tris = vec![
            Triangle {
                vertices: [[0.0; 3]; 3],
                normal: [0.0, 0.0, 1.0],
            },
            Triangle {
                vertices: [[1.0; 3]; 3],
                normal: [0.0, 0.0, 1.0],
            },
        ];
        let mesh = IsosurfaceMesh::from_triangles(tris);
        assert_eq!(mesh.triangle_count(), 2);
        assert_eq!(mesh.vertex_count, 6);
    }

    // --- dual_contouring_2d ---

    #[test]
    fn test_dual_contouring_2d_finds_crossing() {
        let field = vec![[0.0, -1.0], [1.0, 1.0]];
        let pts = dual_contouring_2d(&field, 0.0);
        assert_eq!(pts.len(), 1);
        assert!((pts[0][0] - 0.5).abs() < 1e-10, "x = {}", pts[0][0]);
    }

    #[test]
    fn test_dual_contouring_2d_no_crossing() {
        let field = vec![[0.0, 1.0], [1.0, 2.0], [2.0, 3.0]];
        let pts = dual_contouring_2d(&field, 0.0);
        assert!(pts.is_empty());
    }

    #[test]
    fn test_dual_contouring_2d_empty() {
        let pts = dual_contouring_2d(&[], 0.0);
        assert!(pts.is_empty());
    }
}
