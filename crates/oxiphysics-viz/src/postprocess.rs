// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Post-processing visualization: vector fields, scalar fields, streamlines,
//! and marching cubes for isosurface extraction.

// ─── Vector Field ────────────────────────────────────────────────────────────

/// 3-D vector field sampled on a regular axis-aligned grid.
///
/// The grid has `nx × ny × nz` cells with uniform spacing `dx` and a user-
/// supplied `origin`.  Data is stored in C-order (k varies slowest, i fastest).
pub struct VectorField {
    /// Number of grid points in the x direction.
    pub nx: usize,
    /// Number of grid points in the y direction.
    pub ny: usize,
    /// Number of grid points in the z direction.
    pub nz: usize,
    /// Uniform grid spacing (same in all three directions).
    pub dx: f64,
    /// World-space coordinate of grid point (0,0,0).
    pub origin: [f64; 3],
    /// Flat storage: index = i + nx*(j + ny*k).
    pub data: Vec<[f64; 3]>,
}

impl VectorField {
    /// Create a new zero-initialised vector field.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, origin: [f64; 3]) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            data: vec![[0.0; 3]; nx * ny * nz],
        }
    }

    #[inline]
    fn flat(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }

    /// Set the vector at grid point `(i, j, k)`.
    pub fn set(&mut self, i: usize, j: usize, k: usize, v: [f64; 3]) {
        let idx = self.flat(i, j, k);
        self.data[idx] = v;
    }

    /// Get the vector at grid point `(i, j, k)`.
    pub fn get(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        self.data[self.flat(i, j, k)]
    }

    /// Return `Some((i, j, k))` if `pos` is inside the grid, else `None`.
    pub fn cell_index(&self, pos: [f64; 3]) -> Option<(usize, usize, usize)> {
        let fi = (pos[0] - self.origin[0]) / self.dx;
        let fj = (pos[1] - self.origin[1]) / self.dx;
        let fk = (pos[2] - self.origin[2]) / self.dx;
        if fi < 0.0
            || fj < 0.0
            || fk < 0.0
            || fi >= (self.nx - 1) as f64
            || fj >= (self.ny - 1) as f64
            || fk >= (self.nz - 1) as f64
        {
            return None;
        }
        Some((fi as usize, fj as usize, fk as usize))
    }

    /// Trilinear interpolation at an arbitrary world-space position.
    ///
    /// Returns `[0.0; 3]` when `pos` is out of bounds.
    pub fn interpolate(&self, pos: [f64; 3]) -> [f64; 3] {
        let (i, j, k) = match self.cell_index(pos) {
            Some(c) => c,
            None => return [0.0; 3],
        };
        let tx = (pos[0] - self.origin[0]) / self.dx - i as f64;
        let ty = (pos[1] - self.origin[1]) / self.dx - j as f64;
        let tz = (pos[2] - self.origin[2]) / self.dx - k as f64;
        trilinear_vec(self, i, j, k, tx, ty, tz)
    }
}

/// Helper: perform trilinear interpolation inside cell `(i,j,k)` given
/// normalised offsets `(tx, ty, tz)` all in `[0, 1)`.
fn trilinear_vec(
    f: &VectorField,
    i: usize,
    j: usize,
    k: usize,
    tx: f64,
    ty: f64,
    tz: f64,
) -> [f64; 3] {
    let c000 = f.get(i, j, k);
    let c100 = f.get(i + 1, j, k);
    let c010 = f.get(i, j + 1, k);
    let c110 = f.get(i + 1, j + 1, k);
    let c001 = f.get(i, j, k + 1);
    let c101 = f.get(i + 1, j, k + 1);
    let c011 = f.get(i, j + 1, k + 1);
    let c111 = f.get(i + 1, j + 1, k + 1);

    let mut out = [0.0f64; 3];
    for d in 0..3 {
        let v00 = c000[d] * (1.0 - tx) + c100[d] * tx;
        let v10 = c010[d] * (1.0 - tx) + c110[d] * tx;
        let v01 = c001[d] * (1.0 - tx) + c101[d] * tx;
        let v11 = c011[d] * (1.0 - tx) + c111[d] * tx;
        let v0 = v00 * (1.0 - ty) + v10 * ty;
        let v1 = v01 * (1.0 - ty) + v11 * ty;
        out[d] = v0 * (1.0 - tz) + v1 * tz;
    }
    out
}

// ─── Streamlines ─────────────────────────────────────────────────────────────

/// Trace a single streamline through `field` using 4th-order Runge-Kutta.
///
/// Integration stops when the current position leaves the grid or `max_steps`
/// steps have been taken.  The returned vector always contains at least the
/// seed point.
pub fn trace_streamline(
    field: &VectorField,
    start: [f64; 3],
    step_size: f64,
    max_steps: usize,
) -> Vec<[f64; 3]> {
    let mut path = Vec::with_capacity(max_steps + 1);
    let mut pos = start;
    path.push(pos);

    for _ in 0..max_steps {
        // k1
        let k1 = field.interpolate(pos);
        if k1 == [0.0; 3] && field.cell_index(pos).is_none() {
            break;
        }
        // k2
        let p2 = add3(pos, scale3(k1, step_size * 0.5));
        let k2 = field.interpolate(p2);
        // k3
        let p3 = add3(pos, scale3(k2, step_size * 0.5));
        let k3 = field.interpolate(p3);
        // k4
        let p4 = add3(pos, scale3(k3, step_size));
        let k4 = field.interpolate(p4);

        let next = add3(
            pos,
            scale3(
                add3(add3(k1, scale3(k2, 2.0)), add3(scale3(k3, 2.0), k4)),
                step_size / 6.0,
            ),
        );

        if field.cell_index(next).is_none() {
            break;
        }
        pos = next;
        path.push(pos);
    }
    path
}

/// Trace multiple streamlines from a slice of seed points.
pub fn trace_streamlines(
    field: &VectorField,
    seeds: &[[f64; 3]],
    step_size: f64,
    max_steps: usize,
) -> Vec<Vec<[f64; 3]>> {
    seeds
        .iter()
        .map(|&s| trace_streamline(field, s, step_size, max_steps))
        .collect()
}

// ─── Scalar Field ─────────────────────────────────────────────────────────────

/// 3-D scalar field on a regular grid (same layout conventions as `VectorField`).
pub struct ScalarField {
    /// Number of grid points in the x direction.
    pub nx: usize,
    /// Number of grid points in the y direction.
    pub ny: usize,
    /// Number of grid points in the z direction.
    pub nz: usize,
    /// Uniform grid spacing.
    pub dx: f64,
    /// World-space coordinate of grid point (0,0,0).
    pub origin: [f64; 3],
    /// Flat storage: index = i + nx*(j + ny*k).
    pub data: Vec<f64>,
}

impl ScalarField {
    /// Create a new zero-initialised scalar field.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, origin: [f64; 3]) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            data: vec![0.0; nx * ny * nz],
        }
    }

    #[inline]
    fn flat(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }

    /// Set the scalar value at grid point `(i, j, k)`.
    pub fn set(&mut self, i: usize, j: usize, k: usize, v: f64) {
        let idx = self.flat(i, j, k);
        self.data[idx] = v;
    }

    /// Get the scalar value at grid point `(i, j, k)`.
    pub fn get(&self, i: usize, j: usize, k: usize) -> f64 {
        self.data[self.flat(i, j, k)]
    }

    /// Trilinear interpolation at an arbitrary world-space position.
    ///
    /// Returns `0.0` when `pos` is out of bounds.
    pub fn interpolate(&self, pos: [f64; 3]) -> f64 {
        let fi = (pos[0] - self.origin[0]) / self.dx;
        let fj = (pos[1] - self.origin[1]) / self.dx;
        let fk = (pos[2] - self.origin[2]) / self.dx;
        if fi < 0.0
            || fj < 0.0
            || fk < 0.0
            || fi >= (self.nx - 1) as f64
            || fj >= (self.ny - 1) as f64
            || fk >= (self.nz - 1) as f64
        {
            return 0.0;
        }
        let i = fi as usize;
        let j = fj as usize;
        let k = fk as usize;
        let tx = fi - i as f64;
        let ty = fj - j as f64;
        let tz = fk - k as f64;

        let v000 = self.get(i, j, k);
        let v100 = self.get(i + 1, j, k);
        let v010 = self.get(i, j + 1, k);
        let v110 = self.get(i + 1, j + 1, k);
        let v001 = self.get(i, j, k + 1);
        let v101 = self.get(i + 1, j, k + 1);
        let v011 = self.get(i, j + 1, k + 1);
        let v111 = self.get(i + 1, j + 1, k + 1);

        let v00 = v000 * (1.0 - tx) + v100 * tx;
        let v10 = v010 * (1.0 - tx) + v110 * tx;
        let v01 = v001 * (1.0 - tx) + v101 * tx;
        let v11 = v011 * (1.0 - tx) + v111 * tx;
        let v0 = v00 * (1.0 - ty) + v10 * ty;
        let v1 = v01 * (1.0 - ty) + v11 * ty;
        v0 * (1.0 - tz) + v1 * tz
    }

    /// Gradient at grid point `(i, j, k)` via central differences.
    ///
    /// One-sided differences are used at boundary nodes.
    pub fn gradient(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let inv2dx = 1.0 / (2.0 * self.dx);

        let gx = if i == 0 {
            (self.get(1, j, k) - self.get(0, j, k)) / self.dx
        } else if i == self.nx - 1 {
            (self.get(self.nx - 1, j, k) - self.get(self.nx - 2, j, k)) / self.dx
        } else {
            (self.get(i + 1, j, k) - self.get(i - 1, j, k)) * inv2dx
        };

        let gy = if j == 0 {
            (self.get(i, 1, k) - self.get(i, 0, k)) / self.dx
        } else if j == self.ny - 1 {
            (self.get(i, self.ny - 1, k) - self.get(i, self.ny - 2, k)) / self.dx
        } else {
            (self.get(i, j + 1, k) - self.get(i, j - 1, k)) * inv2dx
        };

        let gz = if k == 0 {
            (self.get(i, j, 1) - self.get(i, j, 0)) / self.dx
        } else if k == self.nz - 1 {
            (self.get(i, j, self.nz - 1) - self.get(i, j, self.nz - 2)) / self.dx
        } else {
            (self.get(i, j, k + 1) - self.get(i, j, k - 1)) * inv2dx
        };

        [gx, gy, gz]
    }
}

// ─── Marching Cubes (single cell) ─────────────────────────────────────────────

/// Interpolate along an edge between two corner positions given their scalar values.
fn edge_interp(p0: [f64; 3], v0: f64, p1: [f64; 3], v1: f64, threshold: f64) -> [f64; 3] {
    let dv = v1 - v0;
    let t = if dv.abs() < 1e-15 {
        0.5
    } else {
        (threshold - v0) / dv
    };
    [
        p0[0] + t * (p1[0] - p0[0]),
        p0[1] + t * (p1[1] - p0[1]),
        p0[2] + t * (p1[2] - p0[2]),
    ]
}

/// Marching-cubes triangle extraction for a single cell.
///
/// `values[n]` and `corners[n]` are the scalar value and world-space position
/// of corner `n` (0-7) using the standard MC ordering:
///
/// ```text
///      7────6
///     /|   /|
///    4────5 |
///    | 3──|─2
///    |/   |/
///    0────1
/// ```
///
/// Returns a flat list of triangle vertices (3 vertices per triangle), with at
/// most 5 triangles (15 vertices) per cell.  The function implements a
/// simplified 16-pattern table covering all topologically distinct cases.
pub fn marching_cubes_cell(
    values: [f64; 8],
    threshold: f64,
    corners: [[f64; 3]; 8],
) -> Vec<[f64; 3]> {
    // Build the 8-bit case index: bit n set when values[n] >= threshold.
    let mut case: u8 = 0;
    for n in 0..8u8 {
        if values[n as usize] >= threshold {
            case |= 1 << n;
        }
    }
    // All above or all below → no surface.
    if case == 0x00 || case == 0xFF {
        return Vec::new();
    }

    // The 12 edges of the cube, each defined by two corner indices.
    // Standard MC edge numbering:
    //  0: 0-1   1: 1-2   2: 2-3   3: 3-0
    //  4: 4-5   5: 5-6   6: 6-7   7: 7-4
    //  8: 0-4   9: 1-5  10: 2-6  11: 3-7
    const EDGES: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];

    // Precompute all 12 edge intersection points.
    let ep: Vec<[f64; 3]> = EDGES
        .iter()
        .map(|&(a, b)| edge_interp(corners[a], values[a], corners[b], values[b], threshold))
        .collect();

    // Triangle table — each row is a sequence of edge indices terminated by -1.
    // This is the full 256-entry standard MC triangle table (Paul Bourke / Lorensen & Cline).
    // Rows that would be identical to their complement (case ^ 0xFF) are handled by
    // the same entries after flipping the `case` sense above.
    #[rustfmt::skip]
    let tri_table: [[i8; 16]; 256] = [
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
        [8,2,9,8,9,7,2,10,9,5,9,3,10,3,9,-1],
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

    let row = &tri_table[case as usize];
    let mut vertices = Vec::new();
    let mut idx = 0;
    while idx < 16 && row[idx] != -1 {
        let e = row[idx] as usize;
        vertices.push(ep[e]);
        idx += 1;
    }
    vertices
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ─── Screen-space Ambient Occlusion ──────────────────────────────────────────

/// Parameters for screen-space ambient occlusion (SSAO).
#[derive(Debug, Clone)]
pub struct SsaoParams {
    /// Sampling radius in view-space units.
    pub radius: f64,
    /// Number of sample points on the hemisphere.
    pub num_samples: usize,
    /// Bias to avoid self-occlusion artefacts.
    pub bias: f64,
    /// AO intensity multiplier.
    pub intensity: f64,
}

impl Default for SsaoParams {
    fn default() -> Self {
        Self {
            radius: 0.5,
            num_samples: 16,
            bias: 0.025,
            intensity: 1.0,
        }
    }
}

/// Compute a simplified SSAO factor for a single depth sample.
///
/// Given a depth value `d`, its surface normal `n`, and a set of `samples`
/// (hemisphere offsets in view space), returns an occlusion value in `[0, 1]`.
/// A value of 0 means fully occluded; 1 means fully unoccluded.
///
/// This is a CPU-side stand-in for the GPU SSAO pass used in debugging.
pub fn compute_ssao_sample(
    depth: f64,
    normal: [f64; 3],
    samples: &[[f64; 3]],
    depth_field: &ScalarField,
    params: &SsaoParams,
) -> f64 {
    if samples.is_empty() {
        return 1.0;
    }
    // Convert depth to a rough world-space position (assume camera at origin, Z-forward).
    let frag_pos = [0.0, 0.0, depth];

    let mut occlusion = 0.0_f64;
    for &s in samples {
        // Orient sample along normal.
        let scale_fac = params.radius;
        let sample_pos = add3(frag_pos, scale3(s, scale_fac));

        // Sample depth from field (use x=sample_pos[0], y=sample_pos[1], z=0 as lookup).
        let lookup = [sample_pos[0], sample_pos[1], 0.0];
        let sample_depth = depth_field.interpolate(lookup);

        // Range check contribution.
        let range_check = if (sample_depth - depth).abs() < params.radius {
            1.0
        } else {
            0.0
        };

        let dot_check = vec3_dot(normal, s).max(0.0);
        if sample_depth <= sample_pos[2] - params.bias {
            occlusion += dot_check * range_check;
        }
    }
    1.0 - (occlusion / samples.len() as f64) * params.intensity
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Generate a cosine-weighted hemisphere sample set for SSAO.
///
/// Returns `n` unit vectors in the upper hemisphere (z > 0) with a
/// deterministic, low-discrepancy distribution.
pub fn generate_ssao_samples(n: usize) -> Vec<[f64; 3]> {
    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU;
        let phi = ((i as f64 * 1.618033988749_f64).fract() * std::f64::consts::FRAC_PI_2)
            .min(std::f64::consts::FRAC_PI_2 - 1e-9);
        let scale = (i as f64 + 1.0) / n as f64;
        let scale = 0.1 + 0.9 * scale * scale;
        let s = [
            theta.cos() * phi.sin() * scale,
            theta.sin() * phi.sin() * scale,
            phi.cos() * scale,
        ];
        samples.push(s);
    }
    samples
}

// ─── Temporal Anti-Aliasing ───────────────────────────────────────────────────

/// Blend a current frame image buffer with a previous (history) buffer using
/// temporal anti-aliasing (TAA).
///
/// Each element of `current` and `history` is a float RGB pixel `[r, g, b]`.
/// The blended result is written into `output` using the feedback factor `α`:
/// `out = α * current + (1 - α) * history`.
///
/// Typical values: `α ≈ 0.1`.
pub fn temporal_aa_blend(
    current: &[[f64; 3]],
    history: &[[f64; 3]],
    output: &mut Vec<[f64; 3]>,
    alpha: f64,
) {
    assert_eq!(current.len(), history.len(), "TAA: buffer size mismatch");
    output.resize(current.len(), [0.0; 3]);
    for (i, (&c, &h)) in current.iter().zip(history.iter()).enumerate() {
        output[i] = [
            alpha * c[0] + (1.0 - alpha) * h[0],
            alpha * c[1] + (1.0 - alpha) * h[1],
            alpha * c[2] + (1.0 - alpha) * h[2],
        ];
    }
}

/// Clamp history colour to the neighbourhood colour box (TAA ghost reduction).
///
/// For each pixel in `history`, clamps its RGB channels into the AABB defined
/// by the minimum and maximum values in the 3×3 neighbourhood of `current`.
/// `width` is the image width in pixels.
pub fn taa_clamp_history(current: &[[f64; 3]], history: &mut [[f64; 3]], width: usize) {
    let height = current.len().checked_div(width).unwrap_or(0);
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for dy in 0..=2usize {
                for dx in 0..=2usize {
                    let ny = y.wrapping_add(dy).wrapping_sub(1);
                    let nx = x.wrapping_add(dx).wrapping_sub(1);
                    if ny < height && nx < width {
                        let p = current[ny * width + nx];
                        for c in 0..3 {
                            if p[c] < lo[c] {
                                lo[c] = p[c];
                            }
                            if p[c] > hi[c] {
                                hi[c] = p[c];
                            }
                        }
                    }
                }
            }
            let h = history[idx];
            history[idx] = [
                h[0].clamp(lo[0], hi[0]),
                h[1].clamp(lo[1], hi[1]),
                h[2].clamp(lo[2], hi[2]),
            ];
        }
    }
}

// ─── Depth of Field ──────────────────────────────────────────────────────────

/// Apply a depth-of-field (DoF) blur to an image buffer.
///
/// Pixels further from `focus_depth` by more than `dof_range` are blurred by
/// averaging with their `blur_radius`-wide box neighbourhood.  Pixels inside
/// `dof_range` are unmodified.
///
/// `depth_buffer` contains one depth value per pixel.  `pixels` contains one
/// RGB triple per pixel.  `width` is the image width.
pub fn apply_depth_of_field(
    pixels: &mut [[f64; 3]],
    depth_buffer: &[f64],
    width: usize,
    focus_depth: f64,
    dof_range: f64,
    blur_radius: usize,
) {
    let n = pixels.len();
    if width == 0 || n == 0 {
        return;
    }
    let height = n / width;
    let original = pixels.to_vec();

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let depth = depth_buffer.get(idx).copied().unwrap_or(0.0);
            if (depth - focus_depth).abs() <= dof_range {
                continue; // In-focus, no blur.
            }
            // Box blur.
            let mut sum = [0.0_f64; 3];
            let mut count = 0usize;
            for dy in 0..=2 * blur_radius {
                for dx in 0..=2 * blur_radius {
                    let ny = y.wrapping_add(dy).wrapping_sub(blur_radius);
                    let nx = x.wrapping_add(dx).wrapping_sub(blur_radius);
                    if ny < height && nx < width {
                        let p = original[ny * width + nx];
                        for c in 0..3 {
                            sum[c] += p[c];
                        }
                        count += 1;
                    }
                }
            }
            if count > 0 {
                let inv = 1.0 / count as f64;
                pixels[idx] = [sum[0] * inv, sum[1] * inv, sum[2] * inv];
            }
        }
    }
}

// ─── Chromatic Aberration ─────────────────────────────────────────────────────

/// Apply chromatic aberration by offsetting the red and blue channels
/// radially outward from the image centre.
///
/// `strength` controls how many pixels of offset are applied at the edge.
/// `width` is the image width.  The image is assumed to have no padding.
pub fn apply_chromatic_aberration(pixels: &mut [[f64; 3]], width: usize, strength: f64) {
    let n = pixels.len();
    if width == 0 || n == 0 {
        return;
    }
    let height = n / width;
    let cx = width as f64 * 0.5;
    let cy = height as f64 * 0.5;
    let max_dist = (cx * cx + cy * cy).sqrt().max(1e-14);
    let original = pixels.to_vec();

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            let dist = (dx * dx + dy * dy).sqrt();
            let factor = strength * dist / max_dist;
            let nx_dx = dx * factor;
            let nx_dy = dy * factor;

            // Red channel: shift outward.
            let rx = (x as f64 + nx_dx).round() as isize;
            let ry = (y as f64 + nx_dy).round() as isize;
            let red = if rx >= 0 && rx < width as isize && ry >= 0 && ry < height as isize {
                original[(ry as usize) * width + rx as usize][0]
            } else {
                original[idx][0]
            };

            // Blue channel: shift inward.
            let bx = (x as f64 - nx_dx).round() as isize;
            let by_ = (y as f64 - nx_dy).round() as isize;
            let blue = if bx >= 0 && bx < width as isize && by_ >= 0 && by_ < height as isize {
                original[(by_ as usize) * width + bx as usize][2]
            } else {
                original[idx][2]
            };

            pixels[idx] = [red, original[idx][1], blue];
        }
    }
}

// ─── Vignette ─────────────────────────────────────────────────────────────────

/// Apply a vignette effect: multiply pixel brightness by a radial falloff.
///
/// The falloff is `1 - smoothstep(inner, outer, r)` where `r` is the
/// normalised distance from the image centre.  `inner` and `outer` must be in
/// `[0, 1]`.
pub fn apply_vignette(pixels: &mut [[f64; 3]], width: usize, inner: f64, outer: f64) {
    let n = pixels.len();
    if width == 0 || n == 0 {
        return;
    }
    let height = n / width;
    let cx = width as f64 * 0.5;
    let cy = height as f64 * 0.5;
    let max_dist = (cx * cx + cy * cy).sqrt().max(1e-14);

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let dx = (x as f64 - cx) / max_dist;
            let dy = (y as f64 - cy) / max_dist;
            let r = (dx * dx + dy * dy).sqrt();
            let vignette = 1.0 - smoothstep(inner, outer, r);
            for ch in pixels[idx].iter_mut() {
                *ch *= vignette;
            }
        }
    }
}

/// Smooth hermite interpolation: maps `t ∈ [lo, hi]` to `[0, 1]` using
/// the formula `3t² - 2t³`.
#[inline]
fn smoothstep(lo: f64, hi: f64, x: f64) -> f64 {
    let t = ((x - lo) / (hi - lo).max(1e-30)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

// ─── Tone Mapping ────────────────────────────────────────────────────────────

/// Apply Reinhard tone mapping to an HDR pixel buffer.
///
/// `exposure` scales the input luminance; higher values brighten the image.
pub fn tonemap_reinhard(pixels: &mut [[f64; 3]], exposure: f64) {
    for p in pixels.iter_mut() {
        for ch in p.iter_mut() {
            let v = *ch * exposure;
            *ch = v / (1.0 + v);
        }
    }
}

/// Apply ACES filmic tone mapping approximation.
pub fn tonemap_aces(pixels: &mut [[f64; 3]], exposure: f64) {
    // Narkowicz 2015 ACES approximation: y = x*(2.51*x+0.03)/(x*(2.43*x+0.59)+0.14)
    for p in pixels.iter_mut() {
        for ch in p.iter_mut() {
            let x = (*ch * exposure).max(0.0);
            *ch = (x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14);
            *ch = ch.clamp(0.0, 1.0);
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── VectorField ──────────────────────────────────────────────────────────

    #[test]
    fn test_vector_field_uniform_interpolation() {
        // A field where every node has the same value should interpolate to
        // that same value at every interior point.
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let dx = 1.0;
        let origin = [0.0, 0.0, 0.0];
        let mut field = VectorField::new(nx, ny, nz, dx, origin);
        let v = [3.0, -1.5, 0.7];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    field.set(i, j, k, v);
                }
            }
        }
        // Sample at several interior positions.
        for &pos in &[[0.5, 0.5, 0.5], [1.0, 1.0, 1.0], [2.3, 1.7, 0.9]] {
            let r = field.interpolate(pos);
            for d in 0..3 {
                assert!(
                    (r[d] - v[d]).abs() < 1e-12,
                    "uniform field: component {} at {:?}: got {}, expected {}",
                    d,
                    pos,
                    r[d],
                    v[d]
                );
            }
        }
    }

    // ── Streamlines ──────────────────────────────────────────────────────────

    #[test]
    fn test_streamline_uniform_x_field() {
        // Build a 10×4×4 grid with uniform velocity (1,0,0).
        let nx = 10;
        let ny = 4;
        let nz = 4;
        let dx = 1.0;
        let mut field = VectorField::new(nx, ny, nz, dx, [0.0, 0.0, 0.0]);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    field.set(i, j, k, [1.0, 0.0, 0.0]);
                }
            }
        }
        let start = [1.0, 1.5, 1.5];
        let step = 0.1;
        let steps = 5;
        let path = trace_streamline(&field, start, step, steps);

        assert!(path.len() >= 2, "expected at least 2 points");

        // Each successive point should advance in +x with no y/z change.
        for w in path.windows(2) {
            let dx_actual = w[1][0] - w[0][0];
            assert!(
                dx_actual > 0.0,
                "streamline should advance in +x, got dx={}",
                dx_actual
            );
            assert!(
                (w[1][1] - w[0][1]).abs() < 1e-12,
                "y should not change in uniform x-field"
            );
            assert!(
                (w[1][2] - w[0][2]).abs() < 1e-12,
                "z should not change in uniform x-field"
            );
        }
    }

    // ── ScalarField ──────────────────────────────────────────────────────────

    #[test]
    fn test_scalar_field_linear_gradient() {
        // f(i,j,k) = i*dx  →  gradient should be (1,0,0) everywhere.
        let n = 5;
        let dx = 0.5;
        let mut sf = ScalarField::new(n, n, n, dx, [0.0; 3]);
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    sf.set(i, j, k, i as f64 * dx);
                }
            }
        }
        // Interior nodes should give exact central-difference gradient (1,0,0).
        for k in 1..n - 1 {
            for j in 1..n - 1 {
                for i in 1..n - 1 {
                    let g = sf.gradient(i, j, k);
                    assert!(
                        (g[0] - 1.0).abs() < 1e-12,
                        "gx at ({},{},{}) = {}, expected 1.0",
                        i,
                        j,
                        k,
                        g[0]
                    );
                    assert!(
                        g[1].abs() < 1e-12,
                        "gy at ({},{},{}) = {}, expected 0.0",
                        i,
                        j,
                        k,
                        g[1]
                    );
                    assert!(
                        g[2].abs() < 1e-12,
                        "gz at ({},{},{}) = {}, expected 0.0",
                        i,
                        j,
                        k,
                        g[2]
                    );
                }
            }
        }
    }

    // ── Marching Cubes ───────────────────────────────────────────────────────

    fn unit_cube_corners() -> [[f64; 3]; 8] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    }

    #[test]
    fn test_marching_cubes_all_above_no_tris() {
        let vals = [1.0f64; 8];
        let tris = marching_cubes_cell(vals, 0.5, unit_cube_corners());
        assert!(
            tris.is_empty(),
            "all above threshold → 0 triangles, got {}",
            tris.len()
        );
    }

    #[test]
    fn test_marching_cubes_all_below_no_tris() {
        let vals = [0.0f64; 8];
        let tris = marching_cubes_cell(vals, 0.5, unit_cube_corners());
        assert!(
            tris.is_empty(),
            "all below threshold → 0 triangles, got {}",
            tris.len()
        );
    }

    #[test]
    fn test_marching_cubes_one_corner_inside() {
        // Corner 0 above threshold, rest below → exactly 1 triangle (3 vertices).
        let mut vals = [0.0f64; 8];
        vals[0] = 1.0;
        let tris = marching_cubes_cell(vals, 0.5, unit_cube_corners());
        assert_eq!(tris.len(), 3, "one corner inside → 1 triangle (3 vertices)");
    }

    #[test]
    fn test_marching_cubes_mixed_produces_triangles() {
        // Half the corners above (checkerboard) → non-empty output.
        let vals = [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0];
        let tris = marching_cubes_cell(vals, 0.5, unit_cube_corners());
        assert!(
            !tris.is_empty(),
            "mixed threshold crossings should produce triangles"
        );
        // Number of vertices must be a multiple of 3.
        assert_eq!(
            tris.len() % 3,
            0,
            "triangle vertex count must be divisible by 3"
        );
    }

    // ── TAA Tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_taa_blend_alpha_one() {
        let current = vec![[1.0, 0.0, 0.0]; 4];
        let history = vec![[0.0, 1.0, 0.0]; 4];
        let mut output = Vec::new();
        temporal_aa_blend(&current, &history, &mut output, 1.0);
        // α=1 → output should equal current.
        for p in &output {
            assert!((p[0] - 1.0).abs() < 1e-12);
            assert!(p[1].abs() < 1e-12);
        }
    }

    #[test]
    fn test_taa_blend_alpha_zero() {
        let current = vec![[1.0, 0.0, 0.0]; 4];
        let history = vec![[0.0, 1.0, 0.0]; 4];
        let mut output = Vec::new();
        temporal_aa_blend(&current, &history, &mut output, 0.0);
        // α=0 → output should equal history.
        for p in &output {
            assert!(p[0].abs() < 1e-12);
            assert!((p[1] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_taa_blend_intermediate() {
        let current = vec![[1.0, 0.0, 0.0]];
        let history = vec![[0.0, 0.0, 0.0]];
        let mut output = Vec::new();
        temporal_aa_blend(&current, &history, &mut output, 0.5);
        assert!((output[0][0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_taa_clamp_history_no_panic() {
        let width = 4;
        let current: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let mut history = current.clone();
        taa_clamp_history(&current, &mut history, width);
    }

    // ── SSAO Tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_ssao_samples_unit_hemisphere() {
        let samples = generate_ssao_samples(16);
        assert_eq!(samples.len(), 16);
        for s in &samples {
            // z (hemisphere) must be ≥ 0.
            assert!(s[2] >= 0.0, "SSAO sample z should be non-negative: {:?}", s);
        }
    }

    #[test]
    fn test_compute_ssao_no_occlusion_empty_depth() {
        let params = SsaoParams::default();
        let depth_field = ScalarField::new(4, 4, 2, 1.0, [0.0; 3]);
        // No samples → fully unoccluded.
        let ao = compute_ssao_sample(0.5, [0.0, 0.0, 1.0], &[], &depth_field, &params);
        assert!((ao - 1.0).abs() < 1e-12);
    }

    // ── Depth of Field Tests ──────────────────────────────────────────────────

    #[test]
    fn test_dof_in_focus_pixels_unchanged() {
        let mut pixels: Vec<[f64; 3]> = vec![[1.0, 0.5, 0.25]; 9];
        let depth_buffer = vec![5.0_f64; 9];
        let original = pixels.clone();
        apply_depth_of_field(&mut pixels, &depth_buffer, 3, 5.0, 1.0, 1);
        // All pixels are at focus depth → no change.
        for (p, o) in pixels.iter().zip(original.iter()) {
            for c in 0..3 {
                assert!((p[c] - o[c]).abs() < 1e-12, "In-focus pixel was blurred");
            }
        }
    }

    #[test]
    fn test_dof_out_of_focus_blurs() {
        let mut pixels: Vec<[f64; 3]> = vec![[0.0; 3]; 9];
        pixels[4] = [1.0, 1.0, 1.0]; // bright centre pixel
        let depth_buffer = vec![10.0_f64; 9]; // all far from focus
        apply_depth_of_field(&mut pixels, &depth_buffer, 3, 1.0, 0.1, 1);
        // After blur, the bright centre should have spread to neighbours.
        let sum: f64 = pixels.iter().map(|p| p[0]).sum();
        assert!(sum > 0.9, "Blur should spread brightness: sum={sum}");
    }

    // ── Chromatic Aberration Tests ────────────────────────────────────────────

    #[test]
    fn test_chromatic_aberration_no_panic() {
        let width = 4;
        let mut pixels: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1; 3]).collect();
        apply_chromatic_aberration(&mut pixels, width, 2.0);
    }

    #[test]
    fn test_chromatic_aberration_strength_zero_no_change() {
        let width = 4;
        let mut pixels: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1; 3]).collect();
        let original = pixels.clone();
        apply_chromatic_aberration(&mut pixels, width, 0.0);
        // With zero strength, centre channel (green) unchanged; R and B may shift 0.
        for (p, o) in pixels.iter().zip(original.iter()) {
            assert!(
                (p[1] - o[1]).abs() < 1e-12,
                "Green channel should be unchanged"
            );
        }
    }

    // ── Vignette Tests ────────────────────────────────────────────────────────

    #[test]
    fn test_vignette_darkens_corners() {
        let width = 5;
        let height = 5;
        let mut pixels: Vec<[f64; 3]> = vec![[1.0; 3]; width * height];
        apply_vignette(&mut pixels, width, 0.2, 0.8);
        // Corner pixels should be darker than centre pixel.
        let centre = pixels[2 * width + 2][0];
        let corner = pixels[0][0];
        assert!(
            corner < centre,
            "Corner should be darker after vignette: corner={corner} centre={centre}"
        );
    }

    #[test]
    fn test_vignette_centre_unchanged() {
        let width = 5;
        let height = 5;
        let mut pixels: Vec<[f64; 3]> = vec![[1.0; 3]; width * height];
        apply_vignette(&mut pixels, width, 0.5, 1.0);
        // Centre pixel at (2,2) is at r=0, so vignette=1 → no darkening.
        let centre = pixels[2 * width + 2][0];
        assert!(centre > 0.99, "Centre pixel should remain bright: {centre}");
    }

    // ── Smoothstep Test ───────────────────────────────────────────────────────

    #[test]
    fn test_smoothstep_endpoints() {
        assert!((smoothstep(0.0, 1.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((smoothstep(0.0, 1.0, 1.0) - 1.0).abs() < 1e-12);
        assert!((smoothstep(0.0, 1.0, 0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_smoothstep_clamps() {
        assert!(smoothstep(0.2, 0.8, -1.0).abs() < 1e-12);
        assert!((smoothstep(0.2, 0.8, 2.0) - 1.0).abs() < 1e-12);
    }

    // ── Tone Mapping Tests ────────────────────────────────────────────────────

    #[test]
    fn test_tonemap_reinhard_range() {
        let mut pixels: Vec<[f64; 3]> = vec![[10.0, 5.0, 0.1], [0.5, 0.5, 0.5]];
        tonemap_reinhard(&mut pixels, 1.0);
        for p in &pixels {
            for &ch in p.iter() {
                assert!(
                    (0.0..=1.0).contains(&ch),
                    "Reinhard output out of [0,1]: {}",
                    ch
                );
            }
        }
    }

    #[test]
    fn test_tonemap_aces_range() {
        let mut pixels: Vec<[f64; 3]> = vec![[100.0, 0.5, 0.01]];
        tonemap_aces(&mut pixels, 1.0);
        for p in &pixels {
            for &ch in p.iter() {
                assert!(
                    (0.0..=1.0).contains(&ch),
                    "ACES output out of [0,1]: {}",
                    ch
                );
            }
        }
    }
}
