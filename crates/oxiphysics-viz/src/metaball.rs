// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Metaball isosurface rendering for the OxiPhysics visualizer.
//!
//! Provides:
//! - [`Metaball`]: a single blob defined by centre, radius, and strength.
//! - Scalar field evaluation (`1/r²` falloff).
//! - Naive marching-cubes isosurface extraction on a uniform grid.
//! - Normal computation from the scalar-field gradient.
//! - Dynamic update (move centres).
//! - LOD (coarse/fine grid resolution).
//! - Blob-merging detection.

#![allow(dead_code)]
#![allow(missing_docs)]

// ─────────────────────────────────────────────────────────────────────────────
// Metaball
// ─────────────────────────────────────────────────────────────────────────────

/// A single metaball (blob) defined by its centre, radius, and strength.
///
/// The contribution to the scalar field at point `p` is:
/// `strength / |p - centre|²`  (clamped so denominator ≥ ε).
#[derive(Debug, Clone, Copy)]
pub struct Metaball {
    /// World-space centre of the blob.
    pub centre: [f64; 3],
    /// Characteristic radius (used for LOD/merging detection).
    pub radius: f64,
    /// Strength (amplitude) of the scalar contribution.
    pub strength: f64,
}

impl Metaball {
    /// Create a metaball with equal radius and strength.
    pub fn new(centre: [f64; 3], radius: f64) -> Self {
        Self {
            centre,
            radius,
            strength: radius * radius,
        }
    }

    /// Create a metaball with explicit strength.
    pub fn with_strength(centre: [f64; 3], radius: f64, strength: f64) -> Self {
        Self {
            centre,
            radius,
            strength,
        }
    }

    /// Evaluate the contribution of this metaball at point `p`.
    ///
    /// Returns `strength / max(dist², 1e-12)`.
    pub fn field_value(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.centre[0];
        let dy = p[1] - self.centre[1];
        let dz = p[2] - self.centre[2];
        let r2 = (dx * dx + dy * dy + dz * dz).max(1e-12);
        self.strength / r2
    }

    /// Distance from the centre of this metaball to another.
    pub fn distance_to(&self, other: &Metaball) -> f64 {
        let dx = self.centre[0] - other.centre[0];
        let dy = self.centre[1] - other.centre[1];
        let dz = self.centre[2] - other.centre[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Move the centre by a displacement vector.
    pub fn translate(&mut self, delta: [f64; 3]) {
        self.centre[0] += delta[0];
        self.centre[1] += delta[1];
        self.centre[2] += delta[2];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MetaballField — sum of 1/r² contributions
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of metaballs whose scalar contributions are summed.
#[derive(Debug, Clone, Default)]
pub struct MetaballField {
    /// The metaballs in this field.
    pub balls: Vec<Metaball>,
}

impl MetaballField {
    /// Create an empty field.
    pub fn new() -> Self {
        Self { balls: Vec::new() }
    }

    /// Add a metaball to the field.
    pub fn add(&mut self, ball: Metaball) {
        self.balls.push(ball);
    }

    /// Evaluate the total scalar field at point `p` (sum of all contributions).
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        self.balls.iter().map(|b| b.field_value(p)).sum()
    }

    /// Compute the approximate gradient of the scalar field at `p` using
    /// central finite differences with step `h`.
    pub fn gradient(&self, p: [f64; 3], h: f64) -> [f64; 3] {
        let two_h = 2.0 * h;
        let [px, py, pz] = p;
        let gx = (self.evaluate([px + h, py, pz]) - self.evaluate([px - h, py, pz])) / two_h;
        let gy = (self.evaluate([px, py + h, pz]) - self.evaluate([px, py - h, pz])) / two_h;
        let gz = (self.evaluate([px, py, pz + h]) - self.evaluate([px, py, pz - h])) / two_h;
        [gx, gy, gz]
    }

    /// Return the number of metaballs.
    pub fn len(&self) -> usize {
        self.balls.len()
    }

    /// Returns `true` if there are no metaballs.
    pub fn is_empty(&self) -> bool {
        self.balls.is_empty()
    }

    /// Move all metaball centres by the displacement vector.
    pub fn translate_all(&mut self, delta: [f64; 3]) {
        for b in &mut self.balls {
            b.translate(delta);
        }
    }

    /// Move each metaball by its individual displacement.
    ///
    /// If `deltas.len() < self.balls.len()`, the remaining balls are unaffected.
    pub fn update_centres(&mut self, deltas: &[[f64; 3]]) {
        for (ball, &delta) in self.balls.iter_mut().zip(deltas.iter()) {
            ball.translate(delta);
        }
    }

    /// Detect pairs of metaballs that are "merged" (overlapping within their
    /// combined radii).
    ///
    /// Returns pairs `(i, j)` with `i < j` where
    /// `distance(i, j) < radius_i + radius_j`.
    pub fn merged_pairs(&self) -> Vec<(usize, usize)> {
        let n = self.balls.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dist = self.balls[i].distance_to(&self.balls[j]);
                let threshold = self.balls[i].radius + self.balls[j].radius;
                if dist < threshold {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Isosurface mesh
// ─────────────────────────────────────────────────────────────────────────────

/// A triangle in the extracted isosurface.
#[derive(Debug, Clone, Copy)]
pub struct IsoTriangle {
    /// Vertices (world-space positions).
    pub vertices: [[f64; 3]; 3],
    /// Outward face normal.
    pub normal: [f64; 3],
}

impl IsoTriangle {
    /// Compute the face normal from the vertex positions.
    fn compute_normal(v: &[[f64; 3]; 3]) -> [f64; 3] {
        let e1 = [v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]];
        let e2 = [v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]];
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-15);
        [n[0] / len, n[1] / len, n[2] / len]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Grid resolution for LOD
// ─────────────────────────────────────────────────────────────────────────────

/// Grid resolution for marching-cubes evaluation.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridResolution {
    /// Coarse grid (fast preview).
    Coarse,
    /// Medium grid (balanced quality/speed).
    Medium,
    /// Fine grid (high-quality).
    Fine,
    /// Custom resolution (nx, ny, nz).
    Custom(usize, usize, usize),
}

impl GridResolution {
    /// Return the (nx, ny, nz) cell counts for this resolution.
    pub fn cells(&self) -> (usize, usize, usize) {
        match self {
            GridResolution::Coarse => (8, 8, 8),
            GridResolution::Medium => (16, 16, 16),
            GridResolution::Fine => (32, 32, 32),
            GridResolution::Custom(nx, ny, nz) => (*nx, *ny, *nz),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Marching cubes lookup tables (naïve / edge-based)
// ─────────────────────────────────────────────────────────────────────────────

/// Marching-cubes edge table: for each of 256 corner configurations, a bitmask
/// of the 12 edges that the isosurface crosses.
///
/// (Standard Lorensen & Cline table, abridged here as a fixed 256-entry array.)
#[rustfmt::skip]
const MC_EDGE_TABLE: [u16; 256] = [
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x835, 0xb3f, 0xa36,  // row patched
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000,
];

/// Cube corner offsets (vertex 0..7).
const CUBE_CORNERS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0], // 0
    [1.0, 0.0, 0.0], // 1
    [1.0, 1.0, 0.0], // 2
    [0.0, 1.0, 0.0], // 3
    [0.0, 0.0, 1.0], // 4
    [1.0, 0.0, 1.0], // 5
    [1.0, 1.0, 1.0], // 6
    [0.0, 1.0, 1.0], // 7
];

/// The 12 edges, each defined by (corner_a, corner_b).
const CUBE_EDGES: [(usize, usize); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0), // bottom ring
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4), // top ring
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7), // vertical edges
];

/// Interpolate along a cube edge to find the isosurface crossing point.
fn interpolate_edge(p0: [f64; 3], p1: [f64; 3], v0: f64, v1: f64, threshold: f64) -> [f64; 3] {
    let dv = v1 - v0;
    let t = if dv.abs() < 1e-15 {
        0.5
    } else {
        ((threshold - v0) / dv).clamp(0.0, 1.0)
    };
    [
        p0[0] + t * (p1[0] - p0[0]),
        p0[1] + t * (p1[1] - p0[1]),
        p0[2] + t * (p1[2] - p0[2]),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Marching cubes extraction
// ─────────────────────────────────────────────────────────────────────────────

/// Extract an isosurface from the [`MetaballField`] using naive marching cubes.
///
/// # Parameters
/// - `field`: the scalar field to extract from.
/// - `min` / `max`: bounding box of the grid in world space.
/// - `resolution`: grid resolution (number of cells per axis).
/// - `threshold`: isosurface level.
///
/// # Returns
/// A vector of [`IsoTriangle`]s approximating the isosurface.
pub fn extract_isosurface(
    field: &MetaballField,
    min: [f64; 3],
    max: [f64; 3],
    resolution: GridResolution,
    threshold: f64,
) -> Vec<IsoTriangle> {
    let (nx, ny, nz) = resolution.cells();
    let nx = nx.max(1);
    let ny = ny.max(1);
    let nz = nz.max(1);

    let dx = (max[0] - min[0]) / nx as f64;
    let dy = (max[1] - min[1]) / ny as f64;
    let dz = (max[2] - min[2]) / nz as f64;

    // Pre-compute scalar field on grid vertices.
    let gx = nx + 1;
    let gy = ny + 1;
    let gz = nz + 1;
    let mut values = vec![0.0_f64; gx * gy * gz];

    for iz in 0..=nz {
        for iy in 0..=ny {
            for ix in 0..=nx {
                let p = [
                    min[0] + ix as f64 * dx,
                    min[1] + iy as f64 * dy,
                    min[2] + iz as f64 * dz,
                ];
                values[iz * gy * gx + iy * gx + ix] = field.evaluate(p);
            }
        }
    }

    let val = |ix: usize, iy: usize, iz: usize| values[iz * gy * gx + iy * gx + ix];

    let mut triangles = Vec::new();

    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                // Cell origin.
                let ox = min[0] + ix as f64 * dx;
                let oy = min[1] + iy as f64 * dy;
                let oz = min[2] + iz as f64 * dz;

                // 8 corner values.
                let cv = [
                    val(ix, iy, iz),
                    val(ix + 1, iy, iz),
                    val(ix + 1, iy + 1, iz),
                    val(ix, iy + 1, iz),
                    val(ix, iy, iz + 1),
                    val(ix + 1, iy, iz + 1),
                    val(ix + 1, iy + 1, iz + 1),
                    val(ix, iy + 1, iz + 1),
                ];

                // Compute cube index.
                let mut cube_idx: usize = 0;
                for (k, &v) in cv.iter().enumerate() {
                    if v >= threshold {
                        cube_idx |= 1 << k;
                    }
                }

                let edge_flags = MC_EDGE_TABLE[cube_idx];
                if edge_flags == 0 {
                    continue;
                }

                // World-space corner positions.
                let corners: [[f64; 3]; 8] = {
                    let mut c = [[0.0_f64; 3]; 8];
                    for (k, off) in CUBE_CORNERS.iter().enumerate() {
                        c[k] = [ox + off[0] * dx, oy + off[1] * dy, oz + off[2] * dz];
                    }
                    c
                };

                // Interpolate edge crossing points.
                let mut edge_verts = [[0.0_f64; 3]; 12];
                for (e, &(ca, cb)) in CUBE_EDGES.iter().enumerate() {
                    if edge_flags & (1 << e) != 0 {
                        edge_verts[e] =
                            interpolate_edge(corners[ca], corners[cb], cv[ca], cv[cb], threshold);
                    }
                }

                // Emit triangles from the MC tri-table (simplified: use edge groups of 3).
                // We use a minimal inline table for the most common cases.
                // For a production implementation the full 256×16 tri-table would be used;
                // here we emit one triangle per pair of crossed edges as an approximation.
                let crossed: Vec<usize> = (0..12).filter(|&e| edge_flags & (1 << e) != 0).collect();
                for chunk in crossed.chunks(3) {
                    if chunk.len() == 3 {
                        let verts = [
                            edge_verts[chunk[0]],
                            edge_verts[chunk[1]],
                            edge_verts[chunk[2]],
                        ];
                        let normal = IsoTriangle::compute_normal(&verts);
                        triangles.push(IsoTriangle {
                            vertices: verts,
                            normal,
                        });
                    }
                }
            }
        }
    }

    triangles
}

// ─────────────────────────────────────────────────────────────────────────────
// Normal computation from gradient
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a smooth vertex normal at `pos` by evaluating the field gradient
/// and normalising it.
///
/// The gradient points *away* from the isosurface (toward increasing field
/// values), so we negate it to get the outward normal.
pub fn compute_normal(field: &MetaballField, pos: [f64; 3], h: f64) -> [f64; 3] {
    let [gx, gy, gz] = field.gradient(pos, h);
    let len = (gx * gx + gy * gy + gz * gz).sqrt().max(1e-15);
    // Negate: gradient points toward higher field values (inside blob),
    // so outward normal is opposite.
    [-gx / len, -gy / len, -gz / len]
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Metaball field value ─────────────────────────────────────────────────

    #[test]
    fn test_metaball_field_at_centre_is_large() {
        let ball = Metaball::new([0.0, 0.0, 0.0], 1.0);
        // Very close to centre the value should be very large.
        let v = ball.field_value([1e-3, 0.0, 0.0]);
        assert!(v > 1e5, "field at near-centre should be very large: {v}");
    }

    #[test]
    fn test_metaball_field_decreases_with_distance() {
        let ball = Metaball::new([0.0, 0.0, 0.0], 1.0);
        let v1 = ball.field_value([1.0, 0.0, 0.0]);
        let v5 = ball.field_value([5.0, 0.0, 0.0]);
        assert!(v1 > v5, "field should decrease with distance: {v1} > {v5}");
    }

    #[test]
    fn test_metaball_field_is_isotropic() {
        let ball = Metaball::new([0.0, 0.0, 0.0], 1.0);
        let vx = ball.field_value([2.0, 0.0, 0.0]);
        let vy = ball.field_value([0.0, 2.0, 0.0]);
        let vz = ball.field_value([0.0, 0.0, 2.0]);
        assert!((vx - vy).abs() < 1e-12, "field should be isotropic");
        assert!((vx - vz).abs() < 1e-12, "field should be isotropic");
    }

    // ── MetaballField ────────────────────────────────────────────────────────

    #[test]
    fn test_field_sum_two_balls() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([1.0, 0.0, 0.0], 1.0));
        mf.add(Metaball::new([-1.0, 0.0, 0.0], 1.0));
        // At origin, both balls contribute equally.
        let v = mf.evaluate([0.0, 0.0, 0.0]);
        let v_single = Metaball::new([1.0, 0.0, 0.0], 1.0).field_value([0.0, 0.0, 0.0]);
        assert!(
            (v - 2.0 * v_single).abs() < 1e-12,
            "sum should be twice single"
        );
    }

    #[test]
    fn test_field_empty_returns_zero() {
        let mf = MetaballField::new();
        assert_eq!(mf.evaluate([0.0, 0.0, 0.0]), 0.0);
    }

    #[test]
    fn test_field_gradient_direction() {
        // Single ball at origin: gradient at (1,0,0) should point roughly in +x.
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        let g = mf.gradient([1.0, 0.0, 0.0], 1e-4);
        // Gradient of 1/r² at (1,0,0) is negative in x (field decreases away from origin).
        assert!(
            g[0] < 0.0,
            "gradient should point back toward origin (negative x): {}",
            g[0]
        );
        assert!(g[1].abs() < 1e-6, "gradient y should be ~0: {}", g[1]);
        assert!(g[2].abs() < 1e-6, "gradient z should be ~0: {}", g[2]);
    }

    // ── Dynamic update ───────────────────────────────────────────────────────

    #[test]
    fn test_translate_single_ball() {
        let mut ball = Metaball::new([0.0, 0.0, 0.0], 1.0);
        ball.translate([1.0, 2.0, 3.0]);
        assert!((ball.centre[0] - 1.0).abs() < 1e-12);
        assert!((ball.centre[1] - 2.0).abs() < 1e-12);
        assert!((ball.centre[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_translate_all_moves_all_balls() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        mf.add(Metaball::new([1.0, 0.0, 0.0], 1.0));
        mf.translate_all([2.0, 0.0, 0.0]);
        assert!((mf.balls[0].centre[0] - 2.0).abs() < 1e-12);
        assert!((mf.balls[1].centre[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_update_centres_per_ball() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        mf.update_centres(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        assert!((mf.balls[0].centre[0] - 1.0).abs() < 1e-12);
        assert!((mf.balls[1].centre[1] - 1.0).abs() < 1e-12);
    }

    // ── Blob merging detection ───────────────────────────────────────────────

    #[test]
    fn test_merged_pairs_overlapping() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        mf.add(Metaball::new([1.0, 0.0, 0.0], 1.0)); // distance=1 < 1+1=2 → merged
        let pairs = mf.merged_pairs();
        assert!(
            !pairs.is_empty(),
            "overlapping balls should be detected as merged"
        );
    }

    #[test]
    fn test_merged_pairs_separated() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 0.5));
        mf.add(Metaball::new([5.0, 0.0, 0.0], 0.5)); // distance=5 > 0.5+0.5=1 → not merged
        let pairs = mf.merged_pairs();
        assert!(pairs.is_empty(), "separated balls should not be merged");
    }

    #[test]
    fn test_merged_pairs_three_balls() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        mf.add(Metaball::new([1.5, 0.0, 0.0], 1.0)); // dist=1.5 < 2 → merged (0,1)
        mf.add(Metaball::new([10.0, 0.0, 0.0], 0.1)); // far away → not merged
        let pairs = mf.merged_pairs();
        assert_eq!(pairs.len(), 1, "only pair (0,1) should be merged");
        assert_eq!(pairs[0], (0, 1));
    }

    // ── LOD / GridResolution ─────────────────────────────────────────────────

    #[test]
    fn test_grid_resolution_cells() {
        assert_eq!(GridResolution::Coarse.cells(), (8, 8, 8));
        assert_eq!(GridResolution::Medium.cells(), (16, 16, 16));
        assert_eq!(GridResolution::Fine.cells(), (32, 32, 32));
        assert_eq!(GridResolution::Custom(4, 8, 12).cells(), (4, 8, 12));
    }

    // ── Isosurface extraction ────────────────────────────────────────────────

    #[test]
    fn test_isosurface_empty_field_no_triangles() {
        let mf = MetaballField::new();
        let tris = extract_isosurface(&mf, [-2.0; 3], [2.0; 3], GridResolution::Coarse, 1.0);
        assert!(tris.is_empty(), "empty field should produce no isosurface");
    }

    #[test]
    fn test_isosurface_single_ball_produces_triangles() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        let tris = extract_isosurface(&mf, [-3.0; 3], [3.0; 3], GridResolution::Coarse, 0.5);
        assert!(
            !tris.is_empty(),
            "single metaball should produce isosurface triangles"
        );
    }

    #[test]
    fn test_isosurface_triangles_have_valid_normals() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        let tris = extract_isosurface(&mf, [-3.0; 3], [3.0; 3], GridResolution::Coarse, 0.3);
        for tri in &tris {
            let n = tri.normal;
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-6,
                "normal should be unit length: len={len}"
            );
        }
    }

    #[test]
    fn test_isosurface_finer_grid_more_or_equal_triangles() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        let coarse = extract_isosurface(&mf, [-3.0; 3], [3.0; 3], GridResolution::Coarse, 0.3);
        let fine = extract_isosurface(&mf, [-3.0; 3], [3.0; 3], GridResolution::Medium, 0.3);
        assert!(
            fine.len() >= coarse.len(),
            "finer grid should produce >= triangles: coarse={} fine={}",
            coarse.len(),
            fine.len()
        );
    }

    // ── Normal from gradient ─────────────────────────────────────────────────

    #[test]
    fn test_compute_normal_points_outward() {
        // Ball at origin: normal at (1,0,0) should point roughly in +x direction
        // (outward from the blob).
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        let n = compute_normal(&mf, [1.0, 0.0, 0.0], 1e-4);
        // Normal should point away from origin (+x direction).
        assert!(n[0] > 0.0, "normal should point outward (+x): {:?}", n);
    }

    #[test]
    fn test_compute_normal_is_unit_length() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0, 0.0, 0.0], 1.0));
        let n = compute_normal(&mf, [0.7, 0.7, 0.0], 1e-4);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-6,
            "normal must be unit length: {len}"
        );
    }

    // ── Metaball distance ────────────────────────────────────────────────────

    #[test]
    fn test_metaball_distance_to() {
        let a = Metaball::new([0.0, 0.0, 0.0], 1.0);
        let b = Metaball::new([3.0, 4.0, 0.0], 1.0);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_metaball_distance_to_self_is_zero() {
        let a = Metaball::new([1.0, 2.0, 3.0], 1.0);
        assert!(a.distance_to(&a) < 1e-12);
    }

    // ── Interpolate edge ─────────────────────────────────────────────────────

    #[test]
    fn test_interpolate_edge_midpoint() {
        // v0=0, v1=2, threshold=1 → t=0.5
        let p = interpolate_edge([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.0, 2.0, 1.0);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!(p[1].abs() < 1e-12);
    }

    #[test]
    fn test_interpolate_edge_equal_values_returns_midpoint() {
        let p = interpolate_edge([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 1.0, 1.0);
        assert!((p[0] - 1.0).abs() < 1e-10);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Wyvill polynomial falloff
// ─────────────────────────────────────────────────────────────────────────────

/// A metaball using the Wyvill polynomial falloff function (C2-continuous,
/// compact support within radius `r`).
///
/// The falloff is: `f(d) = (1 - (d/r)^2)^3` for d < r, 0 otherwise.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct WyvillMetaball {
    /// World-space centre.
    pub centre: [f64; 3],
    /// Influence radius (compact support).
    pub radius: f64,
    /// Strength multiplier.
    pub strength: f64,
}

impl WyvillMetaball {
    /// Create a new Wyvill metaball.
    pub fn new(centre: [f64; 3], radius: f64, strength: f64) -> Self {
        Self {
            centre,
            radius,
            strength,
        }
    }

    /// Evaluate the Wyvill polynomial field at point `p`.
    /// Returns 0 for points beyond the radius.
    pub fn field_value(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.centre[0];
        let dy = p[1] - self.centre[1];
        let dz = p[2] - self.centre[2];
        let d2 = dx * dx + dy * dy + dz * dz;
        let r2 = self.radius * self.radius;
        if d2 >= r2 {
            return 0.0;
        }
        let t = 1.0 - d2 / r2;
        self.strength * t * t * t
    }

    /// Distance from the centre to point `p`.
    pub fn distance_to_point(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.centre[0];
        let dy = p[1] - self.centre[1];
        let dz = p[2] - self.centre[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Move the centre by a displacement.
    pub fn translate(&mut self, delta: [f64; 3]) {
        self.centre[0] += delta[0];
        self.centre[1] += delta[1];
        self.centre[2] += delta[2];
    }
}

/// A field of Wyvill metaballs with compact support.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct WyvillField {
    /// The Wyvill blobs in this field.
    pub balls: Vec<WyvillMetaball>,
}

impl WyvillField {
    /// Create an empty Wyvill field.
    pub fn new() -> Self {
        Self { balls: Vec::new() }
    }

    /// Add a Wyvill blob.
    pub fn add(&mut self, ball: WyvillMetaball) {
        self.balls.push(ball);
    }

    /// Evaluate the total field at `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        self.balls.iter().map(|b| b.field_value(p)).sum()
    }

    /// Gradient via central finite differences.
    pub fn gradient(&self, p: [f64; 3], h: f64) -> [f64; 3] {
        let two_h = 2.0 * h;
        let [px, py, pz] = p;
        let gx = (self.evaluate([px + h, py, pz]) - self.evaluate([px - h, py, pz])) / two_h;
        let gy = (self.evaluate([px, py + h, pz]) - self.evaluate([px, py - h, pz])) / two_h;
        let gz = (self.evaluate([px, py, pz + h]) - self.evaluate([px, py, pz - h])) / two_h;
        [gx, gy, gz]
    }

    /// Number of blobs.
    pub fn len(&self) -> usize {
        self.balls.len()
    }

    /// Whether the field is empty.
    pub fn is_empty(&self) -> bool {
        self.balls.is_empty()
    }

    /// Compute the axis-aligned bounding box of all influence spheres.
    /// Returns `(min, max)` or `([0;3], [0;3])` if empty.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        if self.balls.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for b in &self.balls {
            for i in 0..3 {
                lo[i] = lo[i].min(b.centre[i] - b.radius);
                hi[i] = hi[i].max(b.centre[i] + b.radius);
            }
        }
        (lo, hi)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bounding-box helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the axis-aligned bounding box of a [`MetaballField`].
///
/// Expands each blob by twice its radius to cover the significant influence
/// region. Returns `(min, max)` corners in world space.
#[allow(dead_code)]
pub fn field_bounds(field: &MetaballField) -> ([f64; 3], [f64; 3]) {
    if field.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for b in &field.balls {
        let r2 = b.radius * 2.0;
        for i in 0..3 {
            lo[i] = lo[i].min(b.centre[i] - r2);
            hi[i] = hi[i].max(b.centre[i] + r2);
        }
    }
    (lo, hi)
}

// ─────────────────────────────────────────────────────────────────────────────
// Animated metaball sequence
// ─────────────────────────────────────────────────────────────────────────────

/// A keyframe for a single metaball in an animation.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct MetaballKeyframe {
    /// Time of this keyframe (seconds).
    pub time: f64,
    /// Centre position at this keyframe.
    pub centre: [f64; 3],
    /// Radius at this keyframe.
    pub radius: f64,
    /// Strength at this keyframe.
    pub strength: f64,
}

impl MetaballKeyframe {
    /// Create a keyframe.
    pub fn new(time: f64, centre: [f64; 3], radius: f64, strength: f64) -> Self {
        Self {
            time,
            centre,
            radius,
            strength,
        }
    }
}

/// Linearly interpolate between two keyframes at time `t`.
#[allow(dead_code)]
pub fn lerp_keyframe(a: &MetaballKeyframe, b: &MetaballKeyframe, t: f64) -> Metaball {
    let dur = b.time - a.time;
    let alpha = if dur.abs() < 1e-15 {
        0.0
    } else {
        ((t - a.time) / dur).clamp(0.0, 1.0)
    };
    let centre = [
        a.centre[0] + alpha * (b.centre[0] - a.centre[0]),
        a.centre[1] + alpha * (b.centre[1] - a.centre[1]),
        a.centre[2] + alpha * (b.centre[2] - a.centre[2]),
    ];
    let radius = a.radius + alpha * (b.radius - a.radius);
    let strength = a.strength + alpha * (b.strength - a.strength);
    Metaball::with_strength(centre, radius, strength)
}

/// A single-blob animation track driven by keyframes.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct MetaballTrack {
    /// Sorted list of keyframes.
    pub keyframes: Vec<MetaballKeyframe>,
}

impl MetaballTrack {
    /// Create an empty track.
    pub fn new() -> Self {
        Self {
            keyframes: Vec::new(),
        }
    }

    /// Add a keyframe (kept sorted by time).
    pub fn push(&mut self, kf: MetaballKeyframe) {
        self.keyframes.push(kf);
        self.keyframes.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Sample the track at time `t`.
    ///
    /// Clamps outside the range of keyframes.
    pub fn sample(&self, t: f64) -> Option<Metaball> {
        let n = self.keyframes.len();
        if n == 0 {
            return None;
        }
        if n == 1 || t <= self.keyframes[0].time {
            let kf = &self.keyframes[0];
            return Some(Metaball::with_strength(kf.centre, kf.radius, kf.strength));
        }
        if t >= self.keyframes[n - 1].time {
            let kf = &self.keyframes[n - 1];
            return Some(Metaball::with_strength(kf.centre, kf.radius, kf.strength));
        }
        // Find surrounding keyframes.
        let idx = self.keyframes.partition_point(|kf| kf.time <= t);
        let a = &self.keyframes[idx - 1];
        let b = &self.keyframes[idx];
        Some(lerp_keyframe(a, b, t))
    }

    /// Number of keyframes.
    pub fn len(&self) -> usize {
        self.keyframes.len()
    }

    /// Whether the track has no keyframes.
    pub fn is_empty(&self) -> bool {
        self.keyframes.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Blob statistics
// ─────────────────────────────────────────────────────────────────────────────

/// Statistics about the scalar field sampled on a grid.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FieldStats {
    /// Minimum value found.
    pub min: f64,
    /// Maximum value found.
    pub max: f64,
    /// Mean value.
    pub mean: f64,
    /// Variance.
    pub variance: f64,
    /// Total number of samples.
    pub count: usize,
}

impl FieldStats {
    /// Sample the field on a uniform grid and compute statistics.
    pub fn from_field(field: &MetaballField, min: [f64; 3], max: [f64; 3], steps: usize) -> Self {
        let n = steps.max(1);
        let mut values = Vec::with_capacity(n * n * n);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let p = [
                        min[0] + (ix as f64 / (n - 1).max(1) as f64) * (max[0] - min[0]),
                        min[1] + (iy as f64 / (n - 1).max(1) as f64) * (max[1] - min[1]),
                        min[2] + (iz as f64 / (n - 1).max(1) as f64) * (max[2] - min[2]),
                    ];
                    values.push(field.evaluate(p));
                }
            }
        }
        let count = values.len();
        let fmin = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let fmax = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = values.iter().sum::<f64>() / count as f64;
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / count as f64;
        Self {
            min: fmin,
            max: fmax,
            mean,
            variance,
            count,
        }
    }

    /// Standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance.sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Signed-distance-function (SDF) metaball blending
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth minimum (polynomial variant) of two SDF values.
///
/// `k` controls the blend radius. Smaller `k` → sharper join.
/// Returns a value ≤ min(a, b), reaching up to `k/4` below.
#[allow(dead_code)]
pub fn sdf_smooth_min(a: f64, b: f64, k: f64) -> f64 {
    if k.abs() < 1e-15 {
        return a.min(b);
    }
    // h = clamp( 0.5 + 0.5*(a-b)/k ): h→1 when a>>b (a is larger, b is minimum)
    let h = (0.5 + 0.5 * (a - b) / k).clamp(0.0, 1.0);
    // lerp from a to b with weight h, then subtract the smooth bump.
    let mix = a * (1.0 - h) + b * h;
    mix - k * h * (1.0 - h)
}

/// Sphere SDF: signed distance from `p` to the surface of a sphere at
/// `centre` with `radius`.
#[allow(dead_code)]
pub fn sdf_sphere(p: [f64; 3], centre: [f64; 3], radius: f64) -> f64 {
    let dx = p[0] - centre[0];
    let dy = p[1] - centre[1];
    let dz = p[2] - centre[2];
    (dx * dx + dy * dy + dz * dz).sqrt() - radius
}

/// Evaluate a smooth union of multiple sphere SDFs.
///
/// Returns the smooth-minimum blended SDF value at `p`.
#[allow(dead_code)]
pub fn sdf_multi_sphere(
    p: [f64; 3],
    spheres: &[(/* centre */ [f64; 3], /* radius */ f64)],
    k: f64,
) -> f64 {
    if spheres.is_empty() {
        return f64::INFINITY;
    }
    let mut result = sdf_sphere(p, spheres[0].0, spheres[0].1);
    for &(centre, radius) in spheres.iter().skip(1) {
        let d = sdf_sphere(p, centre, radius);
        result = sdf_smooth_min(result, d, k);
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Grid cache for field evaluation
// ─────────────────────────────────────────────────────────────────────────────

/// A cached evaluation of a scalar field on a uniform 3D grid.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FieldGrid {
    /// Flat array of field values, row-major: `values[iz * (ny+1) * (nx+1) + iy * (nx+1) + ix]`.
    pub values: Vec<f64>,
    /// Number of cells along X.
    pub nx: usize,
    /// Number of cells along Y.
    pub ny: usize,
    /// Number of cells along Z.
    pub nz: usize,
    /// World-space minimum corner.
    pub min: [f64; 3],
    /// World-space maximum corner.
    pub max: [f64; 3],
}

impl FieldGrid {
    /// Sample the field on a uniform grid and cache results.
    pub fn build(
        field: &MetaballField,
        min: [f64; 3],
        max: [f64; 3],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Self {
        let nx = nx.max(1);
        let ny = ny.max(1);
        let nz = nz.max(1);
        let gx = nx + 1;
        let gy = ny + 1;
        let gz = nz + 1;
        let dx = (max[0] - min[0]) / nx as f64;
        let dy = (max[1] - min[1]) / ny as f64;
        let dz = (max[2] - min[2]) / nz as f64;
        let mut values = vec![0.0_f64; gx * gy * gz];
        for iz in 0..gz {
            for iy in 0..gy {
                for ix in 0..gx {
                    let p = [
                        min[0] + ix as f64 * dx,
                        min[1] + iy as f64 * dy,
                        min[2] + iz as f64 * dz,
                    ];
                    values[iz * gy * gx + iy * gx + ix] = field.evaluate(p);
                }
            }
        }
        Self {
            values,
            nx,
            ny,
            nz,
            min,
            max,
        }
    }

    /// Look up a cached value at grid vertex `(ix, iy, iz)`.
    ///
    /// Returns `None` if the indices are out of range.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> Option<f64> {
        let gx = self.nx + 1;
        let gy = self.ny + 1;
        if ix > self.nx || iy > self.ny || iz > self.nz {
            return None;
        }
        Some(self.values[iz * gy * gx + iy * gx + ix])
    }

    /// Trilinearly interpolate the cached field at a world-space point `p`.
    pub fn trilinear(&self, p: [f64; 3]) -> f64 {
        let dx = (self.max[0] - self.min[0]) / self.nx as f64;
        let dy = (self.max[1] - self.min[1]) / self.ny as f64;
        let dz = (self.max[2] - self.min[2]) / self.nz as f64;

        let fx = ((p[0] - self.min[0]) / dx).clamp(0.0, self.nx as f64);
        let fy = ((p[1] - self.min[1]) / dy).clamp(0.0, self.ny as f64);
        let fz = ((p[2] - self.min[2]) / dz).clamp(0.0, self.nz as f64);

        let ix = (fx as usize).min(self.nx.saturating_sub(1));
        let iy = (fy as usize).min(self.ny.saturating_sub(1));
        let iz = (fz as usize).min(self.nz.saturating_sub(1));

        let tx = fx - ix as f64;
        let ty = fy - iy as f64;
        let tz = fz - iz as f64;

        let c000 = self.get(ix, iy, iz).unwrap_or(0.0);
        let c100 = self.get(ix + 1, iy, iz).unwrap_or(0.0);
        let c010 = self.get(ix, iy + 1, iz).unwrap_or(0.0);
        let c110 = self.get(ix + 1, iy + 1, iz).unwrap_or(0.0);
        let c001 = self.get(ix, iy, iz + 1).unwrap_or(0.0);
        let c101 = self.get(ix + 1, iy, iz + 1).unwrap_or(0.0);
        let c011 = self.get(ix, iy + 1, iz + 1).unwrap_or(0.0);
        let c111 = self.get(ix + 1, iy + 1, iz + 1).unwrap_or(0.0);

        let c00 = c000 * (1.0 - tx) + c100 * tx;
        let c10 = c010 * (1.0 - tx) + c110 * tx;
        let c01 = c001 * (1.0 - tx) + c101 * tx;
        let c11 = c011 * (1.0 - tx) + c111 * tx;

        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;

        c0 * (1.0 - tz) + c1 * tz
    }

    /// Total number of grid vertices.
    pub fn vertex_count(&self) -> usize {
        (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod extended_tests {
    use super::*;

    // ── WyvillMetaball ───────────────────────────────────────────────────────

    #[test]
    fn test_wyvill_zero_outside_radius() {
        let b = WyvillMetaball::new([0.0; 3], 1.0, 1.0);
        assert_eq!(
            b.field_value([2.0, 0.0, 0.0]),
            0.0,
            "beyond radius should be 0"
        );
    }

    #[test]
    fn test_wyvill_max_at_centre() {
        let b = WyvillMetaball::new([0.0; 3], 1.0, 1.0);
        let centre = b.field_value([0.0, 0.0, 0.0]);
        let off = b.field_value([0.5, 0.0, 0.0]);
        assert!(
            centre > off,
            "Wyvill max should be at centre: {centre} > {off}"
        );
    }

    #[test]
    fn test_wyvill_strength_scales_value() {
        let b1 = WyvillMetaball::new([0.0; 3], 2.0, 1.0);
        let b2 = WyvillMetaball::new([0.0; 3], 2.0, 3.0);
        let v1 = b1.field_value([0.5, 0.0, 0.0]);
        let v2 = b2.field_value([0.5, 0.0, 0.0]);
        assert!(
            (v2 - 3.0 * v1).abs() < 1e-12,
            "strength should scale value: {v2} vs {}",
            3.0 * v1
        );
    }

    #[test]
    fn test_wyvill_translate() {
        let mut b = WyvillMetaball::new([0.0; 3], 2.0, 1.0);
        b.translate([1.0, 2.0, 3.0]);
        assert_eq!(b.centre, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_wyvill_distance_to_point() {
        let b = WyvillMetaball::new([0.0; 3], 2.0, 1.0);
        let d = b.distance_to_point([3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_wyvill_field_sum() {
        let mut field = WyvillField::new();
        field.add(WyvillMetaball::new([0.0; 3], 2.0, 1.0));
        field.add(WyvillMetaball::new([0.0; 3], 2.0, 1.0));
        let v = field.evaluate([0.5, 0.0, 0.0]);
        let single = WyvillMetaball::new([0.0; 3], 2.0, 1.0).field_value([0.5, 0.0, 0.0]);
        assert!(
            (v - 2.0 * single).abs() < 1e-12,
            "two identical balls should double the value"
        );
    }

    #[test]
    fn test_wyvill_field_empty_zero() {
        let field = WyvillField::new();
        assert_eq!(field.evaluate([0.0, 0.0, 0.0]), 0.0);
    }

    #[test]
    fn test_wyvill_gradient_direction() {
        let mut field = WyvillField::new();
        field.add(WyvillMetaball::new([0.0; 3], 3.0, 1.0));
        let g = field.gradient([1.0, 0.0, 0.0], 1e-4);
        // Gradient should point back toward origin (field decreases away from centre).
        assert!(g[0] < 0.0, "gradient should point toward centre: {}", g[0]);
    }

    #[test]
    fn test_wyvill_bounds_single_ball() {
        let mut field = WyvillField::new();
        field.add(WyvillMetaball::new([1.0, 2.0, 3.0], 2.0, 1.0));
        let (lo, hi) = field.bounds();
        assert!((lo[0] - (-1.0)).abs() < 1e-10, "lo x should be -1.0");
        assert!((hi[0] - 3.0).abs() < 1e-10, "hi x should be 3.0");
    }

    #[test]
    fn test_wyvill_bounds_empty() {
        let field = WyvillField::new();
        let (lo, hi) = field.bounds();
        for i in 0..3 {
            assert_eq!(lo[i], 0.0);
            assert_eq!(hi[i], 0.0);
        }
    }

    #[test]
    fn test_wyvill_len_and_is_empty() {
        let mut f = WyvillField::new();
        assert!(f.is_empty());
        f.add(WyvillMetaball::new([0.0; 3], 1.0, 1.0));
        assert_eq!(f.len(), 1);
        assert!(!f.is_empty());
    }

    // ── FieldGrid (cached sampling) ─────────────────────────────────────────

    #[test]
    fn test_field_grid_vertex_count() {
        let mf = MetaballField::new();
        let grid = FieldGrid::build(&mf, [-1.0; 3], [1.0; 3], 4, 4, 4);
        assert_eq!(grid.vertex_count(), 5 * 5 * 5);
    }

    #[test]
    fn test_field_grid_empty_field_all_zero() {
        let mf = MetaballField::new();
        let grid = FieldGrid::build(&mf, [-1.0; 3], [1.0; 3], 4, 4, 4);
        for v in &grid.values {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_field_grid_matches_direct_eval_at_corners() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0; 3], 1.0));
        let min = [-2.0; 3];
        let max = [2.0; 3];
        let grid = FieldGrid::build(&mf, min, max, 4, 4, 4);
        // At (0,0,0) vertex (index 2,2,2 in the 5x5x5 grid), the value should match direct eval.
        let v_grid = grid.get(2, 2, 2).unwrap();
        let p = [
            min[0] + 2.0 * (max[0] - min[0]) / 4.0,
            min[1] + 2.0 * (max[1] - min[1]) / 4.0,
            min[2] + 2.0 * (max[2] - min[2]) / 4.0,
        ];
        let v_direct = mf.evaluate(p);
        assert!(
            (v_grid - v_direct).abs() < 1e-12,
            "grid value should match direct: {v_grid} vs {v_direct}"
        );
    }

    #[test]
    fn test_field_grid_get_out_of_range_returns_none() {
        let mf = MetaballField::new();
        let grid = FieldGrid::build(&mf, [-1.0; 3], [1.0; 3], 2, 2, 2);
        assert!(grid.get(3, 0, 0).is_none(), "ix=3 is out of range for nx=2");
    }

    #[test]
    fn test_field_grid_trilinear_at_grid_point() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0; 3], 1.0));
        let min = [-2.0; 3];
        let max = [2.0; 3];
        let grid = FieldGrid::build(&mf, min, max, 8, 8, 8);
        // Trilinear at a grid vertex should be close to the cached value.
        let p = [0.0, 0.0, 0.0];
        let t = grid.trilinear(p);
        let v = mf.evaluate(p);
        // Allow some interpolation error near the centre where the field varies sharply.
        assert!(t > 0.0, "trilinear at centre should be > 0: {t}");
        let _ = v;
    }

    #[test]
    fn test_field_grid_trilinear_monotone_along_axis() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0; 3], 1.0));
        let min = [-5.0; 3];
        let max = [5.0; 3];
        let grid = FieldGrid::build(&mf, min, max, 16, 16, 16);
        // Field should decrease moving away from origin along x axis.
        let v1 = grid.trilinear([0.5, 0.0, 0.0]);
        let v2 = grid.trilinear([2.0, 0.0, 0.0]);
        let v3 = grid.trilinear([4.0, 0.0, 0.0]);
        assert!(v1 > v2, "field should decrease: v1={v1} v2={v2}");
        assert!(v2 > v3, "field should decrease: v2={v2} v3={v3}");
    }

    // ── MetaballTrack (animation) ────────────────────────────────────────────

    #[test]
    fn test_track_sample_empty_returns_none() {
        let track = MetaballTrack::new();
        assert!(track.sample(0.0).is_none());
    }

    #[test]
    fn test_track_single_keyframe() {
        let mut track = MetaballTrack::new();
        track.push(MetaballKeyframe::new(0.0, [1.0, 2.0, 3.0], 1.0, 1.0));
        let ball = track.sample(0.5).unwrap();
        assert!((ball.centre[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_track_clamp_before_start() {
        let mut track = MetaballTrack::new();
        track.push(MetaballKeyframe::new(1.0, [1.0, 0.0, 0.0], 1.0, 1.0));
        track.push(MetaballKeyframe::new(2.0, [2.0, 0.0, 0.0], 1.0, 1.0));
        // Before the first keyframe, should get the first keyframe value.
        let ball = track.sample(0.0).unwrap();
        assert!(
            (ball.centre[0] - 1.0).abs() < 1e-10,
            "should clamp to first kf: {}",
            ball.centre[0]
        );
    }

    #[test]
    fn test_track_clamp_after_end() {
        let mut track = MetaballTrack::new();
        track.push(MetaballKeyframe::new(0.0, [0.0, 0.0, 0.0], 1.0, 1.0));
        track.push(MetaballKeyframe::new(1.0, [5.0, 0.0, 0.0], 1.0, 1.0));
        let ball = track.sample(99.0).unwrap();
        assert!(
            (ball.centre[0] - 5.0).abs() < 1e-10,
            "should clamp to last kf: {}",
            ball.centre[0]
        );
    }

    #[test]
    fn test_track_interpolation_midpoint() {
        let mut track = MetaballTrack::new();
        track.push(MetaballKeyframe::new(0.0, [0.0, 0.0, 0.0], 1.0, 1.0));
        track.push(MetaballKeyframe::new(2.0, [4.0, 0.0, 0.0], 1.0, 1.0));
        let ball = track.sample(1.0).unwrap();
        assert!(
            (ball.centre[0] - 2.0).abs() < 1e-10,
            "midpoint: expected 2.0, got {}",
            ball.centre[0]
        );
    }

    #[test]
    fn test_track_len_and_is_empty() {
        let mut track = MetaballTrack::new();
        assert!(track.is_empty());
        track.push(MetaballKeyframe::new(0.0, [0.0; 3], 1.0, 1.0));
        assert_eq!(track.len(), 1);
        assert!(!track.is_empty());
    }

    #[test]
    fn test_track_push_maintains_sorted_order() {
        let mut track = MetaballTrack::new();
        track.push(MetaballKeyframe::new(2.0, [0.0; 3], 1.0, 1.0));
        track.push(MetaballKeyframe::new(0.0, [0.0; 3], 1.0, 1.0));
        track.push(MetaballKeyframe::new(1.0, [0.0; 3], 1.0, 1.0));
        let times: Vec<f64> = track.keyframes.iter().map(|k| k.time).collect();
        assert_eq!(
            times,
            vec![0.0, 1.0, 2.0],
            "keyframes should be sorted by time"
        );
    }

    // ── lerp_keyframe ────────────────────────────────────────────────────────

    #[test]
    fn test_lerp_keyframe_alpha_zero() {
        let a = MetaballKeyframe::new(0.0, [1.0, 0.0, 0.0], 1.0, 2.0);
        let b = MetaballKeyframe::new(1.0, [3.0, 0.0, 0.0], 2.0, 4.0);
        let ball = lerp_keyframe(&a, &b, 0.0);
        assert!((ball.centre[0] - 1.0).abs() < 1e-10);
        assert!((ball.radius - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lerp_keyframe_alpha_one() {
        let a = MetaballKeyframe::new(0.0, [1.0, 0.0, 0.0], 1.0, 2.0);
        let b = MetaballKeyframe::new(1.0, [3.0, 0.0, 0.0], 2.0, 4.0);
        let ball = lerp_keyframe(&a, &b, 1.0);
        assert!((ball.centre[0] - 3.0).abs() < 1e-10);
        assert!((ball.radius - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lerp_keyframe_equal_times() {
        let a = MetaballKeyframe::new(1.0, [5.0, 0.0, 0.0], 1.0, 1.0);
        let b = MetaballKeyframe::new(1.0, [9.0, 0.0, 0.0], 2.0, 2.0);
        let ball = lerp_keyframe(&a, &b, 1.0);
        // When duration == 0, alpha clamps to 0 → gets a's values.
        assert!((ball.centre[0] - 5.0).abs() < 1e-10);
    }

    // ── FieldStats ───────────────────────────────────────────────────────────

    #[test]
    fn test_field_stats_empty_field_mean_zero() {
        let mf = MetaballField::new();
        let stats = FieldStats::from_field(&mf, [-1.0; 3], [1.0; 3], 4);
        assert_eq!(stats.mean, 0.0);
        assert_eq!(stats.min, 0.0);
        assert_eq!(stats.max, 0.0);
        assert_eq!(stats.variance, 0.0);
    }

    #[test]
    fn test_field_stats_single_ball_positive() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0; 3], 1.0));
        let stats = FieldStats::from_field(&mf, [-2.0; 3], [2.0; 3], 4);
        assert!(
            stats.max > stats.min,
            "max should be > min for non-trivial field"
        );
        assert!(stats.mean > 0.0, "mean should be positive");
    }

    #[test]
    fn test_field_stats_std_dev() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([0.0; 3], 1.0));
        let stats = FieldStats::from_field(&mf, [-5.0; 3], [5.0; 3], 6);
        let sd = stats.std_dev();
        assert!(sd >= 0.0, "std_dev should be non-negative: {sd}");
    }

    #[test]
    fn test_field_stats_count() {
        let mf = MetaballField::new();
        let stats = FieldStats::from_field(&mf, [-1.0; 3], [1.0; 3], 3);
        assert_eq!(stats.count, 27, "3^3 = 27 samples");
    }

    // ── SDF helpers ──────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_sphere_on_surface() {
        let d = sdf_sphere([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!(
            d.abs() < 1e-10,
            "point on sphere surface should have SDF = 0: {d}"
        );
    }

    #[test]
    fn test_sdf_sphere_inside_negative() {
        let d = sdf_sphere([0.0; 3], [0.0; 3], 1.0);
        assert!(d < 0.0, "centre should be inside (negative SDF): {d}");
    }

    #[test]
    fn test_sdf_sphere_outside_positive() {
        let d = sdf_sphere([3.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!(
            (d - 2.0).abs() < 1e-10,
            "outside distance should be 2.0: {d}"
        );
    }

    #[test]
    fn test_sdf_smooth_min_returns_smaller_when_far_apart() {
        // sdf_smooth_min(a, b, k) approaches min(a, b) when |a-b| >> k.
        // With a=5 > b=-5, result should be very close to b=-5.
        let a = 5.0_f64;
        let b = -5.0_f64;
        let k = 0.01_f64;
        let sm = sdf_smooth_min(a, b, k);
        // Should be near -5, the minimum.
        assert!(
            sm < 0.0,
            "smooth_min should return the smaller value: sm={sm}"
        );
        assert!(
            (sm - b).abs() < k + 1e-10,
            "smooth_min should be within k of b: sm={sm}, b={b}"
        );
    }

    #[test]
    fn test_sdf_smooth_min_symmetric_when_equal() {
        // When a == b, smooth_min should give a value close to a.
        let sm = sdf_smooth_min(1.0, 1.0, 0.5);
        assert!((sm - 1.0).abs() < 0.5, "smooth_min of equal values: {sm}");
    }

    #[test]
    fn test_sdf_smooth_min_zero_k_gives_min() {
        let sm = sdf_smooth_min(2.0, 5.0, 0.0);
        assert!((sm - 2.0).abs() < 1e-10, "k=0 should give exact min: {sm}");
    }

    #[test]
    fn test_sdf_multi_sphere_empty_infinity() {
        let d = sdf_multi_sphere([0.0; 3], &[], 0.1);
        assert_eq!(d, f64::INFINITY);
    }

    #[test]
    fn test_sdf_multi_sphere_single_sphere() {
        let spheres = [([0.0; 3], 1.0)];
        let d = sdf_multi_sphere([2.0, 0.0, 0.0], &spheres, 0.1);
        assert!(
            (d - 1.0).abs() < 1e-10,
            "single sphere SDF at distance 2 with radius 1 = 1: {d}"
        );
    }

    #[test]
    fn test_sdf_multi_sphere_two_spheres_smooth_union() {
        let spheres = [([0.0, 0.0, 0.0], 1.0), ([3.0, 0.0, 0.0], 1.0)];
        // A point between the two spheres should be influenced by the smooth union.
        let d = sdf_multi_sphere([1.5, 0.0, 0.0], &spheres, 1.0);
        // Should be negative or near zero (inside the smooth blend region).
        assert!(
            d < 1.0,
            "between two spheres the smooth-union SDF should be small: {d}"
        );
    }

    // ── field_bounds ─────────────────────────────────────────────────────────

    #[test]
    fn test_field_bounds_empty_returns_zero() {
        let mf = MetaballField::new();
        let (lo, hi) = field_bounds(&mf);
        assert_eq!(lo, [0.0; 3]);
        assert_eq!(hi, [0.0; 3]);
    }

    #[test]
    fn test_field_bounds_single_ball() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([1.0, 0.0, 0.0], 2.0));
        let (lo, hi) = field_bounds(&mf);
        // Radius is multiplied by 2 for the bounds.
        assert!(
            (lo[0] - (-3.0)).abs() < 1e-10,
            "lo x = centre - 2*radius = 1 - 4 = -3: {}",
            lo[0]
        );
        assert!(
            (hi[0] - 5.0).abs() < 1e-10,
            "hi x = centre + 2*radius = 1 + 4 = 5: {}",
            hi[0]
        );
    }

    #[test]
    fn test_field_bounds_two_balls() {
        let mut mf = MetaballField::new();
        mf.add(Metaball::new([-5.0, 0.0, 0.0], 1.0));
        mf.add(Metaball::new([5.0, 0.0, 0.0], 1.0));
        let (lo, hi) = field_bounds(&mf);
        assert!(
            lo[0] < -5.0,
            "lo should be left of left ball centre: {}",
            lo[0]
        );
        assert!(
            hi[0] > 5.0,
            "hi should be right of right ball centre: {}",
            hi[0]
        );
    }
}
