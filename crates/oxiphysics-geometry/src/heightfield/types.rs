//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::{Real, Vec3};

use super::functions::{heightfield_tessellate, ray_aabb_xz, ray_triangle, tri_area};

/// A grid-based height field for terrain representation.
///
/// Heights are stored in row-major order. The field spans from
/// (0,0) to (scale_x * (cols-1), scale_z * (rows-1)) in the XZ plane.
#[derive(Debug, Clone)]
pub struct HeightField {
    /// Height values in row-major order.
    pub heights: Vec<Real>,
    /// Number of rows (Z direction).
    pub rows: usize,
    /// Number of columns (X direction).
    pub cols: usize,
    /// Spacing in the X direction.
    pub scale_x: Real,
    /// Spacing in the Z direction.
    pub scale_z: Real,
}
impl HeightField {
    /// Create a new height field.
    pub fn new(heights: Vec<Real>, rows: usize, cols: usize, scale_x: Real, scale_z: Real) -> Self {
        assert_eq!(heights.len(), rows * cols);
        Self {
            heights,
            rows,
            cols,
            scale_x,
            scale_z,
        }
    }
    /// Get the height at grid position (row, col).
    pub fn height_at(&self, row: usize, col: usize) -> Real {
        self.heights[row * self.cols + col]
    }
    /// Create a height field from a function.
    ///
    /// The function `f(col, row)` returns the height at each grid point.
    pub fn from_fn(
        cols: usize,
        rows: usize,
        cell_size: Real,
        f: impl Fn(usize, usize) -> Real,
    ) -> Self {
        let mut heights = Vec::with_capacity(rows * cols);
        for row in 0..rows {
            for col in 0..cols {
                heights.push(f(col, row));
            }
        }
        Self {
            heights,
            rows,
            cols,
            scale_x: cell_size,
            scale_z: cell_size,
        }
    }
    /// Bilinear interpolation of height at normalized coordinates (u, v).
    ///
    /// `u` and `v` are in `[0, 1]`, mapping to the full extent of the grid.
    pub fn height_at_uv(&self, u: Real, v: Real) -> Real {
        let u = u.clamp(0.0, 1.0);
        let v = v.clamp(0.0, 1.0);
        let fx = u * (self.cols - 1) as Real;
        let fz = v * (self.rows - 1) as Real;
        let col0 = (fx as usize).min(self.cols - 2);
        let row0 = (fz as usize).min(self.rows - 2);
        let col1 = col0 + 1;
        let row1 = row0 + 1;
        let tx = fx - col0 as Real;
        let tz = fz - row0 as Real;
        let h00 = self.height_at(row0, col0);
        let h10 = self.height_at(row0, col1);
        let h01 = self.height_at(row1, col0);
        let h11 = self.height_at(row1, col1);
        let h0 = h00 * (1.0 - tx) + h10 * tx;
        let h1 = h01 * (1.0 - tx) + h11 * tx;
        h0 * (1.0 - tz) + h1 * tz
    }
    /// Compute the surface normal at grid point (col, row) via finite differences.
    pub fn normal_at_grid(&self, col: usize, row: usize) -> [Real; 3] {
        let dh_dx = if col == 0 {
            (self.height_at(row, col + 1) - self.height_at(row, col)) / self.scale_x
        } else if col == self.cols - 1 {
            (self.height_at(row, col) - self.height_at(row, col - 1)) / self.scale_x
        } else {
            (self.height_at(row, col + 1) - self.height_at(row, col - 1)) / (2.0 * self.scale_x)
        };
        let dh_dz = if row == 0 {
            (self.height_at(row + 1, col) - self.height_at(row, col)) / self.scale_z
        } else if row == self.rows - 1 {
            (self.height_at(row, col) - self.height_at(row - 1, col)) / self.scale_z
        } else {
            (self.height_at(row + 1, col) - self.height_at(row - 1, col)) / (2.0 * self.scale_z)
        };
        let nx = -dh_dx;
        let ny = 1.0;
        let nz = -dh_dz;
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        [nx / len, ny / len, nz / len]
    }
    /// Tessellate the height field into a triangle mesh.
    ///
    /// Returns (vertices, triangles) where each vertex is `[x, y, z]`
    /// and each triangle is 3 vertex indices.
    pub fn to_triangle_mesh(&self) -> (Vec<[Real; 3]>, Vec<[usize; 3]>) {
        let mut vertices = Vec::with_capacity(self.rows * self.cols);
        for row in 0..self.rows {
            for col in 0..self.cols {
                vertices.push([
                    col as Real * self.scale_x,
                    self.height_at(row, col),
                    row as Real * self.scale_z,
                ]);
            }
        }
        let mut triangles = Vec::with_capacity((self.rows - 1) * (self.cols - 1) * 2);
        for row in 0..(self.rows - 1) {
            for col in 0..(self.cols - 1) {
                let i00 = row * self.cols + col;
                let i10 = row * self.cols + col + 1;
                let i01 = (row + 1) * self.cols + col;
                let i11 = (row + 1) * self.cols + col + 1;
                triangles.push([i00, i10, i11]);
                triangles.push([i00, i11, i01]);
            }
        }
        (vertices, triangles)
    }
    /// Return the minimum height in the field.
    pub fn min_height(&self) -> Real {
        self.height_bounds().0
    }
    /// Return the maximum height in the field.
    pub fn max_height(&self) -> Real {
        self.height_bounds().1
    }
    /// Compute the total surface area by summing triangle areas from tessellation.
    pub fn surface_area(&self) -> Real {
        if self.rows < 2 || self.cols < 2 {
            return 0.0;
        }
        let mut area = 0.0;
        for row in 0..(self.rows - 1) {
            for col in 0..(self.cols - 1) {
                let v00 = [
                    col as Real * self.scale_x,
                    self.height_at(row, col),
                    row as Real * self.scale_z,
                ];
                let v10 = [
                    (col + 1) as Real * self.scale_x,
                    self.height_at(row, col + 1),
                    row as Real * self.scale_z,
                ];
                let v01 = [
                    col as Real * self.scale_x,
                    self.height_at(row + 1, col),
                    (row + 1) as Real * self.scale_z,
                ];
                let v11 = [
                    (col + 1) as Real * self.scale_x,
                    self.height_at(row + 1, col + 1),
                    (row + 1) as Real * self.scale_z,
                ];
                area += tri_area(&v00, &v10, &v11);
                area += tri_area(&v00, &v11, &v01);
            }
        }
        area
    }
    /// Apply Laplacian smoothing to the height field.
    ///
    /// Each interior height is replaced by the average of its 4-connected neighbors.
    /// Boundary heights are unchanged.
    pub fn smooth(&mut self, iterations: usize) {
        for _ in 0..iterations {
            let mut new_heights = self.heights.clone();
            for row in 1..(self.rows - 1) {
                for col in 1..(self.cols - 1) {
                    let up = self.height_at(row - 1, col);
                    let down = self.height_at(row + 1, col);
                    let left = self.height_at(row, col - 1);
                    let right = self.height_at(row, col + 1);
                    new_heights[row * self.cols + col] = (up + down + left + right) / 4.0;
                }
            }
            self.heights = new_heights;
        }
    }
    /// Ray cast using grid DDA traversal and per-cell triangle intersection.
    ///
    /// Returns `Some((toi, [nx, ny, nz]))` for the closest hit within `max_toi`,
    /// or `None` if no intersection.
    pub fn ray_cast_grid(
        &self,
        origin: [Real; 3],
        dir: [Real; 3],
        max_toi: Real,
    ) -> Option<(Real, [Real; 3])> {
        if self.rows < 2 || self.cols < 2 {
            return None;
        }
        let grid_w = (self.cols - 1) as Real * self.scale_x;
        let grid_h = (self.rows - 1) as Real * self.scale_z;
        let (t_enter, t_exit) = ray_aabb_xz(
            origin[0], origin[2], dir[0], dir[2], grid_w, grid_h, max_toi,
        )?;
        let eps = 1e-6;
        let mut t = t_enter.max(0.0);
        let step = self.scale_x.min(self.scale_z) * 0.5;
        let mut best: Option<(Real, [Real; 3])> = None;
        while t <= t_exit.min(max_toi) {
            let px = origin[0] + dir[0] * t;
            let pz = origin[2] + dir[2] * t;
            let col = ((px / self.scale_x).floor() as isize)
                .max(0)
                .min((self.cols - 2) as isize) as usize;
            let row = ((pz / self.scale_z).floor() as isize)
                .max(0)
                .min((self.rows - 2) as isize) as usize;
            let v00 = Vec3::new(
                col as Real * self.scale_x,
                self.height_at(row, col),
                row as Real * self.scale_z,
            );
            let v10 = Vec3::new(
                (col + 1) as Real * self.scale_x,
                self.height_at(row, col + 1),
                row as Real * self.scale_z,
            );
            let v01 = Vec3::new(
                col as Real * self.scale_x,
                self.height_at(row + 1, col),
                (row + 1) as Real * self.scale_z,
            );
            let v11 = Vec3::new(
                (col + 1) as Real * self.scale_x,
                self.height_at(row + 1, col + 1),
                (row + 1) as Real * self.scale_z,
            );
            let o = Vec3::new(origin[0], origin[1], origin[2]);
            let d = Vec3::new(dir[0], dir[1], dir[2]);
            for (a, b, c) in [(&v00, &v10, &v11), (&v00, &v11, &v01)] {
                if let Some(hit) = ray_triangle(&o, &d, max_toi, a, b, c)
                    && best.as_ref().is_none_or(|(bt, _)| hit.toi < *bt)
                {
                    let dot = hit.normal.dot(&d);
                    let n = if dot > 0.0 {
                        [-hit.normal.x, -hit.normal.y, -hit.normal.z]
                    } else {
                        [hit.normal.x, hit.normal.y, hit.normal.z]
                    };
                    best = Some((hit.toi, n));
                }
            }
            if best.is_some() {
                return best;
            }
            t += step + eps;
        }
        best
    }
    /// Get min and max heights.
    pub(super) fn height_bounds(&self) -> (Real, Real) {
        let mut min = Real::INFINITY;
        let mut max = Real::NEG_INFINITY;
        for &h in &self.heights {
            if h < min {
                min = h;
            }
            if h > max {
                max = h;
            }
        }
        (min, max)
    }
}
impl HeightField {
    /// Generate a downsampled (LOD) version of the height field.
    ///
    /// `factor` must be ≥ 1.  A factor of 2 halves the resolution by averaging
    /// 2×2 blocks of cells.  If the grid is too small the original is returned.
    pub fn lod_downsample(&self, factor: usize) -> HeightField {
        if factor <= 1 {
            return self.clone();
        }
        let new_cols = self.cols.div_ceil(factor);
        let new_rows = self.rows.div_ceil(factor);
        if new_cols < 2 || new_rows < 2 {
            return self.clone();
        }
        let mut heights = Vec::with_capacity(new_rows * new_cols);
        for r in 0..new_rows {
            for c in 0..new_cols {
                let r0 = r * factor;
                let c0 = c * factor;
                let r1 = (r0 + factor).min(self.rows);
                let c1 = (c0 + factor).min(self.cols);
                let mut sum = 0.0;
                let mut count = 0usize;
                for rr in r0..r1 {
                    for cc in c0..c1 {
                        sum += self.height_at(rr, cc);
                        count += 1;
                    }
                }
                heights.push(if count > 0 { sum / count as Real } else { 0.0 });
            }
        }
        HeightField::new(
            heights,
            new_rows,
            new_cols,
            self.scale_x * factor as Real,
            self.scale_z * factor as Real,
        )
    }
    /// Generate multiple LOD levels.
    ///
    /// Returns a Vec where index 0 is the original, index 1 is 2× downsampled,
    /// index 2 is 4× downsampled, etc., until the grid is too small.
    pub fn lod_pyramid(&self, levels: usize) -> Vec<HeightField> {
        let mut result = vec![self.clone()];
        let mut factor = 2;
        for _ in 1..levels {
            let prev = result.last().expect("collection should not be empty");
            if prev.rows < 4 || prev.cols < 4 {
                break;
            }
            result.push(prev.lod_downsample(2));
            factor *= 2;
        }
        let _ = factor;
        result
    }
    /// Compute the normal at an arbitrary world-space XZ position via bilinear
    /// interpolation of the four surrounding grid normals.
    pub fn normal_at_world(&self, x: Real, z: Real) -> [Real; 3] {
        let fx = (x / self.scale_x).clamp(0.0, (self.cols - 1) as Real);
        let fz = (z / self.scale_z).clamp(0.0, (self.rows - 1) as Real);
        let c0 = (fx as usize).min(self.cols - 1);
        let r0 = (fz as usize).min(self.rows - 1);
        let c1 = (c0 + 1).min(self.cols - 1);
        let r1 = (r0 + 1).min(self.rows - 1);
        let tx = fx - c0 as Real;
        let tz = fz - r0 as Real;
        let n00 = self.normal_at_grid(c0, r0);
        let n10 = self.normal_at_grid(c1, r0);
        let n01 = self.normal_at_grid(c0, r1);
        let n11 = self.normal_at_grid(c1, r1);
        let mut n = [0.0f64; 3];
        for i in 0..3 {
            let a = n00[i] * (1.0 - tx) + n10[i] * tx;
            let b = n01[i] * (1.0 - tx) + n11[i] * tx;
            n[i] = a * (1.0 - tz) + b * tz;
        }
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len > 1e-10 {
            [n[0] / len, n[1] / len, n[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        }
    }
    /// Serialize the height field to a flat `Vec`f64` prefixed by metadata.
    ///
    /// Format: `\[rows, cols, scale_x, scale_z, h0, h1, ...\]`
    pub fn serialize(&self) -> Vec<Real> {
        let mut data = Vec::with_capacity(4 + self.heights.len());
        data.push(self.rows as Real);
        data.push(self.cols as Real);
        data.push(self.scale_x);
        data.push(self.scale_z);
        data.extend_from_slice(&self.heights);
        data
    }
    /// Deserialize a height field from a flat `Vec`f64` (inverse of `serialize`).
    pub fn deserialize(data: &[Real]) -> Option<HeightField> {
        if data.len() < 4 {
            return None;
        }
        let rows = data[0] as usize;
        let cols = data[1] as usize;
        let scale_x = data[2];
        let scale_z = data[3];
        if data.len() != 4 + rows * cols {
            return None;
        }
        Some(HeightField::new(
            data[4..].to_vec(),
            rows,
            cols,
            scale_x,
            scale_z,
        ))
    }
    /// Compute per-vertex normals for all grid vertices.
    ///
    /// Returns a flat `Vec<[Real; 3]>` in row-major order.
    pub fn compute_all_normals(&self) -> Vec<[Real; 3]> {
        let mut normals = Vec::with_capacity(self.rows * self.cols);
        for row in 0..self.rows {
            for col in 0..self.cols {
                normals.push(self.normal_at_grid(col, row));
            }
        }
        normals
    }
    /// Height at arbitrary world-space XZ position (bilinear interpolation).
    ///
    /// Clamps to grid extents.
    pub fn height_at_world(&self, x: Real, z: Real) -> Real {
        let u = (x / (self.scale_x * (self.cols - 1) as Real)).clamp(0.0, 1.0);
        let v = (z / (self.scale_z * (self.rows - 1) as Real)).clamp(0.0, 1.0);
        self.height_at_uv(u, v)
    }
    /// Ray cast with DDA and a maximum number of steps (bounded version).
    ///
    /// Avoids unbounded loops on very large grids.
    pub fn ray_cast_bounded(
        &self,
        origin: [Real; 3],
        dir: [Real; 3],
        max_toi: Real,
        max_steps: usize,
    ) -> Option<(Real, [Real; 3])> {
        let step = self.scale_x.min(self.scale_z) * 0.25;
        let mut t = 0.0;
        let mut best: Option<(Real, [Real; 3])> = None;
        for _ in 0..max_steps {
            if t > max_toi {
                break;
            }
            let px = origin[0] + dir[0] * t;
            let pz = origin[2] + dir[2] * t;
            if px < 0.0
                || px > (self.cols - 1) as Real * self.scale_x
                || pz < 0.0
                || pz > (self.rows - 1) as Real * self.scale_z
            {
                t += step;
                continue;
            }
            let py = origin[1] + dir[1] * t;
            let terrain_y = self.height_at_world(px, pz);
            if py <= terrain_y {
                let refined_t = if t > step { t - step * 0.5 } else { 0.0 };
                let col = ((px / self.scale_x) as usize).min(self.cols.saturating_sub(2));
                let row = ((pz / self.scale_z) as usize).min(self.rows.saturating_sub(2));
                let normal = self.normal_at_grid(col, row);
                best = Some((refined_t.max(0.0), normal));
                break;
            }
            t += step;
        }
        best
    }
}
impl HeightField {
    /// Axis-aligned bounding box of this height field in world space.
    ///
    /// Returns `(min, max)` where each is `[x, y, z]`.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let (min_h, max_h) = self.height_bounds();
        let max_x = (self.cols - 1) as f64 * self.scale_x;
        let max_z = (self.rows - 1) as f64 * self.scale_z;
        ([0.0, min_h, 0.0], [max_x, max_h, max_z])
    }
    /// Bilinear height interpolation at world-space position `(x, z)`.
    ///
    /// Unlike `height_at(row, col)` (grid indices) and `height_at_world`,
    /// this method accepts `(x, z)` world coordinates and returns the
    /// interpolated height, clamping to the grid extents.
    pub fn height_at_xz(&self, x: f64, z: f64) -> f64 {
        self.height_at_world(x, z)
    }
    /// Per-vertex normals for the entire grid in row-major order.
    ///
    /// Each normal is computed via finite differences over the 4-connected
    /// neighborhood and returned as a unit vector `[nx, ny, nz]`.
    pub fn normals(&self) -> Vec<[f64; 3]> {
        self.compute_all_normals()
    }
    /// Tessellate the height field into a triangle mesh.
    ///
    /// Returns `(vertices, indices)` where:
    /// - `vertices` is a list of world-space `[x, y, z]` positions,
    ///   one per grid point, in row-major order.
    /// - `indices` is a list of `[i0, i1, i2]` triangles (two per quad cell).
    pub fn tessellate(&self) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let vertices: Vec<[f64; 3]> = (0..self.rows)
            .flat_map(|row| {
                (0..self.cols).map(move |col| {
                    [
                        col as f64 * self.scale_x,
                        self.height_at(row, col),
                        row as f64 * self.scale_z,
                    ]
                })
            })
            .collect();
        let indices = heightfield_tessellate(self);
        (vertices, indices)
    }
    /// DDA ray intersection returning `(t, normal)` for the first hit.
    ///
    /// Uses grid-cell DDA traversal in the XZ plane, testing the two
    /// triangles of each cell with Möller–Trumbore intersection.
    /// Returns `None` if the ray misses the terrain or is going away from it.
    ///
    /// The search is limited to a reasonable distance (diagonal of the AABB).
    pub fn ray_intersect(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<(f64, [f64; 3])> {
        let (_, max_h) = self.height_bounds();
        if dir[1] > 0.0 && origin[1] >= max_h {
            return None;
        }
        let grid_w = (self.cols - 1) as f64 * self.scale_x;
        let grid_h = (self.rows - 1) as f64 * self.scale_z;
        let diag = (grid_w * grid_w + grid_h * grid_h + (max_h - 0.0).abs().powi(2)).sqrt();
        let max_toi = (origin[1] - max_h).abs() / dir[1].abs().max(1e-12) + diag * 2.0;
        let max_toi = max_toi.min(1e9);
        self.ray_cast_grid(origin, dir, max_toi)
    }
}
impl HeightField {
    /// Compute the slope magnitude at grid point `(col, row)`.
    ///
    /// Slope is defined as `sqrt((dh/dx)² + (dh/dz)²)` using the same
    /// finite-difference stencil as `normal_at_grid`.
    pub fn slope_at(&self, col: usize, row: usize) -> f64 {
        let dh_dx = if col == 0 {
            (self.height_at(row, col + 1) - self.height_at(row, col)) / self.scale_x
        } else if col == self.cols - 1 {
            (self.height_at(row, col) - self.height_at(row, col - 1)) / self.scale_x
        } else {
            (self.height_at(row, col + 1) - self.height_at(row, col - 1)) / (2.0 * self.scale_x)
        };
        let dh_dz = if row == 0 {
            (self.height_at(row + 1, col) - self.height_at(row, col)) / self.scale_z
        } else if row == self.rows - 1 {
            (self.height_at(row, col) - self.height_at(row - 1, col)) / self.scale_z
        } else {
            (self.height_at(row + 1, col) - self.height_at(row - 1, col)) / (2.0 * self.scale_z)
        };
        (dh_dx * dh_dx + dh_dz * dh_dz).sqrt()
    }
    /// Compute the mean curvature at grid point `(col, row)`.
    ///
    /// Uses a second-order central difference approximation:
    /// `κ ≈ (d²h/dx² + d²h/dz²) / 2`.
    pub fn curvature_at(&self, col: usize, row: usize) -> f64 {
        let h = self.height_at(row, col);
        let d2h_dx2 = if col == 0 || col == self.cols - 1 {
            0.0
        } else {
            let hl = self.height_at(row, col - 1);
            let hr = self.height_at(row, col + 1);
            (hl - 2.0 * h + hr) / (self.scale_x * self.scale_x)
        };
        let d2h_dz2 = if row == 0 || row == self.rows - 1 {
            0.0
        } else {
            let hd = self.height_at(row - 1, col);
            let hu = self.height_at(row + 1, col);
            (hd - 2.0 * h + hu) / (self.scale_z * self.scale_z)
        };
        (d2h_dx2 + d2h_dz2) * 0.5
    }
    /// Return a flat `Vec`f64` of slope magnitudes in row-major order.
    pub fn slope_map(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.rows * self.cols);
        for row in 0..self.rows {
            for col in 0..self.cols {
                out.push(self.slope_at(col, row));
            }
        }
        out
    }
    /// Return a flat `Vec`f64` of mean curvatures in row-major order.
    pub fn curvature_map(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.rows * self.cols);
        for row in 0..self.rows {
            for col in 0..self.cols {
                out.push(self.curvature_at(col, row));
            }
        }
        out
    }
    /// Find the closest surface point to query point `q` by brute-force search
    /// over all grid vertices.
    ///
    /// Returns `(closest_world_point, row, col)`.
    pub fn closest_vertex(&self, q: [f64; 3]) -> ([f64; 3], usize, usize) {
        let mut best_dist2 = f64::INFINITY;
        let mut best_pt = [0.0f64; 3];
        let mut best_row = 0usize;
        let mut best_col = 0usize;
        for row in 0..self.rows {
            for col in 0..self.cols {
                let vx = col as f64 * self.scale_x;
                let vy = self.height_at(row, col);
                let vz = row as f64 * self.scale_z;
                let dx = q[0] - vx;
                let dy = q[1] - vy;
                let dz = q[2] - vz;
                let d2 = dx * dx + dy * dy + dz * dz;
                if d2 < best_dist2 {
                    best_dist2 = d2;
                    best_pt = [vx, vy, vz];
                    best_row = row;
                    best_col = col;
                }
            }
        }
        (best_pt, best_row, best_col)
    }
    /// Resample the height field to new grid dimensions by bilinear interpolation.
    ///
    /// `new_cols` and `new_rows` must be ≥ 2.  Returns a new `HeightField` with
    /// the same world extents but a different resolution.
    pub fn resample(&self, new_cols: usize, new_rows: usize) -> HeightField {
        assert!(new_cols >= 2 && new_rows >= 2);
        let new_scale_x = (self.scale_x * (self.cols - 1) as f64) / (new_cols - 1) as f64;
        let new_scale_z = (self.scale_z * (self.rows - 1) as f64) / (new_rows - 1) as f64;
        let mut heights = Vec::with_capacity(new_rows * new_cols);
        for nr in 0..new_rows {
            for nc in 0..new_cols {
                let x = nc as f64 * new_scale_x;
                let z = nr as f64 * new_scale_z;
                heights.push(self.height_at_world(x, z));
            }
        }
        HeightField::new(heights, new_rows, new_cols, new_scale_x, new_scale_z)
    }
    /// Simple hydraulic erosion step: for each cell, water flows to the
    /// steepest-descent neighbor and carries a fraction `sediment_rate` of
    /// height difference.
    ///
    /// `iterations` controls how many passes are performed.  For correctness
    /// the caller should keep `sediment_rate` ≤ 0.5.
    pub fn hydraulic_erode(&mut self, sediment_rate: f64, iterations: usize) {
        for _ in 0..iterations {
            let mut delta = vec![0.0f64; self.rows * self.cols];
            for row in 0..self.rows {
                for col in 0..self.cols {
                    let h = self.height_at(row, col);
                    let neighbors: [(isize, isize); 4] = [
                        (row as isize - 1, col as isize),
                        (row as isize + 1, col as isize),
                        (row as isize, col as isize - 1),
                        (row as isize, col as isize + 1),
                    ];
                    let mut steepest_diff = 0.0f64;
                    let mut steepest_nr = -1isize;
                    let mut steepest_nc = -1isize;
                    for (nr, nc) in neighbors {
                        if nr >= 0 && nr < self.rows as isize && nc >= 0 && nc < self.cols as isize
                        {
                            let nh = self.height_at(nr as usize, nc as usize);
                            let diff = h - nh;
                            if diff > steepest_diff {
                                steepest_diff = diff;
                                steepest_nr = nr;
                                steepest_nc = nc;
                            }
                        }
                    }
                    if steepest_nr >= 0 && steepest_diff > 0.0 {
                        let transfer = sediment_rate * steepest_diff;
                        delta[row * self.cols + col] -= transfer;
                        delta[steepest_nr as usize * self.cols + steepest_nc as usize] += transfer;
                    }
                }
            }
            for (h, d) in self.heights.iter_mut().zip(delta.iter()) {
                *h += d;
            }
        }
    }
    /// Compute a simple flow-accumulation map (watershed proxy).
    ///
    /// Each cell accumulates 1 unit plus the total of all upstream cells that
    /// drain into it following steepest descent.  Returns a flat `Vec`f64`
    /// in row-major order; larger values indicate valleys/channels.
    pub fn flow_accumulation(&self) -> Vec<f64> {
        let n = self.rows * self.cols;
        let mut drain: Vec<Option<usize>> = vec![None; n];
        for row in 0..self.rows {
            for col in 0..self.cols {
                let h = self.height_at(row, col);
                let neighbors: [(isize, isize); 4] = [
                    (row as isize - 1, col as isize),
                    (row as isize + 1, col as isize),
                    (row as isize, col as isize - 1),
                    (row as isize, col as isize + 1),
                ];
                let mut best_diff = 0.0f64;
                let mut best_idx: Option<usize> = None;
                for (nr, nc) in neighbors {
                    if nr >= 0 && nr < self.rows as isize && nc >= 0 && nc < self.cols as isize {
                        let nh = self.height_at(nr as usize, nc as usize);
                        let diff = h - nh;
                        if diff > best_diff {
                            best_diff = diff;
                            best_idx = Some(nr as usize * self.cols + nc as usize);
                        }
                    }
                }
                drain[row * self.cols + col] = best_idx;
            }
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let ha = self.heights[a];
            let hb = self.heights[b];
            hb.partial_cmp(&ha).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut accum = vec![1.0f64; n];
        for &idx in &order {
            if let Some(target) = drain[idx] {
                let val = accum[idx];
                accum[target] += val;
            }
        }
        accum
    }
    /// Clamp all height values to the range `\[min_h, max_h\]`.
    pub fn clamp_heights(&mut self, min_h: f64, max_h: f64) {
        for h in &mut self.heights {
            *h = h.clamp(min_h, max_h);
        }
    }
    /// Scale all heights by a uniform factor.
    pub fn scale_heights(&mut self, factor: f64) {
        for h in &mut self.heights {
            *h *= factor;
        }
    }
    /// Offset all heights by a constant value.
    pub fn offset_heights(&mut self, offset: f64) {
        for h in &mut self.heights {
            *h += offset;
        }
    }
    /// Normalize heights to the range `\[0, 1\]`.  If all heights are equal, all
    /// become 0.
    pub fn normalize_heights(&mut self) {
        let (min_h, max_h) = self.height_bounds();
        let range = max_h - min_h;
        if range < 1e-15 {
            for h in &mut self.heights {
                *h = 0.0;
            }
        } else {
            for h in &mut self.heights {
                *h = (*h - min_h) / range;
            }
        }
    }
    /// Invert heights: each height `h` becomes `max_height - (h - min_height)`.
    pub fn invert_heights(&mut self) {
        let (min_h, max_h) = self.height_bounds();
        for h in &mut self.heights {
            *h = max_h - (*h - min_h);
        }
    }
    /// Compute the average height over the entire grid.
    pub fn mean_height(&self) -> f64 {
        if self.heights.is_empty() {
            return 0.0;
        }
        self.heights.iter().sum::<f64>() / self.heights.len() as f64
    }
    /// Compute the variance of heights.
    pub fn height_variance(&self) -> f64 {
        if self.heights.is_empty() {
            return 0.0;
        }
        let mean = self.mean_height();
        self.heights
            .iter()
            .map(|&h| (h - mean) * (h - mean))
            .sum::<f64>()
            / self.heights.len() as f64
    }
    /// Count the number of local maxima (peaks) in the grid.
    ///
    /// A cell is a peak if it is strictly higher than all 4-connected neighbors.
    pub fn count_peaks(&self) -> usize {
        let mut count = 0usize;
        for row in 0..self.rows {
            for col in 0..self.cols {
                let h = self.height_at(row, col);
                let neighbors: [(isize, isize); 4] = [
                    (row as isize - 1, col as isize),
                    (row as isize + 1, col as isize),
                    (row as isize, col as isize - 1),
                    (row as isize, col as isize + 1),
                ];
                let is_peak = neighbors.iter().all(|&(nr, nc)| {
                    if nr < 0 || nr >= self.rows as isize || nc < 0 || nc >= self.cols as isize {
                        true
                    } else {
                        self.height_at(nr as usize, nc as usize) < h
                    }
                });
                if is_peak {
                    count += 1;
                }
            }
        }
        count
    }
    /// Approximate volume under the height field (above y=0) using the trapezoidal rule.
    ///
    /// Each grid cell contributes `(average_height) * cell_area`.
    pub fn volume(&self) -> f64 {
        if self.rows < 2 || self.cols < 2 {
            return 0.0;
        }
        let mut vol = 0.0_f64;
        let cell_area = self.scale_x * self.scale_z;
        for row in 0..self.rows - 1 {
            for col in 0..self.cols - 1 {
                let h00 = self.height_at(row, col);
                let h10 = self.height_at(row, col + 1);
                let h01 = self.height_at(row + 1, col);
                let h11 = self.height_at(row + 1, col + 1);
                let avg = (h00 + h10 + h01 + h11) * 0.25;
                vol += avg * cell_area;
            }
        }
        vol
    }
    /// Cast a ray against the height field. Returns hit information if intersection found.
    ///
    /// Uses the tessellated triangle mesh for accurate intersection.
    pub fn ray_cast(
        &self,
        ray_origin: &oxiphysics_core::math::Vec3,
        ray_direction: &oxiphysics_core::math::Vec3,
        max_toi: f64,
    ) -> Option<HeightfieldRayHit> {
        let (verts, tris) = self.to_triangle_mesh();
        let mut closest_toi = max_toi;
        let mut hit = None;
        for tri_idx in &tris {
            let v0 = oxiphysics_core::math::Vec3::new(
                verts[tri_idx[0]][0],
                verts[tri_idx[0]][1],
                verts[tri_idx[0]][2],
            );
            let v1 = oxiphysics_core::math::Vec3::new(
                verts[tri_idx[1]][0],
                verts[tri_idx[1]][1],
                verts[tri_idx[1]][2],
            );
            let v2 = oxiphysics_core::math::Vec3::new(
                verts[tri_idx[2]][0],
                verts[tri_idx[2]][1],
                verts[tri_idx[2]][2],
            );
            if let Some(ray_hit) =
                ray_triangle(ray_origin, ray_direction, closest_toi, &v0, &v1, &v2)
                && ray_hit.toi >= 0.0
                && ray_hit.toi < closest_toi
            {
                closest_toi = ray_hit.toi;
                let point = oxiphysics_core::math::Vec3::new(
                    ray_origin.x + ray_direction.x * ray_hit.toi,
                    ray_origin.y + ray_direction.y * ray_hit.toi,
                    ray_origin.z + ray_direction.z * ray_hit.toi,
                );
                hit = Some(HeightfieldRayHit {
                    toi: ray_hit.toi,
                    point,
                });
            }
        }
        hit
    }
}
/// A ray-hit result from `HeightField::ray_cast`.
#[derive(Debug, Clone)]
pub struct HeightfieldRayHit {
    /// Time of impact (distance along ray).
    pub toi: f64,
    /// World-space hit point.
    pub point: oxiphysics_core::math::Vec3,
}
/// Result of a heightfield ray traversal.
#[derive(Debug, Clone)]
pub struct HeightfieldRaycast {
    /// Column index of the cell that was hit.
    pub cell_ix: usize,
    /// Row index of the cell that was hit.
    pub cell_iz: usize,
    /// Ray parameter t at the intersection.
    pub t: f64,
    /// Surface normal at the intersection.
    pub normal: [f64; 3],
}
/// State for a DDA grid-traversal ray cast over a height field.
///
/// Call [`HeightfieldRayTraversal::new`] to set up the traversal, then call
/// [`HeightfieldRayTraversal::next_cell`] repeatedly to visit each cell along
/// the ray in XZ-projection order.
#[derive(Debug, Clone)]
pub struct HeightfieldRayTraversal {
    /// Current cell column index.
    pub cell_col: isize,
    /// Current cell row index.
    pub cell_row: isize,
    /// Direction step in the column dimension (+1 or -1).
    pub(super) step_col: isize,
    /// Direction step in the row dimension (+1 or -1).
    pub(super) step_row: isize,
    /// Parameter t at the next X-boundary crossing.
    pub(super) t_max_x: f64,
    /// Parameter t at the next Z-boundary crossing.
    pub(super) t_max_z: f64,
    /// Change in t per cell step in X.
    pub(super) t_delta_x: f64,
    /// Change in t per cell step in Z.
    pub(super) t_delta_z: f64,
    /// Maximum t to traverse.
    pub(super) t_max: f64,
    /// True once the traversal has finished.
    pub done: bool,
}
impl HeightfieldRayTraversal {
    /// Initialize a DDA traversal over `hf` starting at `ray_origin` in direction
    /// `ray_dir` (need not be normalised).  `max_t` bounds the traversal.
    ///
    /// Returns `None` if the ray does not enter the XZ footprint of the grid.
    pub fn new(
        hf: &HeightField,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_t: f64,
    ) -> Option<Self> {
        if hf.rows < 2 || hf.cols < 2 {
            return None;
        }
        let grid_w = (hf.cols - 1) as f64 * hf.scale_x;
        let grid_h = (hf.rows - 1) as f64 * hf.scale_z;
        let (t_enter, t_exit) = ray_aabb_xz(
            ray_origin[0],
            ray_origin[2],
            ray_dir[0],
            ray_dir[2],
            grid_w,
            grid_h,
            max_t,
        )?;
        let t_start = t_enter.max(0.0);
        if t_start > t_exit {
            return None;
        }
        let px = ray_origin[0] + ray_dir[0] * t_start;
        let pz = ray_origin[2] + ray_dir[2] * t_start;
        let cell_col = ((px / hf.scale_x).floor() as isize).clamp(0, (hf.cols - 2) as isize);
        let cell_row = ((pz / hf.scale_z).floor() as isize).clamp(0, (hf.rows - 2) as isize);
        let step_col = if ray_dir[0] >= 0.0 { 1_isize } else { -1 };
        let step_row = if ray_dir[2] >= 0.0 { 1_isize } else { -1 };
        let t_delta_x = if ray_dir[0].abs() < 1e-12 {
            f64::INFINITY
        } else {
            hf.scale_x / ray_dir[0].abs()
        };
        let t_delta_z = if ray_dir[2].abs() < 1e-12 {
            f64::INFINITY
        } else {
            hf.scale_z / ray_dir[2].abs()
        };
        let x_boundary = if step_col > 0 {
            (cell_col + 1) as f64 * hf.scale_x
        } else {
            cell_col as f64 * hf.scale_x
        };
        let z_boundary = if step_row > 0 {
            (cell_row + 1) as f64 * hf.scale_z
        } else {
            cell_row as f64 * hf.scale_z
        };
        let t_max_x = if ray_dir[0].abs() < 1e-12 {
            f64::INFINITY
        } else {
            t_start + (x_boundary - px) / ray_dir[0]
        };
        let t_max_z = if ray_dir[2].abs() < 1e-12 {
            f64::INFINITY
        } else {
            t_start + (z_boundary - pz) / ray_dir[2]
        };
        Some(Self {
            cell_col,
            cell_row,
            step_col,
            step_row,
            t_max_x,
            t_max_z,
            t_delta_x,
            t_delta_z,
            t_max: t_exit.min(max_t),
            done: false,
        })
    }
    /// Advance to the next cell along the ray.
    ///
    /// Returns `(col, row, t_entry)` for the cell just entered, or `None`
    /// when the traversal exits the grid.
    pub fn next_cell(&mut self) -> Option<(usize, usize, f64)> {
        if self.done {
            return None;
        }
        let col = self.cell_col;
        let row = self.cell_row;
        let t_entry = self.t_max_x.min(self.t_max_z);
        if t_entry > self.t_max {
            self.done = true;
            return None;
        }
        if self.t_max_x < self.t_max_z {
            self.cell_col += self.step_col;
            self.t_max_x += self.t_delta_x;
        } else {
            self.cell_row += self.step_row;
            self.t_max_z += self.t_delta_z;
        }
        if col < 0 || row < 0 {
            self.done = true;
            return None;
        }
        Some((col as usize, row as usize, t_entry))
    }
}
