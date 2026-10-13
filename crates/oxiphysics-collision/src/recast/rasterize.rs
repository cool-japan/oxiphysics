// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Triangle soup rasterization into an open heightfield.

use super::RecastConfig;

/// A single height span in one xz column.
#[derive(Debug, Clone)]
pub struct Span {
    /// Min y in voxel units.
    pub smin: u16,
    /// Max y in voxel units.
    pub smax: u16,
    /// True if this span is walkable (slope within max_slope).
    pub walkable: bool,
    /// Area tag (0 = generic walkable / uncategorized).
    pub area: u8,
}

/// Open heightfield: one column per (ix, iz) cell.
pub struct Heightfield {
    pub width: usize,  // number of columns in X
    pub height: usize, // number of columns in Z
    pub bmin: [f64; 3],
    pub bmax: [f64; 3],
    pub cs: f64, // cell size (XZ)
    pub ch: f64, // cell height (Y)
    /// Column-major storage: spans[iz * width + ix]
    pub spans: Vec<Vec<Span>>,
}

/// Upper bound on the number of heightfield columns (`width * height`).
///
/// Every column owns a `Vec<Span>` header, so a grid this large already needs
/// well over a gigabyte; larger requests (typically caused by a tiny
/// `cell_size` or a huge input extent) are rejected instead of overflowing.
pub const MAX_HEIGHTFIELD_CELLS: usize = 1 << 26;

impl Heightfield {
    /// Create an empty heightfield of `width x height` columns.
    ///
    /// Returns an error when `width * height` overflows or exceeds
    /// [`MAX_HEIGHTFIELD_CELLS`], or when the cell dimensions are not finite
    /// and positive.
    pub fn new(
        width: usize,
        height: usize,
        bmin: [f64; 3],
        bmax: [f64; 3],
        cs: f64,
        ch: f64,
    ) -> Result<Self, String> {
        if !(cs.is_finite() && cs > 0.0 && ch.is_finite() && ch > 0.0) {
            return Err(format!(
                "cell size ({cs}) and cell height ({ch}) must be finite and > 0"
            ));
        }
        let n = width
            .checked_mul(height)
            .filter(|&n| n <= MAX_HEIGHTFIELD_CELLS)
            .ok_or_else(|| {
                format!(
                    "heightfield of {width} x {height} columns exceeds the limit of \
                     {MAX_HEIGHTFIELD_CELLS} columns"
                )
            })?;
        Ok(Self {
            width,
            height,
            bmin,
            bmax,
            cs,
            ch,
            spans: vec![Vec::new(); n],
        })
    }

    #[inline]
    pub fn col_idx(&self, ix: usize, iz: usize) -> usize {
        iz * self.width + ix
    }

    pub fn add_span(&mut self, ix: usize, iz: usize, smin: u16, smax: u16, walkable: bool) {
        let idx = self.col_idx(ix, iz);

        // Check whether this span overlaps or is adjacent to an existing span.
        // If so, merge them rather than creating a duplicate.
        for existing in self.spans[idx].iter_mut() {
            // Overlapping or touching: merge into one span, preserving the wider range.
            if smin <= existing.smax.saturating_add(1) && smax >= existing.smin.saturating_sub(1) {
                existing.smin = existing.smin.min(smin);
                existing.smax = existing.smax.max(smax);
                // Walkable if either source is walkable
                if walkable {
                    existing.walkable = true;
                }
                return;
            }
        }

        // No overlap — insert new span and keep sorted
        self.spans[idx].push(Span {
            smin,
            smax,
            walkable,
            area: 0,
        });
        self.spans[idx].sort_by_key(|s| s.smin);
    }
}

/// Rasterize a triangle soup into a heightfield.
pub fn rasterize_triangles(tris: &[[f64; 9]], cfg: &RecastConfig) -> Result<Heightfield, String> {
    if tris.iter().flatten().any(|c| !c.is_finite()) {
        return Err("input triangles contain non-finite coordinates".to_string());
    }
    // Compute bounding box over all vertices
    let mut bmin = [f64::INFINITY; 3];
    let mut bmax = [f64::NEG_INFINITY; 3];
    for tri in tris {
        for v in 0..3 {
            for ax in 0..3 {
                let c = tri[v * 3 + ax];
                if c < bmin[ax] {
                    bmin[ax] = c;
                }
                if c > bmax[ax] {
                    bmax[ax] = c;
                }
            }
        }
    }
    if bmin.iter().chain(bmax.iter()).any(|c| !c.is_finite()) {
        return Err("input triangle soup has no finite bounding box".to_string());
    }
    if !(cfg.cell_size.is_finite() && cfg.cell_size > 0.0) {
        return Err(format!(
            "cell_size must be finite and > 0, got {}",
            cfg.cell_size
        ));
    }
    if !(cfg.cell_height.is_finite() && cfg.cell_height > 0.0) {
        return Err(format!(
            "cell_height must be finite and > 0, got {}",
            cfg.cell_height
        ));
    }
    // Small padding
    for ax in 0..3 {
        bmin[ax] -= cfg.cell_size;
        bmax[ax] += cfg.cell_size;
    }

    let width = grid_dimension((bmax[0] - bmin[0]) / cfg.cell_size, "X")?;
    let height = grid_dimension((bmax[2] - bmin[2]) / cfg.cell_size, "Z")?;

    // Span heights are stored as u16 voxel indices; reject vertical extents
    // that would silently wrap.
    let vertical_voxels = ((bmax[1] - bmin[1]) / cfg.cell_height).ceil();
    if !vertical_voxels.is_finite() || vertical_voxels + 1.0 >= f64::from(u16::MAX) {
        return Err(format!(
            "vertical extent {} at cell_height {} needs {vertical_voxels} voxels, \
             more than the u16 span range",
            bmax[1] - bmin[1],
            cfg.cell_height
        ));
    }

    let mut hf = Heightfield::new(width, height, bmin, bmax, cfg.cell_size, cfg.cell_height)?;
    let slope_limit = cfg.agent_max_slope.to_radians().cos();

    for tri in tris {
        rasterize_one_triangle(tri, &mut hf, slope_limit);
    }
    Ok(hf)
}

/// Number of grid columns needed to cover `extent_in_cells` cells along one axis.
fn grid_dimension(extent_in_cells: f64, axis: &str) -> Result<usize, String> {
    let cells = extent_in_cells.ceil();
    if !cells.is_finite() || cells > MAX_HEIGHTFIELD_CELLS as f64 {
        return Err(format!(
            "heightfield {axis} dimension of {cells} cells exceeds the limit of \
             {MAX_HEIGHTFIELD_CELLS} columns"
        ));
    }
    Ok((cells as usize).max(1))
}

fn rasterize_one_triangle(tri: &[f64; 9], hf: &mut Heightfield, slope_cos_limit: f64) {
    let v = [
        [tri[0], tri[1], tri[2]],
        [tri[3], tri[4], tri[5]],
        [tri[6], tri[7], tri[8]],
    ];

    // Determine if walkable: check normal Y component against slope limit
    let e0 = [v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]];
    let e1 = [v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]];
    let nx = e0[1] * e1[2] - e0[2] * e1[1];
    let ny = e0[2] * e1[0] - e0[0] * e1[2];
    let nz = e0[0] * e1[1] - e0[1] * e1[0];
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    let walkable = if len > 1e-10 {
        (ny / len).abs() >= slope_cos_limit
    } else {
        false
    };

    // AABB of triangle in XZ
    let xmin = v[0][0].min(v[1][0]).min(v[2][0]);
    let xmax = v[0][0].max(v[1][0]).max(v[2][0]);
    let zmin = v[0][2].min(v[1][2]).min(v[2][2]);
    let zmax = v[0][2].max(v[1][2]).max(v[2][2]);

    let ix0 = (((xmin - hf.bmin[0]) / hf.cs).floor() as isize).max(0) as usize;
    let ix1 = (((xmax - hf.bmin[0]) / hf.cs).ceil() as isize).min(hf.width as isize - 1) as usize;
    let iz0 = (((zmin - hf.bmin[2]) / hf.cs).floor() as isize).max(0) as usize;
    let iz1 = (((zmax - hf.bmin[2]) / hf.cs).ceil() as isize).min(hf.height as isize - 1) as usize;

    for iz in iz0..=iz1 {
        for ix in ix0..=ix1 {
            let col_xmin = hf.bmin[0] + ix as f64 * hf.cs;
            let col_zmin = hf.bmin[2] + iz as f64 * hf.cs;
            let col_xmax = col_xmin + hf.cs;
            let col_zmax = col_zmin + hf.cs;

            if let Some((ymin, ymax)) =
                triangle_y_range_in_column(&v, col_xmin, col_xmax, col_zmin, col_zmax)
            {
                let smin = (((ymin - hf.bmin[1]) / hf.ch).floor() as i32).max(0) as u16;
                let smax =
                    (((ymax - hf.bmin[1]) / hf.ch).ceil() as i32).max(smin as i32 + 1) as u16;
                hf.add_span(ix, iz, smin, smax, walkable);
            }
        }
    }
}

/// Compute the Y range of the triangle projected onto the XZ column [xmin,xmax]×[zmin,zmax].
///
/// Uses Sutherland-Hodgman polygon clipping to correctly handle the case where the
/// column is completely inside the triangle (no edge intersection, but overlap exists).
///
/// Returns `None` if the triangle doesn't overlap the column at all.
fn triangle_y_range_in_column(
    v: &[[f64; 3]; 3],
    col_xmin: f64,
    col_xmax: f64,
    col_zmin: f64,
    col_zmax: f64,
) -> Option<(f64, f64)> {
    // Quick AABB reject
    let tri_xmin = v[0][0].min(v[1][0]).min(v[2][0]);
    let tri_xmax = v[0][0].max(v[1][0]).max(v[2][0]);
    let tri_zmin = v[0][2].min(v[1][2]).min(v[2][2]);
    let tri_zmax = v[0][2].max(v[1][2]).max(v[2][2]);
    if tri_xmax < col_xmin || tri_xmin > col_xmax {
        return None;
    }
    if tri_zmax < col_zmin || tri_zmin > col_zmax {
        return None;
    }

    // Sutherland-Hodgman clipping in XZ.
    // Each polygon vertex is (x, y, z).
    let mut poly: Vec<[f64; 3]> = v.to_vec();

    // Clip against each of the 4 half-planes in XZ.
    // We clip one half-plane at a time and update `poly`.
    for (axis, lo, inside_sign) in [
        (0usize, col_xmin, 1.0f64), // x >= col_xmin
        (0, col_xmax, -1.0),        // x <= col_xmax  →  -x >= -col_xmax
        (2usize, col_zmin, 1.0),    // z >= col_zmin
        (2, col_zmax, -1.0),        // z <= col_zmax
    ] {
        if poly.is_empty() {
            return None;
        }
        let boundary = if inside_sign > 0.0 { lo } else { -lo };
        let n = poly.len();
        let mut clipped: Vec<[f64; 3]> = Vec::with_capacity(n + 1);
        for i in 0..n {
            let a = poly[i];
            let b = poly[(i + 1) % n];
            let da = inside_sign * a[axis] - boundary;
            let db = inside_sign * b[axis] - boundary;
            if da >= 0.0 {
                clipped.push(a);
            }
            if (da < 0.0) != (db < 0.0) {
                // Edge crosses boundary
                let t = da / (da - db);
                clipped.push([
                    a[0] + t * (b[0] - a[0]),
                    a[1] + t * (b[1] - a[1]),
                    a[2] + t * (b[2] - a[2]),
                ]);
            }
        }
        poly = clipped;
    }

    if poly.is_empty() {
        return None;
    }

    let ymin = poly.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
    let ymax = poly.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
    Some((ymin, ymax))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heightfield_new_rejects_overflowing_dimensions() {
        let res = Heightfield::new(usize::MAX, 2, [0.0; 3], [1.0; 3], 0.3, 0.2);
        assert!(res.is_err(), "width * height overflow must be an error");
        let res = Heightfield::new(1 << 14, 1 << 13, [0.0; 3], [1.0; 3], 0.3, 0.2);
        assert!(res.is_err(), "2^27 columns exceed MAX_HEIGHTFIELD_CELLS");
        let hf = Heightfield::new(4, 3, [0.0; 3], [1.0; 3], 0.3, 0.2).expect("small grid");
        assert_eq!(hf.spans.len(), 12);
    }

    #[test]
    fn test_heightfield_new_rejects_bad_cell_dimensions() {
        assert!(Heightfield::new(2, 2, [0.0; 3], [1.0; 3], 0.0, 0.2).is_err());
        assert!(Heightfield::new(2, 2, [0.0; 3], [1.0; 3], 0.3, f64::NAN).is_err());
    }

    #[test]
    fn test_rasterize_rejects_vertical_extent_beyond_u16() {
        let cfg = RecastConfig::default();
        // 1e5 m tall sliver at 0.2 m voxels needs 5e5 voxels > u16::MAX.
        let tris = [[0.0, 0.0, 0.0, 1.0, 1.0e5, 0.0, 1.0, 0.0, 1.0]];
        assert!(rasterize_triangles(&tris, &cfg).is_err());
    }

    #[test]
    fn test_rasterize_rejects_zero_cell_size() {
        let cfg = RecastConfig {
            cell_size: 0.0,
            ..Default::default()
        };
        let tris = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]];
        assert!(rasterize_triangles(&tris, &cfg).is_err());
    }
}
