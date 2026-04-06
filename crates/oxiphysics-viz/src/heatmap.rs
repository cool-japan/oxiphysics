#![allow(clippy::too_many_arguments)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Heatmap and 2D field visualization.
//!
//! Provides data structures and algorithms for rendering 2D scalar/vector fields:
//!
//! - [`HeatmapGrid`] — 2D grid of `f32` values with bilinear sampling
//! - [`HeatmapRenderer`] — render grid to RGBA pixel buffer using a colormap
//! - [`ContourLine`] — isoline at a given level via marching squares
//! - [`ContourSet`] — collection of isolines at multiple levels
//! - [`VectorField2D`] — 2D velocity/vector field; magnitude, divergence, curl
//! - [`ArrowGlyph`] — single arrow glyph at a position
//! - [`StreamlineIntegrator`] — RK4 streamline integration in a vector field
//! - [`SeedPoints`] — uniform / random seed point generation
//! - [`GradientField`] — gradient of a scalar field via central differences
//! - [`HeatmapAnimation`] — sequence of `HeatmapGrid` frames with interpolation

#![allow(dead_code)]

use crate::colormap::{Colormap, map_scalar};

// ---------------------------------------------------------------------------
// HeatmapGrid
// ---------------------------------------------------------------------------

/// A 2D grid of `f32` scalar values of size `nx × ny`.
///
/// Values are stored row-major: `data[y * nx + x]`.
/// The grid spans `[x_min, x_max] × [y_min, y_max]`.
#[derive(Debug, Clone)]
pub struct HeatmapGrid {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// Scalar data, row-major (`data[y * nx + x]`).
    pub data: Vec<f32>,
    /// Left boundary of the grid domain.
    pub x_min: f32,
    /// Right boundary of the grid domain.
    pub x_max: f32,
    /// Bottom boundary of the grid domain.
    pub y_min: f32,
    /// Top boundary of the grid domain.
    pub y_max: f32,
}

impl HeatmapGrid {
    /// Create a new heatmap grid, filling all values with `fill`.
    pub fn new(
        nx: usize,
        ny: usize,
        x_min: f32,
        x_max: f32,
        y_min: f32,
        y_max: f32,
        fill: f32,
    ) -> Self {
        Self {
            nx,
            ny,
            data: vec![fill; nx * ny],
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    /// Create a grid from existing data. Panics if `data.len() != nx * ny`.
    pub fn from_data(
        nx: usize,
        ny: usize,
        data: Vec<f32>,
        x_min: f32,
        x_max: f32,
        y_min: f32,
        y_max: f32,
    ) -> Self {
        assert_eq!(data.len(), nx * ny, "data length must equal nx * ny");
        Self {
            nx,
            ny,
            data,
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    /// Get value at grid cell (ix, iy). Returns 0 if out of bounds.
    pub fn get(&self, ix: usize, iy: usize) -> f32 {
        if ix < self.nx && iy < self.ny {
            self.data[iy * self.nx + ix]
        } else {
            0.0
        }
    }

    /// Set value at grid cell (ix, iy). No-op if out of bounds.
    pub fn set(&mut self, ix: usize, iy: usize, v: f32) {
        if ix < self.nx && iy < self.ny {
            self.data[iy * self.nx + ix] = v;
        }
    }

    /// Minimum value in the grid.
    pub fn min_val(&self) -> f32 {
        self.data.iter().copied().fold(f32::INFINITY, f32::min)
    }

    /// Maximum value in the grid.
    pub fn max_val(&self) -> f32 {
        self.data.iter().copied().fold(f32::NEG_INFINITY, f32::max)
    }

    /// Bilinearly sample the grid at world coordinates (x, y).
    /// Clamps to grid boundaries.
    pub fn sample(&self, x: f32, y: f32) -> f32 {
        if self.nx == 0 || self.ny == 0 {
            return 0.0;
        }
        let dx = (self.x_max - self.x_min) / (self.nx as f32 - 1.0).max(1.0);
        let dy = (self.y_max - self.y_min) / (self.ny as f32 - 1.0).max(1.0);
        let fx = ((x - self.x_min) / dx).clamp(0.0, (self.nx - 1) as f32);
        let fy = ((y - self.y_min) / dy).clamp(0.0, (self.ny - 1) as f32);
        let ix0 = fx.floor() as usize;
        let iy0 = fy.floor() as usize;
        let ix1 = (ix0 + 1).min(self.nx - 1);
        let iy1 = (iy0 + 1).min(self.ny - 1);
        let tx = fx - fx.floor();
        let ty = fy - fy.floor();
        let v00 = self.get(ix0, iy0);
        let v10 = self.get(ix1, iy0);
        let v01 = self.get(ix0, iy1);
        let v11 = self.get(ix1, iy1);
        // Bilinear interpolation
        let v0 = v00 * (1.0 - tx) + v10 * tx;
        let v1 = v01 * (1.0 - tx) + v11 * tx;
        v0 * (1.0 - ty) + v1 * ty
    }

    /// Cell width (spacing between adjacent columns).
    pub fn dx(&self) -> f32 {
        if self.nx <= 1 {
            1.0
        } else {
            (self.x_max - self.x_min) / (self.nx - 1) as f32
        }
    }

    /// Cell height (spacing between adjacent rows).
    pub fn dy(&self) -> f32 {
        if self.ny <= 1 {
            1.0
        } else {
            (self.y_max - self.y_min) / (self.ny - 1) as f32
        }
    }
}

// ---------------------------------------------------------------------------
// HeatmapRenderer
// ---------------------------------------------------------------------------

/// Renders a [`HeatmapGrid`] to an RGBA pixel buffer using a colormap.
///
/// Output buffer is `width × height` pixels, each 4 bytes (R, G, B, A).
#[derive(Debug, Clone)]
pub struct HeatmapRenderer {
    /// Output image width in pixels.
    pub width: usize,
    /// Output image height in pixels.
    pub height: usize,
    /// Colormap to use for scalar mapping.
    pub colormap: Colormap,
    /// Override minimum value for color mapping. If `None`, uses grid min.
    pub vmin: Option<f32>,
    /// Override maximum value for color mapping. If `None`, uses grid max.
    pub vmax: Option<f32>,
}

impl HeatmapRenderer {
    /// Create a new renderer.
    pub fn new(width: usize, height: usize, colormap: Colormap) -> Self {
        Self {
            width,
            height,
            colormap,
            vmin: None,
            vmax: None,
        }
    }

    /// Set explicit value range for color mapping.
    pub fn with_range(mut self, vmin: f32, vmax: f32) -> Self {
        self.vmin = Some(vmin);
        self.vmax = Some(vmax);
        self
    }

    /// Render the grid to an RGBA pixel buffer (length = width * height * 4).
    pub fn render(&self, grid: &HeatmapGrid) -> Vec<u8> {
        let vmin = self.vmin.unwrap_or_else(|| grid.min_val()) as f64;
        let vmax = self.vmax.unwrap_or_else(|| grid.max_val()) as f64;
        let mut pixels = vec![0u8; self.width * self.height * 4];
        for py in 0..self.height {
            for px in 0..self.width {
                // Map pixel to world coordinates
                let x = grid.x_min
                    + (px as f32 / (self.width as f32 - 1.0).max(1.0)) * (grid.x_max - grid.x_min);
                let y = grid.y_min
                    + (py as f32 / (self.height as f32 - 1.0).max(1.0)) * (grid.y_max - grid.y_min);
                let val = grid.sample(x, y) as f64;
                let color = map_scalar(val, vmin, vmax, self.colormap);
                let idx = (py * self.width + px) * 4;
                pixels[idx] = (color.r * 255.0).clamp(0.0, 255.0) as u8;
                pixels[idx + 1] = (color.g * 255.0).clamp(0.0, 255.0) as u8;
                pixels[idx + 2] = (color.b * 255.0).clamp(0.0, 255.0) as u8;
                pixels[idx + 3] = (color.a * 255.0).clamp(0.0, 255.0) as u8;
            }
        }
        pixels
    }
}

// ---------------------------------------------------------------------------
// ContourLine — marching squares
// ---------------------------------------------------------------------------

/// A single isoline (contour) at a given scalar level.
///
/// Generated via the marching squares algorithm on a [`HeatmapGrid`].
#[derive(Debug, Clone)]
pub struct ContourLine {
    /// Scalar level for this contour.
    pub level: f32,
    /// Sequence of 2D points `[x, y]` forming the isoline.
    pub points: Vec<[f32; 2]>,
}

impl ContourLine {
    /// Extract an isoline at `level` from `grid` using marching squares.
    ///
    /// Returns line segments as pairs of endpoints. Each segment is two consecutive
    /// entries in `points`: `[p0, p1, p2, p3, ...]` where `(p0,p1)`, `(p2,p3)` are segments.
    pub fn extract(grid: &HeatmapGrid, level: f32) -> Self {
        let mut points = Vec::new();
        if grid.nx < 2 || grid.ny < 2 {
            return Self { level, points };
        }
        let dx = grid.dx();
        let dy = grid.dy();

        for iy in 0..(grid.ny - 1) {
            for ix in 0..(grid.nx - 1) {
                let v00 = grid.get(ix, iy);
                let v10 = grid.get(ix + 1, iy);
                let v01 = grid.get(ix, iy + 1);
                let v11 = grid.get(ix + 1, iy + 1);

                // Marching squares: compute case index
                let b00 = if v00 >= level { 1u8 } else { 0u8 };
                let b10 = if v10 >= level { 2u8 } else { 0u8 };
                let b11 = if v11 >= level { 4u8 } else { 0u8 };
                let b01 = if v01 >= level { 8u8 } else { 0u8 };
                let case = b00 | b10 | b11 | b01;

                // World coordinates of cell corners
                let x0 = grid.x_min + ix as f32 * dx;
                let x1 = x0 + dx;
                let y0 = grid.y_min + iy as f32 * dy;
                let y1 = y0 + dy;

                // Linear interpolation along an edge
                let interp = |va: f32, vb: f32, pa: [f32; 2], pb: [f32; 2]| -> [f32; 2] {
                    let t = if (vb - va).abs() < 1e-30 {
                        0.5
                    } else {
                        (level - va) / (vb - va)
                    };
                    [pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1])]
                };

                let p00 = [x0, y0];
                let p10 = [x1, y0];
                let p01 = [x0, y1];
                let p11 = [x1, y1];

                // Edge midpoints via interpolation
                let e_bottom = || interp(v00, v10, p00, p10); // bottom edge
                let e_right = || interp(v10, v11, p10, p11); // right edge
                let e_top = || interp(v01, v11, p01, p11); // top edge
                let e_left = || interp(v00, v01, p00, p01); // left edge

                // Generate line segments based on marching squares case
                match case {
                    0 | 15 => {}
                    1 => {
                        points.push(e_bottom());
                        points.push(e_left());
                    }
                    2 => {
                        points.push(e_bottom());
                        points.push(e_right());
                    }
                    3 => {
                        points.push(e_left());
                        points.push(e_right());
                    }
                    4 => {
                        points.push(e_right());
                        points.push(e_top());
                    }
                    5 => {
                        points.push(e_bottom());
                        points.push(e_right());
                        points.push(e_top());
                        points.push(e_left());
                    }
                    6 => {
                        points.push(e_bottom());
                        points.push(e_top());
                    }
                    7 => {
                        points.push(e_top());
                        points.push(e_left());
                    }
                    8 => {
                        points.push(e_top());
                        points.push(e_left());
                    }
                    9 => {
                        points.push(e_bottom());
                        points.push(e_top());
                    }
                    10 => {
                        points.push(e_bottom());
                        points.push(e_left());
                        points.push(e_right());
                        points.push(e_top());
                    }
                    11 => {
                        points.push(e_right());
                        points.push(e_top());
                    }
                    12 => {
                        points.push(e_left());
                        points.push(e_right());
                    }
                    13 => {
                        points.push(e_bottom());
                        points.push(e_right());
                    }
                    14 => {
                        points.push(e_bottom());
                        points.push(e_left());
                    }
                    _ => {}
                }
            }
        }
        Self { level, points }
    }

    /// Number of line segments (each segment = 2 consecutive points).
    pub fn segment_count(&self) -> usize {
        self.points.len() / 2
    }

    /// Total arc length of all segments.
    pub fn arc_length(&self) -> f32 {
        let mut total = 0.0f32;
        let mut i = 0;
        while i + 1 < self.points.len() {
            let [x0, y0] = self.points[i];
            let [x1, y1] = self.points[i + 1];
            let dx = x1 - x0;
            let dy = y1 - y0;
            total += (dx * dx + dy * dy).sqrt();
            i += 2;
        }
        total
    }
}

// ---------------------------------------------------------------------------
// ContourSet
// ---------------------------------------------------------------------------

/// A collection of [`ContourLine`]s at multiple levels.
#[derive(Debug, Clone)]
pub struct ContourSet {
    /// All contour lines in this set.
    pub contours: Vec<ContourLine>,
}

impl ContourSet {
    /// Generate contours from a grid at multiple levels.
    pub fn extract(grid: &HeatmapGrid, levels: &[f32]) -> Self {
        let contours = levels
            .iter()
            .map(|&lv| ContourLine::extract(grid, lv))
            .collect();
        Self { contours }
    }

    /// Generate `n` equally-spaced contour levels between `vmin` and `vmax`.
    pub fn extract_n(grid: &HeatmapGrid, n: usize, vmin: f32, vmax: f32) -> Self {
        if n == 0 {
            return Self { contours: vec![] };
        }
        let levels: Vec<f32> = (0..n)
            .map(|i| vmin + (i as f32 / (n - 1).max(1) as f32) * (vmax - vmin))
            .collect();
        Self::extract(grid, &levels)
    }

    /// Find the contour line closest to the given world point `[x, y]`.
    /// Returns the index of the closest contour in `self.contours`.
    pub fn closest_contour_idx(&self, point: [f32; 2]) -> Option<usize> {
        let [px, py] = point;
        let mut best_idx = None;
        let mut best_dist = f32::INFINITY;
        for (idx, contour) in self.contours.iter().enumerate() {
            let mut i = 0;
            while i + 1 < contour.points.len() {
                let [x0, y0] = contour.points[i];
                let [x1, y1] = contour.points[i + 1];
                // Distance from point to segment
                let d = point_to_segment_dist(px, py, x0, y0, x1, y1);
                if d < best_dist {
                    best_dist = d;
                    best_idx = Some(idx);
                }
                i += 2;
            }
        }
        best_idx
    }

    /// Total number of contour line segments across all levels.
    pub fn total_segment_count(&self) -> usize {
        self.contours.iter().map(|c| c.segment_count()).sum()
    }
}

/// Distance from point (px, py) to segment (x0,y0)-(x1,y1).
fn point_to_segment_dist(px: f32, py: f32, x0: f32, y0: f32, x1: f32, y1: f32) -> f32 {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-30 {
        let ex = px - x0;
        let ey = py - y0;
        return (ex * ex + ey * ey).sqrt();
    }
    let t = ((px - x0) * dx + (py - y0) * dy) / len2;
    let t = t.clamp(0.0, 1.0);
    let cx = x0 + t * dx - px;
    let cy = y0 + t * dy - py;
    (cx * cx + cy * cy).sqrt()
}

// ---------------------------------------------------------------------------
// VectorField2D
// ---------------------------------------------------------------------------

/// A 2D vector field over a uniform grid.
///
/// Each cell stores a `(u, v)` velocity vector.
/// Grid layout: `u_data[y * nx + x]`, `v_data[y * nx + x]`.
#[derive(Debug, Clone)]
pub struct VectorField2D {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// U-component (x-direction) of vectors, row-major.
    pub u_data: Vec<f32>,
    /// V-component (y-direction) of vectors, row-major.
    pub v_data: Vec<f32>,
    /// Left boundary.
    pub x_min: f32,
    /// Right boundary.
    pub x_max: f32,
    /// Bottom boundary.
    pub y_min: f32,
    /// Top boundary.
    pub y_max: f32,
}

impl VectorField2D {
    /// Create a vector field from u and v arrays.
    pub fn new(
        nx: usize,
        ny: usize,
        u_data: Vec<f32>,
        v_data: Vec<f32>,
        x_min: f32,
        x_max: f32,
        y_min: f32,
        y_max: f32,
    ) -> Self {
        assert_eq!(u_data.len(), nx * ny);
        assert_eq!(v_data.len(), nx * ny);
        Self {
            nx,
            ny,
            u_data,
            v_data,
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    /// Get `(u, v)` at grid cell `(ix, iy)`.
    pub fn get(&self, ix: usize, iy: usize) -> (f32, f32) {
        if ix < self.nx && iy < self.ny {
            (
                self.u_data[iy * self.nx + ix],
                self.v_data[iy * self.nx + ix],
            )
        } else {
            (0.0, 0.0)
        }
    }

    /// Bilinearly sample `(u, v)` at world position `(x, y)`.
    pub fn sample(&self, x: f32, y: f32) -> (f32, f32) {
        if self.nx == 0 || self.ny == 0 {
            return (0.0, 0.0);
        }
        let dx = if self.nx > 1 {
            (self.x_max - self.x_min) / (self.nx - 1) as f32
        } else {
            1.0
        };
        let dy = if self.ny > 1 {
            (self.y_max - self.y_min) / (self.ny - 1) as f32
        } else {
            1.0
        };
        let fx = ((x - self.x_min) / dx).clamp(0.0, (self.nx - 1) as f32);
        let fy = ((y - self.y_min) / dy).clamp(0.0, (self.ny - 1) as f32);
        let ix0 = fx.floor() as usize;
        let iy0 = fy.floor() as usize;
        let ix1 = (ix0 + 1).min(self.nx - 1);
        let iy1 = (iy0 + 1).min(self.ny - 1);
        let tx = fx - fx.floor();
        let ty = fy - fy.floor();

        let (u00, v00) = self.get(ix0, iy0);
        let (u10, v10) = self.get(ix1, iy0);
        let (u01, v01) = self.get(ix0, iy1);
        let (u11, v11) = self.get(ix1, iy1);

        let u = (u00 * (1.0 - tx) + u10 * tx) * (1.0 - ty) + (u01 * (1.0 - tx) + u11 * tx) * ty;
        let v = (v00 * (1.0 - tx) + v10 * tx) * (1.0 - ty) + (v01 * (1.0 - tx) + v11 * tx) * ty;
        (u, v)
    }

    /// Compute magnitude field (as `HeatmapGrid`).
    pub fn magnitude_grid(&self) -> HeatmapGrid {
        let data: Vec<f32> = self
            .u_data
            .iter()
            .zip(self.v_data.iter())
            .map(|(&u, &v)| (u * u + v * v).sqrt())
            .collect();
        HeatmapGrid::from_data(
            self.nx, self.ny, data, self.x_min, self.x_max, self.y_min, self.y_max,
        )
    }

    /// Compute discrete divergence: ∂u/∂x + ∂v/∂y (central differences).
    pub fn divergence(&self) -> HeatmapGrid {
        let dx = if self.nx > 2 {
            (self.x_max - self.x_min) / (self.nx - 1) as f32
        } else {
            1.0
        };
        let dy = if self.ny > 2 {
            (self.y_max - self.y_min) / (self.ny - 1) as f32
        } else {
            1.0
        };
        let mut div = vec![0.0f32; self.nx * self.ny];
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                let (u_xp, _) = self.get(ix.saturating_add(1).min(self.nx - 1), iy);
                let (u_xm, _) = self.get(ix.saturating_sub(1), iy);
                let (_, v_yp) = self.get(ix, iy.saturating_add(1).min(self.ny - 1));
                let (_, v_ym) = self.get(ix, iy.saturating_sub(1));
                let dscale_x = if ix == 0 || ix == self.nx - 1 {
                    1.0
                } else {
                    2.0
                };
                let dscale_y = if iy == 0 || iy == self.ny - 1 {
                    1.0
                } else {
                    2.0
                };
                div[iy * self.nx + ix] =
                    (u_xp - u_xm) / (dscale_x * dx) + (v_yp - v_ym) / (dscale_y * dy);
            }
        }
        HeatmapGrid::from_data(
            self.nx, self.ny, div, self.x_min, self.x_max, self.y_min, self.y_max,
        )
    }

    /// Compute discrete curl (z-component): ∂v/∂x - ∂u/∂y (central differences).
    pub fn curl(&self) -> HeatmapGrid {
        let dx = if self.nx > 2 {
            (self.x_max - self.x_min) / (self.nx - 1) as f32
        } else {
            1.0
        };
        let dy = if self.ny > 2 {
            (self.y_max - self.y_min) / (self.ny - 1) as f32
        } else {
            1.0
        };
        let mut curl_data = vec![0.0f32; self.nx * self.ny];
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                let (_, v_xp) = self.get(ix.saturating_add(1).min(self.nx - 1), iy);
                let (_, v_xm) = self.get(ix.saturating_sub(1), iy);
                let (u_yp, _) = self.get(ix, iy.saturating_add(1).min(self.ny - 1));
                let (u_ym, _) = self.get(ix, iy.saturating_sub(1));
                let dscale_x = if ix == 0 || ix == self.nx - 1 {
                    1.0
                } else {
                    2.0
                };
                let dscale_y = if iy == 0 || iy == self.ny - 1 {
                    1.0
                } else {
                    2.0
                };
                curl_data[iy * self.nx + ix] =
                    (v_xp - v_xm) / (dscale_x * dx) - (u_yp - u_ym) / (dscale_y * dy);
            }
        }
        HeatmapGrid::from_data(
            self.nx, self.ny, curl_data, self.x_min, self.x_max, self.y_min, self.y_max,
        )
    }
}

// ---------------------------------------------------------------------------
// ArrowGlyph
// ---------------------------------------------------------------------------

/// An arrow glyph at a 2D position with direction and magnitude.
#[derive(Debug, Clone, Copy)]
pub struct ArrowGlyph {
    /// Position of the arrow base `[x, y]`.
    pub position: [f32; 2],
    /// Direction vector `[dx, dy]` (normalized internally when rendering).
    pub direction: [f32; 2],
    /// Magnitude (length of the arrow before clamping).
    pub magnitude: f32,
    /// Maximum arrow length (clamps the drawn length).
    pub max_length: f32,
}

impl ArrowGlyph {
    /// Create an arrow glyph.
    pub fn new(position: [f32; 2], direction: [f32; 2], magnitude: f32, max_length: f32) -> Self {
        Self {
            position,
            direction,
            magnitude,
            max_length,
        }
    }

    /// Compute arrow tip position, clamping length to `max_length`.
    pub fn tip(&self) -> [f32; 2] {
        let [dx, dy] = self.direction;
        let dir_len = (dx * dx + dy * dy).sqrt();
        if dir_len < 1e-30 {
            return self.position;
        }
        let length = self.magnitude.min(self.max_length);
        let [bx, by] = self.position;
        [bx + (dx / dir_len) * length, by + (dy / dir_len) * length]
    }

    /// Actual drawn length (clamped).
    pub fn drawn_length(&self) -> f32 {
        self.magnitude.min(self.max_length)
    }

    /// Generate arrow as a line segment: `[base, tip]`.
    pub fn as_segment(&self) -> [[f32; 2]; 2] {
        [self.position, self.tip()]
    }
}

/// Generate arrow glyphs for a vector field on a coarse subgrid.
///
/// `step_x` and `step_y` define the sampling stride in grid cells.
pub fn arrow_glyphs_from_field(
    field: &VectorField2D,
    step_x: usize,
    step_y: usize,
    max_length: f32,
) -> Vec<ArrowGlyph> {
    let mut glyphs = Vec::new();
    let dx = field.nx.max(2) - 1;
    let dy = field.ny.max(2) - 1;
    let cell_x = (field.x_max - field.x_min) / dx as f32;
    let cell_y = (field.y_max - field.y_min) / dy as f32;

    let step_x = step_x.max(1);
    let step_y = step_y.max(1);
    let mut iy = 0;
    while iy < field.ny {
        let mut ix = 0;
        while ix < field.nx {
            let x = field.x_min + ix as f32 * cell_x;
            let y = field.y_min + iy as f32 * cell_y;
            let (u, v) = field.get(ix, iy);
            let mag = (u * u + v * v).sqrt();
            glyphs.push(ArrowGlyph::new([x, y], [u, v], mag, max_length));
            ix += step_x;
        }
        iy += step_y;
    }
    glyphs
}

// ---------------------------------------------------------------------------
// StreamlineIntegrator
// ---------------------------------------------------------------------------

/// Integrates streamlines in a [`VectorField2D`] using RK4.
#[derive(Debug, Clone)]
pub struct StreamlineIntegrator {
    /// Integration step size.
    pub dt: f32,
    /// Maximum number of steps per streamline.
    pub max_steps: usize,
    /// Stop when velocity magnitude is below this threshold.
    pub min_velocity: f32,
}

/// A single integrated streamline.
#[derive(Debug, Clone)]
pub struct Streamline2D {
    /// Path points `[x, y]`.
    pub points: Vec<[f32; 2]>,
}

impl StreamlineIntegrator {
    /// Create a new integrator.
    pub fn new(dt: f32, max_steps: usize, min_velocity: f32) -> Self {
        Self {
            dt,
            max_steps,
            min_velocity,
        }
    }

    /// Integrate a streamline from seed point `[x, y]` in `field`.
    pub fn integrate(&self, field: &VectorField2D, seed: [f32; 2]) -> Streamline2D {
        let mut points = Vec::with_capacity(self.max_steps + 1);
        let [mut x, mut y] = seed;
        points.push([x, y]);

        for _ in 0..self.max_steps {
            // RK4
            let (u1, v1) = field.sample(x, y);
            if (u1 * u1 + v1 * v1).sqrt() < self.min_velocity {
                break;
            }
            let (u2, v2) = field.sample(x + 0.5 * self.dt * u1, y + 0.5 * self.dt * v1);
            let (u3, v3) = field.sample(x + 0.5 * self.dt * u2, y + 0.5 * self.dt * v2);
            let (u4, v4) = field.sample(x + self.dt * u3, y + self.dt * v3);
            x += self.dt * (u1 + 2.0 * u2 + 2.0 * u3 + u4) / 6.0;
            y += self.dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0;

            // Clamp to field domain
            if x < field.x_min || x > field.x_max || y < field.y_min || y > field.y_max {
                break;
            }
            points.push([x, y]);
        }
        Streamline2D { points }
    }

    /// Integrate multiple streamlines from a list of seed points.
    pub fn integrate_all(&self, field: &VectorField2D, seeds: &[[f32; 2]]) -> Vec<Streamline2D> {
        seeds
            .iter()
            .map(|&seed| self.integrate(field, seed))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// SeedPoints
// ---------------------------------------------------------------------------

/// Seed point generation for streamline integration.
#[derive(Debug, Clone)]
pub struct SeedPoints {
    /// Generated seed points `[x, y]`.
    pub points: Vec<[f32; 2]>,
}

impl SeedPoints {
    /// Generate `nx × ny` uniformly-spaced seed points in the domain
    /// `[x_min, x_max] × [y_min, y_max]`.
    pub fn uniform(nx: usize, ny: usize, x_min: f32, x_max: f32, y_min: f32, y_max: f32) -> Self {
        let mut points = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            for ix in 0..nx {
                let x = if nx <= 1 {
                    (x_min + x_max) * 0.5
                } else {
                    x_min + (ix as f32 / (nx - 1) as f32) * (x_max - x_min)
                };
                let y = if ny <= 1 {
                    (y_min + y_max) * 0.5
                } else {
                    y_min + (iy as f32 / (ny - 1) as f32) * (y_max - y_min)
                };
                points.push([x, y]);
            }
        }
        Self { points }
    }

    /// Generate `n` random seed points in the domain using a simple LCG RNG.
    /// Uses a deterministic seed for reproducibility.
    pub fn random(n: usize, x_min: f32, x_max: f32, y_min: f32, y_max: f32, seed: u64) -> Self {
        let mut lcg = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut points = Vec::with_capacity(n);
        for _ in 0..n {
            lcg = lcg
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let rx = (lcg >> 33) as f32 / (u32::MAX as f32);
            lcg = lcg
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let ry = (lcg >> 33) as f32 / (u32::MAX as f32);
            let x = x_min + rx * (x_max - x_min);
            let y = y_min + ry * (y_max - y_min);
            points.push([x, y]);
        }
        Self { points }
    }

    /// Generate equidistant seed points along a horizontal line at `y`.
    pub fn line_x(n: usize, y: f32, x_min: f32, x_max: f32) -> Self {
        let points = (0..n)
            .map(|i| {
                let x = if n <= 1 {
                    (x_min + x_max) * 0.5
                } else {
                    x_min + (i as f32 / (n - 1) as f32) * (x_max - x_min)
                };
                [x, y]
            })
            .collect();
        Self { points }
    }

    /// Generate equidistant seed points along a vertical line at `x`.
    pub fn line_y(n: usize, x: f32, y_min: f32, y_max: f32) -> Self {
        let points = (0..n)
            .map(|i| {
                let y = if n <= 1 {
                    (y_min + y_max) * 0.5
                } else {
                    y_min + (i as f32 / (n - 1) as f32) * (y_max - y_min)
                };
                [x, y]
            })
            .collect();
        Self { points }
    }
}

// ---------------------------------------------------------------------------
// GradientField
// ---------------------------------------------------------------------------

/// Gradient of a scalar field computed via central differences.
///
/// Produces a [`VectorField2D`] where `(u, v) = (∂f/∂x, ∂f/∂y)`.
#[derive(Debug, Clone)]
pub struct GradientField;

impl GradientField {
    /// Compute the gradient of `grid` via central differences.
    /// At boundaries, forward/backward differences are used.
    pub fn compute(grid: &HeatmapGrid) -> VectorField2D {
        let dx = grid.dx();
        let dy = grid.dy();
        let mut u = vec![0.0f32; grid.nx * grid.ny];
        let mut v = vec![0.0f32; grid.nx * grid.ny];

        for iy in 0..grid.ny {
            for ix in 0..grid.nx {
                // ∂f/∂x
                let grad_x = if ix == 0 {
                    (grid.get(ix + 1, iy) - grid.get(ix, iy)) / dx
                } else if ix == grid.nx - 1 {
                    (grid.get(ix, iy) - grid.get(ix - 1, iy)) / dx
                } else {
                    (grid.get(ix + 1, iy) - grid.get(ix - 1, iy)) / (2.0 * dx)
                };
                // ∂f/∂y
                let grad_y = if iy == 0 {
                    (grid.get(ix, iy + 1) - grid.get(ix, iy)) / dy
                } else if iy == grid.ny - 1 {
                    (grid.get(ix, iy) - grid.get(ix, iy - 1)) / dy
                } else {
                    (grid.get(ix, iy + 1) - grid.get(ix, iy - 1)) / (2.0 * dy)
                };
                u[iy * grid.nx + ix] = grad_x;
                v[iy * grid.nx + ix] = grad_y;
            }
        }
        VectorField2D::new(
            grid.nx, grid.ny, u, v, grid.x_min, grid.x_max, grid.y_min, grid.y_max,
        )
    }
}

// ---------------------------------------------------------------------------
// HeatmapAnimation
// ---------------------------------------------------------------------------

/// A sequence of [`HeatmapGrid`] frames representing a time-varying scalar field.
#[derive(Debug, Clone)]
pub struct HeatmapAnimation {
    /// Sequence of frames (must all have the same nx, ny, and domain).
    pub frames: Vec<HeatmapGrid>,
    /// Time stamp for each frame (seconds).
    pub times: Vec<f32>,
}

impl HeatmapAnimation {
    /// Create an animation from frames with uniform time spacing `dt`.
    pub fn new_uniform(frames: Vec<HeatmapGrid>, dt: f32) -> Self {
        let times = (0..frames.len()).map(|i| i as f32 * dt).collect();
        Self { frames, times }
    }

    /// Create an animation with explicit per-frame time stamps.
    pub fn new(frames: Vec<HeatmapGrid>, times: Vec<f32>) -> Self {
        assert_eq!(
            frames.len(),
            times.len(),
            "frames and times must have equal length"
        );
        Self { frames, times }
    }

    /// Number of frames.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Total duration (time of last frame).
    pub fn duration(&self) -> f32 {
        self.times.last().copied().unwrap_or(0.0)
    }

    /// Interpolate between frames at time `t`.
    ///
    /// Linearly interpolates cell values between the two nearest frames.
    /// Clamps `t` to `[times[0\], times[last]]`.
    pub fn sample_at(&self, t: f32) -> HeatmapGrid {
        if self.frames.is_empty() {
            panic!("cannot sample empty animation");
        }
        if self.frames.len() == 1 {
            return self.frames[0].clone();
        }
        let t = t.clamp(
            self.times[0],
            *self.times.last().expect("collection should not be empty"),
        );

        // Find surrounding frames
        let idx = self.times.partition_point(|&ti| ti <= t).saturating_sub(1);
        let idx = idx.min(self.frames.len() - 2);

        let t0 = self.times[idx];
        let t1 = self.times[idx + 1];
        let alpha = if (t1 - t0).abs() < 1e-30 {
            0.0
        } else {
            ((t - t0) / (t1 - t0)).clamp(0.0, 1.0)
        };

        let f0 = &self.frames[idx];
        let f1 = &self.frames[idx + 1];
        let data: Vec<f32> = f0
            .data
            .iter()
            .zip(f1.data.iter())
            .map(|(&a, &b)| a * (1.0 - alpha) + b * alpha)
            .collect();

        HeatmapGrid::from_data(f0.nx, f0.ny, data, f0.x_min, f0.x_max, f0.y_min, f0.y_max)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- HeatmapGrid tests ----

    #[test]
    fn test_grid_new_fill_value() {
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 3.125);
        for &v in &grid.data {
            assert!((v - 3.125).abs() < 1e-6);
        }
    }

    #[test]
    fn test_grid_get_set() {
        let mut grid = HeatmapGrid::new(3, 3, 0.0, 1.0, 0.0, 1.0, 0.0);
        grid.set(1, 2, 5.0);
        assert!((grid.get(1, 2) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_grid_out_of_bounds_returns_zero() {
        let grid = HeatmapGrid::new(3, 3, 0.0, 1.0, 0.0, 1.0, 1.0);
        assert_eq!(grid.get(10, 10), 0.0);
    }

    #[test]
    fn test_grid_min_max() {
        let data = vec![1.0f32, 5.0, 3.0, -2.0];
        let grid = HeatmapGrid::from_data(2, 2, data, 0.0, 1.0, 0.0, 1.0);
        assert!((grid.min_val() + 2.0).abs() < 1e-6);
        assert!((grid.max_val() - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_grid_sample_corners() {
        let data = vec![0.0f32, 1.0, 2.0, 3.0]; // 2x2
        let grid = HeatmapGrid::from_data(2, 2, data, 0.0, 1.0, 0.0, 1.0);
        assert!((grid.sample(0.0, 0.0) - 0.0).abs() < 1e-5);
        assert!((grid.sample(1.0, 0.0) - 1.0).abs() < 1e-5);
        assert!((grid.sample(0.0, 1.0) - 2.0).abs() < 1e-5);
        assert!((grid.sample(1.0, 1.0) - 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_grid_sample_bilinear_center() {
        let data = vec![0.0f32, 1.0, 1.0, 2.0]; // 2x2
        let grid = HeatmapGrid::from_data(2, 2, data, 0.0, 1.0, 0.0, 1.0);
        // Center should be average = (0+1+1+2)/4 = 1.0
        let val = grid.sample(0.5, 0.5);
        assert!(
            (val - 1.0).abs() < 1e-5,
            "bilinear center = 1.0, got {}",
            val
        );
    }

    #[test]
    fn test_grid_dx_dy() {
        let grid = HeatmapGrid::new(5, 3, 0.0, 4.0, 0.0, 2.0, 0.0);
        assert!((grid.dx() - 1.0).abs() < 1e-6);
        assert!((grid.dy() - 1.0).abs() < 1e-6);
    }

    // ---- HeatmapRenderer tests ----

    #[test]
    fn test_renderer_output_size() {
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 0.5);
        let renderer = HeatmapRenderer::new(8, 8, Colormap::Viridis);
        let pixels = renderer.render(&grid);
        assert_eq!(pixels.len(), 8 * 8 * 4);
    }

    #[test]
    fn test_renderer_alpha_is_255() {
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 0.5);
        let renderer = HeatmapRenderer::new(4, 4, Colormap::Jet);
        let pixels = renderer.render(&grid);
        for chunk in pixels.chunks(4) {
            assert_eq!(chunk[3], 255, "alpha channel should be 255");
        }
    }

    #[test]
    fn test_renderer_uniform_grid_uniform_color() {
        // Uniform value = 0.5: all pixels should have the same color
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 0.5);
        let renderer = HeatmapRenderer::new(4, 4, Colormap::Viridis);
        let pixels = renderer.render(&grid);
        let first = [pixels[0], pixels[1], pixels[2]];
        for chunk in pixels.chunks(4) {
            assert_eq!([chunk[0], chunk[1], chunk[2]], first);
        }
    }

    #[test]
    fn test_renderer_with_range() {
        let data = vec![0.0f32, 1.0, 0.0, 1.0];
        let grid = HeatmapGrid::from_data(2, 2, data, 0.0, 1.0, 0.0, 1.0);
        let renderer = HeatmapRenderer::new(2, 2, Colormap::Jet).with_range(0.0, 1.0);
        let pixels = renderer.render(&grid);
        assert_eq!(pixels.len(), 2 * 2 * 4);
    }

    // ---- ContourLine tests ----

    #[test]
    fn test_contour_no_crossings() {
        // All values above level → no contour
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 2.0);
        let contour = ContourLine::extract(&grid, 1.0);
        assert_eq!(contour.segment_count(), 0);
    }

    #[test]
    fn test_contour_all_below_level() {
        // All values below level → no contour
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 0.0);
        let contour = ContourLine::extract(&grid, 1.0);
        assert_eq!(contour.segment_count(), 0);
    }

    #[test]
    fn test_contour_simple_vertical_crossing() {
        // Left half < 0.5, right half > 0.5: vertical contour through middle
        let mut grid = HeatmapGrid::new(4, 4, 0.0, 3.0, 0.0, 3.0, 0.0);
        for iy in 0..4 {
            grid.set(0, iy, 0.0);
            grid.set(1, iy, 0.0);
            grid.set(2, iy, 1.0);
            grid.set(3, iy, 1.0);
        }
        let contour = ContourLine::extract(&grid, 0.5);
        assert!(
            contour.segment_count() > 0,
            "should have vertical contour segments"
        );
    }

    #[test]
    fn test_contour_arc_length_positive() {
        let mut grid = HeatmapGrid::new(4, 4, 0.0, 3.0, 0.0, 3.0, 0.0);
        for iy in 0..4 {
            grid.set(0, iy, 0.0);
            grid.set(1, iy, 0.0);
            grid.set(2, iy, 1.0);
            grid.set(3, iy, 1.0);
        }
        let contour = ContourLine::extract(&grid, 0.5);
        if contour.segment_count() > 0 {
            assert!(contour.arc_length() > 0.0);
        }
    }

    // ---- ContourSet tests ----

    #[test]
    fn test_contour_set_n_levels() {
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 0.5);
        let cs = ContourSet::extract_n(&grid, 5, 0.0, 1.0);
        assert_eq!(cs.contours.len(), 5);
    }

    #[test]
    fn test_contour_set_empty_levels() {
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 0.5);
        let cs = ContourSet::extract_n(&grid, 0, 0.0, 1.0);
        assert!(cs.contours.is_empty());
    }

    #[test]
    fn test_contour_set_closest_contour() {
        let mut grid = HeatmapGrid::new(4, 4, 0.0, 3.0, 0.0, 3.0, 0.0);
        for iy in 0..4 {
            grid.set(2, iy, 1.0);
            grid.set(3, iy, 1.0);
        }
        let levels = vec![0.1, 0.5, 0.9];
        let cs = ContourSet::extract(&grid, &levels);
        let idx = cs.closest_contour_idx([1.5, 1.5]);
        assert!(idx.is_some());
    }

    // ---- VectorField2D tests ----

    #[test]
    fn test_vector_field_get() {
        let u = vec![1.0f32, 2.0, 3.0, 4.0];
        let v = vec![0.0f32; 4];
        let field = VectorField2D::new(2, 2, u, v, 0.0, 1.0, 0.0, 1.0);
        let (u00, _) = field.get(0, 0);
        assert!((u00 - 1.0).abs() < 1e-6);
        let (u11, _) = field.get(1, 1);
        assert!((u11 - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_vector_field_magnitude_grid() {
        let u = vec![3.0f32, 0.0, 0.0, 0.0];
        let v = vec![4.0f32, 0.0, 0.0, 0.0];
        let field = VectorField2D::new(2, 2, u, v, 0.0, 1.0, 0.0, 1.0);
        let mag = field.magnitude_grid();
        assert!((mag.get(0, 0) - 5.0).abs() < 1e-5, "magnitude should be 5");
    }

    #[test]
    fn test_vector_field_divergence_uniform() {
        // Uniform field: div = 0
        let u = vec![1.0f32; 9];
        let v = vec![0.0f32; 9];
        let field = VectorField2D::new(3, 3, u, v, 0.0, 1.0, 0.0, 1.0);
        let div = field.divergence();
        let center = div.get(1, 1);
        assert!(
            center.abs() < 1e-5,
            "divergence of uniform field should be 0"
        );
    }

    #[test]
    fn test_vector_field_curl_irrotational() {
        // Field (x, y) has curl = ∂y/∂x - ∂x/∂y = 0 - 0 = 0
        let mut u = vec![0.0f32; 9];
        let mut v = vec![0.0f32; 9];
        for iy in 0..3usize {
            for ix in 0..3usize {
                u[iy * 3 + ix] = ix as f32;
                v[iy * 3 + ix] = iy as f32;
            }
        }
        let field = VectorField2D::new(3, 3, u, v, 0.0, 2.0, 0.0, 2.0);
        let curl = field.curl();
        let center = curl.get(1, 1);
        assert!(
            center.abs() < 0.1,
            "irrotational field curl should be ~0, got {}",
            center
        );
    }

    #[test]
    fn test_vector_field_sample_at_corner() {
        let u = vec![1.0f32, 2.0, 3.0, 4.0];
        let v = vec![0.0f32; 4];
        let field = VectorField2D::new(2, 2, u, v, 0.0, 1.0, 0.0, 1.0);
        let (us, _) = field.sample(0.0, 0.0);
        assert!((us - 1.0).abs() < 1e-5);
        let (us2, _) = field.sample(1.0, 1.0);
        assert!((us2 - 4.0).abs() < 1e-5);
    }

    // ---- ArrowGlyph tests ----

    #[test]
    fn test_arrow_glyph_tip_direction() {
        let arrow = ArrowGlyph::new([0.0, 0.0], [1.0, 0.0], 2.0, 10.0);
        let tip = arrow.tip();
        assert!((tip[0] - 2.0).abs() < 1e-5);
        assert!(tip[1].abs() < 1e-5);
    }

    #[test]
    fn test_arrow_glyph_clamp_length() {
        let arrow = ArrowGlyph::new([0.0, 0.0], [1.0, 0.0], 100.0, 5.0);
        assert!((arrow.drawn_length() - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_arrow_glyph_zero_direction() {
        let arrow = ArrowGlyph::new([1.0, 2.0], [0.0, 0.0], 1.0, 5.0);
        let tip = arrow.tip();
        // Zero direction → tip = position
        assert!((tip[0] - 1.0).abs() < 1e-6);
        assert!((tip[1] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_arrow_glyph_as_segment() {
        let arrow = ArrowGlyph::new([0.0, 0.0], [0.0, 1.0], 3.0, 10.0);
        let seg = arrow.as_segment();
        assert_eq!(seg[0], [0.0f32, 0.0]);
        assert!((seg[1][1] - 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_arrow_glyphs_from_field_count() {
        let u = vec![1.0f32; 9];
        let v = vec![0.0f32; 9];
        let field = VectorField2D::new(3, 3, u, v, 0.0, 1.0, 0.0, 1.0);
        let glyphs = arrow_glyphs_from_field(&field, 1, 1, 1.0);
        assert_eq!(glyphs.len(), 9);
    }

    // ---- StreamlineIntegrator tests ----

    #[test]
    fn test_streamline_uniform_field_direction() {
        // Uniform field (1, 0): streamline should go right
        let u = vec![1.0f32; 16];
        let v = vec![0.0f32; 16];
        let field = VectorField2D::new(4, 4, u, v, 0.0, 3.0, 0.0, 3.0);
        let integrator = StreamlineIntegrator::new(0.1, 10, 1e-6);
        let sl = integrator.integrate(&field, [0.5, 1.5]);
        // All points should have increasing x
        for i in 1..sl.points.len() {
            assert!(sl.points[i][0] >= sl.points[i - 1][0]);
        }
    }

    #[test]
    fn test_streamline_stops_at_boundary() {
        let u = vec![1.0f32; 16]; // all going right
        let v = vec![0.0f32; 16];
        let field = VectorField2D::new(4, 4, u, v, 0.0, 1.0, 0.0, 1.0);
        let integrator = StreamlineIntegrator::new(0.1, 100, 1e-6);
        let sl = integrator.integrate(&field, [0.5, 0.5]);
        // Should stop when x > 1.0
        if let Some(&last) = sl.points.last() {
            assert!(last[0] <= 1.0 + 1e-5, "should not exceed x_max");
        }
    }

    #[test]
    fn test_streamline_min_velocity_stops() {
        let u = vec![0.0f32; 9]; // zero velocity everywhere
        let v = vec![0.0f32; 9];
        let field = VectorField2D::new(3, 3, u, v, 0.0, 1.0, 0.0, 1.0);
        let integrator = StreamlineIntegrator::new(0.1, 50, 0.1);
        let sl = integrator.integrate(&field, [0.5, 0.5]);
        // Should stop immediately due to low velocity
        assert!(sl.points.len() <= 2);
    }

    #[test]
    fn test_streamline_seed_is_first_point() {
        let u = vec![1.0f32; 16];
        let v = vec![0.0f32; 16];
        let field = VectorField2D::new(4, 4, u, v, 0.0, 3.0, 0.0, 3.0);
        let integrator = StreamlineIntegrator::new(0.1, 5, 1e-6);
        let sl = integrator.integrate(&field, [1.0, 1.5]);
        assert!((sl.points[0][0] - 1.0).abs() < 1e-6);
        assert!((sl.points[0][1] - 1.5).abs() < 1e-6);
    }

    // ---- SeedPoints tests ----

    #[test]
    fn test_seed_points_uniform_count() {
        let seeds = SeedPoints::uniform(3, 4, 0.0, 1.0, 0.0, 1.0);
        assert_eq!(seeds.points.len(), 12);
    }

    #[test]
    fn test_seed_points_uniform_corners() {
        let seeds = SeedPoints::uniform(2, 2, 0.0, 1.0, 0.0, 1.0);
        assert_eq!(seeds.points.len(), 4);
        assert!((seeds.points[0][0]).abs() < 1e-6);
        assert!((seeds.points[3][0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_seed_points_random_in_bounds() {
        let seeds = SeedPoints::random(20, 0.0, 5.0, -1.0, 1.0, 42);
        assert_eq!(seeds.points.len(), 20);
        for &[x, y] in &seeds.points {
            assert!((0.0..=5.0).contains(&x), "x={} out of bounds", x);
            assert!((-1.0..=1.0).contains(&y), "y={} out of bounds", y);
        }
    }

    #[test]
    fn test_seed_points_line_x() {
        let seeds = SeedPoints::line_x(5, 2.0, 0.0, 4.0);
        assert_eq!(seeds.points.len(), 5);
        for &[_x, y] in &seeds.points {
            assert!((y - 2.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_seed_points_line_y() {
        let seeds = SeedPoints::line_y(5, 3.0, 0.0, 4.0);
        assert_eq!(seeds.points.len(), 5);
        for &[x, _y] in &seeds.points {
            assert!((x - 3.0).abs() < 1e-6);
        }
    }

    // ---- GradientField tests ----

    #[test]
    fn test_gradient_linear_field_x() {
        // f(x, y) = x → ∂f/∂x = 1, ∂f/∂y = 0
        let mut grid = HeatmapGrid::new(3, 3, 0.0, 2.0, 0.0, 2.0, 0.0);
        for iy in 0..3 {
            for ix in 0..3 {
                grid.set(ix, iy, ix as f32);
            }
        }
        let grad = GradientField::compute(&grid);
        let (gx, gy) = grad.get(1, 1); // central cell
        assert!((gx - 1.0).abs() < 0.1, "∂f/∂x should be ~1, got {}", gx);
        assert!(gy.abs() < 0.1, "∂f/∂y should be ~0, got {}", gy);
    }

    #[test]
    fn test_gradient_linear_field_y() {
        // f(x, y) = y → ∂f/∂x = 0, ∂f/∂y = 1
        let mut grid = HeatmapGrid::new(3, 3, 0.0, 2.0, 0.0, 2.0, 0.0);
        for iy in 0..3 {
            for ix in 0..3 {
                grid.set(ix, iy, iy as f32);
            }
        }
        let grad = GradientField::compute(&grid);
        let (gx, gy) = grad.get(1, 1);
        assert!(gx.abs() < 0.1, "∂f/∂x should be ~0, got {}", gx);
        assert!((gy - 1.0).abs() < 0.1, "∂f/∂y should be ~1, got {}", gy);
    }

    #[test]
    fn test_gradient_constant_field_zero_gradient() {
        let grid = HeatmapGrid::new(4, 4, 0.0, 1.0, 0.0, 1.0, 7.0);
        let grad = GradientField::compute(&grid);
        for iy in 0..4 {
            for ix in 0..4 {
                let (gx, gy) = grad.get(ix, iy);
                assert!(gx.abs() < 1e-5, "constant field gradient x should be 0");
                assert!(gy.abs() < 1e-5, "constant field gradient y should be 0");
            }
        }
    }

    // ---- HeatmapAnimation tests ----

    #[test]
    fn test_animation_frame_count() {
        let frames: Vec<HeatmapGrid> = (0..5)
            .map(|_| HeatmapGrid::new(3, 3, 0.0, 1.0, 0.0, 1.0, 0.0))
            .collect();
        let anim = HeatmapAnimation::new_uniform(frames, 0.1);
        assert_eq!(anim.frame_count(), 5);
    }

    #[test]
    fn test_animation_duration() {
        let frames: Vec<HeatmapGrid> = (0..5)
            .map(|_| HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 0.0))
            .collect();
        let anim = HeatmapAnimation::new_uniform(frames, 0.1);
        assert!((anim.duration() - 0.4).abs() < 1e-6);
    }

    #[test]
    fn test_animation_sample_at_frame_0() {
        let frame0 = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 1.0);
        let frame1 = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 3.0);
        let anim = HeatmapAnimation::new_uniform(vec![frame0, frame1], 1.0);
        let s = anim.sample_at(0.0);
        assert!((s.data[0] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_animation_sample_at_midpoint() {
        let frame0 = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 0.0);
        let frame1 = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 2.0);
        let anim = HeatmapAnimation::new_uniform(vec![frame0, frame1], 1.0);
        let s = anim.sample_at(0.5);
        assert!(
            (s.data[0] - 1.0).abs() < 1e-5,
            "midpoint interpolation should give 1.0"
        );
    }

    #[test]
    fn test_animation_sample_clamps_below() {
        let frame0 = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 5.0);
        let frame1 = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 10.0);
        let anim = HeatmapAnimation::new_uniform(vec![frame0, frame1], 1.0);
        let s = anim.sample_at(-100.0); // clamped to t=0
        assert!((s.data[0] - 5.0).abs() < 1e-5);
    }

    #[test]
    fn test_animation_single_frame() {
        let frame = HeatmapGrid::new(2, 2, 0.0, 1.0, 0.0, 1.0, 42.0);
        let anim = HeatmapAnimation::new_uniform(vec![frame], 1.0);
        let s = anim.sample_at(0.0);
        assert!((s.data[0] - 42.0).abs() < 1e-5);
    }

    // ---- Point to segment distance test ----

    #[test]
    fn test_point_to_segment_dist_on_segment() {
        let d = point_to_segment_dist(1.0, 0.0, 0.0, 0.0, 2.0, 0.0);
        assert!(d.abs() < 1e-5);
    }

    #[test]
    fn test_point_to_segment_dist_off_segment() {
        let d = point_to_segment_dist(1.0, 1.0, 0.0, 0.0, 2.0, 0.0);
        assert!((d - 1.0).abs() < 1e-5);
    }
}
