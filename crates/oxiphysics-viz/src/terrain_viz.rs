#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Terrain and landscape visualization for the OxiPhysics engine.
//!
//! This module provides:
//!
//! - **[`Heightfield`]** — regular grid terrain with LOD mesh generation
//! - **[`TextureSplat`]** — multi-layer texture blending by slope and altitude
//! - **[`ContourExtractor`]** — marching-squares contour-line extraction
//! - **[`WaterSurface`]** — flat water level with shoreline computation
//! - **[`ErosionMap`]** — erosion overlay visualization
//! - **[`VegetationDensityMap`]** — vegetation density heat-map
//! - **[`SlopeAspectMap`]** — slope angle and cardinal aspect computation
//! - **[`TerrainShadow`]** — ray-marching terrain shadow casting
//! - **[`CrossSection`]** — terrain profile along an arbitrary polyline
//! - **[`ElevationColorRamp`]** — configurable elevation-to-color mapping

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// Heightfield — core terrain data structure
// ─────────────────────────────────────────────────────────────────────────────

/// A regular-grid heightfield with `nx × ny` sample points.
///
/// Vertices are stored in row-major order: the sample at column `i`, row `j`
/// is located at index `j * nx + i`.  The world-space extents are
/// `[0, (nx-1)·dx] × [0, (ny-1)·dy]` with height along the Y axis.
#[derive(Debug, Clone)]
pub struct Heightfield {
    /// Number of columns (X-axis sample count).
    pub nx: usize,
    /// Number of rows (Z-axis sample count).
    pub ny: usize,
    /// Column spacing in world units.
    pub dx: f64,
    /// Row spacing in world units.
    pub dz: f64,
    /// Height values (Y), length `nx * ny`.
    pub heights: Vec<f64>,
}

impl Heightfield {
    /// Construct a zero-filled heightfield.
    pub fn new(nx: usize, ny: usize, dx: f64, dz: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            dz,
            heights: vec![0.0; nx * ny],
        }
    }

    /// Return the height at column `i`, row `j`, or `None` if out of bounds.
    pub fn height(&self, i: usize, j: usize) -> Option<f64> {
        if i >= self.nx || j >= self.ny {
            return None;
        }
        Some(self.heights[j * self.nx + i])
    }

    /// Set the height at column `i`, row `j`.  Returns `false` if out of bounds.
    pub fn set_height(&mut self, i: usize, j: usize, h: f64) -> bool {
        if i >= self.nx || j >= self.ny {
            return false;
        }
        self.heights[j * self.nx + i] = h;
        true
    }

    /// Bilinear interpolation of height at continuous grid coordinates `(u, v)`.
    ///
    /// `u ∈ [0, nx-1]`, `v ∈ [0, ny-1]`.  Returns `None` when out of range.
    pub fn sample(&self, u: f64, v: f64) -> Option<f64> {
        if u < 0.0 || v < 0.0 {
            return None;
        }
        let i0 = u.floor() as usize;
        let j0 = v.floor() as usize;
        if i0 + 1 >= self.nx || j0 + 1 >= self.ny {
            return None;
        }
        let fu = u - i0 as f64;
        let fv = v - j0 as f64;
        let h00 = self.heights[j0 * self.nx + i0];
        let h10 = self.heights[j0 * self.nx + i0 + 1];
        let h01 = self.heights[(j0 + 1) * self.nx + i0];
        let h11 = self.heights[(j0 + 1) * self.nx + i0 + 1];
        Some(
            h00 * (1.0 - fu) * (1.0 - fv)
                + h10 * fu * (1.0 - fv)
                + h01 * (1.0 - fu) * fv
                + h11 * fu * fv,
        )
    }

    /// World-space position `[x, y, z]` for grid sample `(i, j)`.
    pub fn world_pos(&self, i: usize, j: usize) -> Option<[f64; 3]> {
        let h = self.height(i, j)?;
        Some([i as f64 * self.dx, h, j as f64 * self.dz])
    }

    /// Minimum height in the heightfield.
    pub fn min_height(&self) -> f64 {
        self.heights.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Maximum height in the heightfield.
    pub fn max_height(&self) -> f64 {
        self.heights
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Height range `(min, max)`.
    pub fn height_range(&self) -> (f64, f64) {
        (self.min_height(), self.max_height())
    }

    /// Finite-difference surface normal at `(i, j)` (un-normalised).
    ///
    /// The normal points in the direction of steepest ascent (upward component
    /// is always positive for flat or gentle terrain).
    pub fn normal(&self, i: usize, j: usize) -> [f64; 3] {
        let h = |ci: usize, cj: usize| self.heights[cj * self.nx + ci];
        let ci = i.clamp(1, self.nx - 2);
        let cj = j.clamp(1, self.ny - 2);
        let dhdx = (h(ci + 1, cj) - h(ci - 1, cj)) / (2.0 * self.dx);
        let dhdz = (h(ci, cj + 1) - h(ci, cj - 1)) / (2.0 * self.dz);
        let len = (dhdx * dhdx + 1.0 + dhdz * dhdz).sqrt();
        [-dhdx / len, 1.0 / len, -dhdz / len]
    }

    /// Slope angle in radians at `(i, j)` (angle with the horizontal plane).
    pub fn slope_rad(&self, i: usize, j: usize) -> f64 {
        let h = |ci: usize, cj: usize| self.heights[cj * self.nx + ci];
        let ci = i.clamp(1, self.nx - 2);
        let cj = j.clamp(1, self.ny - 2);
        let dhdx = (h(ci + 1, cj) - h(ci - 1, cj)) / (2.0 * self.dx);
        let dhdz = (h(ci, cj + 1) - h(ci, cj - 1)) / (2.0 * self.dz);
        (dhdx * dhdx + dhdz * dhdz).sqrt().atan()
    }

    /// Generate an LOD triangle mesh as `(vertices, indices)`.
    ///
    /// `step` controls the spacing between sampled vertices.  `step = 1` gives
    /// the full-resolution mesh; larger values reduce triangle count.  The
    /// returned vertex positions are `[x, y, z]` world-space triples.
    pub fn lod_mesh(&self, step: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let s = step.max(1);
        let cols: Vec<usize> = (0..self.nx).step_by(s).collect();
        let rows: Vec<usize> = (0..self.ny).step_by(s).collect();
        let mut verts: Vec<[f64; 3]> = Vec::with_capacity(cols.len() * rows.len());
        for &j in &rows {
            for &i in &cols {
                let h = self.heights[j * self.nx + i];
                verts.push([i as f64 * self.dx, h, j as f64 * self.dz]);
            }
        }
        let cw = cols.len();
        let rh = rows.len();
        let mut tris: Vec<[u32; 3]> = Vec::with_capacity((cw - 1) * (rh - 1) * 2);
        for rj in 0..rh.saturating_sub(1) {
            for ci in 0..cw.saturating_sub(1) {
                let tl = (rj * cw + ci) as u32;
                let tr = tl + 1;
                let bl = tl + cw as u32;
                let br = bl + 1;
                tris.push([tl, bl, tr]);
                tris.push([tr, bl, br]);
            }
        }
        (verts, tris)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TextureSplat — multi-layer texture blending
// ─────────────────────────────────────────────────────────────────────────────

/// A single texture layer in a terrain splat map.
#[derive(Debug, Clone)]
pub struct SplatLayer {
    /// Human-readable name (e.g., `"rock"`, `"grass"`, `"snow"`).
    pub name: String,
    /// Minimum slope angle (radians) at which this layer begins to appear.
    pub slope_min_rad: f64,
    /// Maximum slope angle (radians) beyond which this layer is fully blended out.
    pub slope_max_rad: f64,
    /// Minimum altitude (world-Y) at which this layer is active.
    pub altitude_min: f64,
    /// Maximum altitude at which this layer is active.
    pub altitude_max: f64,
    /// Base RGBA colour for this layer.
    pub color: [f32; 4],
}

impl SplatLayer {
    /// Create a new splat layer.
    pub fn new(
        name: impl Into<String>,
        slope_min_rad: f64,
        slope_max_rad: f64,
        altitude_min: f64,
        altitude_max: f64,
        color: [f32; 4],
    ) -> Self {
        Self {
            name: name.into(),
            slope_min_rad,
            slope_max_rad,
            altitude_min,
            altitude_max,
            color,
        }
    }

    /// Blend weight for a surface point at given slope (rad) and altitude.
    ///
    /// Returns a value in `[0, 1]` where 1 means fully active.
    pub fn weight(&self, slope_rad: f64, altitude: f64) -> f32 {
        let slope_w = if self.slope_max_rad <= self.slope_min_rad {
            if slope_rad >= self.slope_min_rad {
                1.0
            } else {
                0.0
            }
        } else {
            let t = (slope_rad - self.slope_min_rad) / (self.slope_max_rad - self.slope_min_rad);
            t.clamp(0.0, 1.0)
        };
        let alt_w = if self.altitude_max <= self.altitude_min {
            if altitude >= self.altitude_min {
                1.0
            } else {
                0.0
            }
        } else {
            let t = (altitude - self.altitude_min) / (self.altitude_max - self.altitude_min);
            // bell: peaks at 0.5 of the altitude range
            let t2 = t.clamp(0.0, 1.0);
            (1.0 - (2.0 * t2 - 1.0).abs()).max(0.0)
        };
        (slope_w as f32) * (alt_w as f32)
    }
}

/// A collection of [`SplatLayer`]s that together define terrain surface appearance.
#[derive(Debug, Clone, Default)]
pub struct TextureSplat {
    /// The layers; typically 3–6 for real terrain.
    pub layers: Vec<SplatLayer>,
}

impl TextureSplat {
    /// Create an empty splat.
    pub fn new() -> Self {
        Self { layers: Vec::new() }
    }

    /// Add a layer.
    pub fn add_layer(&mut self, layer: SplatLayer) {
        self.layers.push(layer);
    }

    /// Blend all layers at a given surface point and return a normalised RGBA colour.
    ///
    /// Weights are normalised so they sum to 1.  If all weights are zero the
    /// first layer (or `[0.5, 0.5, 0.5, 1.0]` if empty) is returned.
    pub fn blend(&self, slope_rad: f64, altitude: f64) -> [f32; 4] {
        if self.layers.is_empty() {
            return [0.5, 0.5, 0.5, 1.0];
        }
        let weights: Vec<f32> = self
            .layers
            .iter()
            .map(|l| l.weight(slope_rad, altitude))
            .collect();
        let total: f32 = weights.iter().sum();
        if total < 1e-6 {
            return self.layers[0].color;
        }
        let mut out = [0.0f32; 4];
        for (w, l) in weights.iter().zip(self.layers.iter()) {
            let wn = w / total;
            for k in 0..4 {
                out[k] += wn * l.color[k];
            }
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ContourExtractor — marching-squares iso-contour lines
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D line segment in world space (XZ plane).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContourSegment {
    /// First endpoint `[x, z]`.
    pub a: [f64; 2],
    /// Second endpoint `[x, z]`.
    pub b: [f64; 2],
    /// Elevation of the iso-contour that generated this segment.
    pub elevation: f64,
}

/// Extracts iso-contour lines from a [`Heightfield`] at specified elevations.
#[derive(Debug, Clone)]
pub struct ContourExtractor<'a> {
    /// Reference to the terrain.
    pub field: &'a Heightfield,
}

impl<'a> ContourExtractor<'a> {
    /// Create a new extractor for the given heightfield.
    pub fn new(field: &'a Heightfield) -> Self {
        Self { field }
    }

    /// Extract contour segments for a single elevation `h`.
    ///
    /// Uses the marching-squares algorithm on every 2×2 cell.
    pub fn extract_contour(&self, h: f64) -> Vec<ContourSegment> {
        let hf = self.field;
        let mut segs = Vec::new();
        for j in 0..hf.ny.saturating_sub(1) {
            for i in 0..hf.nx.saturating_sub(1) {
                let h00 = hf.heights[j * hf.nx + i];
                let h10 = hf.heights[j * hf.nx + i + 1];
                let h01 = hf.heights[(j + 1) * hf.nx + i];
                let h11 = hf.heights[(j + 1) * hf.nx + i + 1];

                let above = |v: f64| v >= h;
                let code = (above(h00) as u8)
                    | ((above(h10) as u8) << 1)
                    | ((above(h11) as u8) << 2)
                    | ((above(h01) as u8) << 3);

                // Skip uniform cells
                if code == 0 || code == 15 {
                    continue;
                }

                // Linear interpolation along an edge
                let interp = |va: f64, vb: f64, xa: f64, xb: f64| -> f64 {
                    if (vb - va).abs() < 1e-12 {
                        (xa + xb) * 0.5
                    } else {
                        xa + (h - va) / (vb - va) * (xb - xa)
                    }
                };

                let x0 = i as f64 * hf.dx;
                let x1 = (i + 1) as f64 * hf.dx;
                let z0 = j as f64 * hf.dz;
                let z1 = (j + 1) as f64 * hf.dz;

                // Bottom edge (j=j, i→i+1)
                let bottom = || [interp(h00, h10, x0, x1), z0];
                // Right edge (i=i+1, j→j+1)
                let right = || [x1, interp(h10, h11, z0, z1)];
                // Top edge (j=j+1, i→i+1)
                let top = || [interp(h01, h11, x0, x1), z1];
                // Left edge (i=i, j→j+1)
                let left = || [x0, interp(h00, h01, z0, z1)];

                let push = |a: [f64; 2], b: [f64; 2], segs: &mut Vec<ContourSegment>| {
                    segs.push(ContourSegment { a, b, elevation: h });
                };

                match code {
                    1 | 14 => push(bottom(), left(), &mut segs),
                    2 | 13 => push(bottom(), right(), &mut segs),
                    3 | 12 => push(left(), right(), &mut segs),
                    4 | 11 => push(top(), right(), &mut segs),
                    6 | 9 => push(bottom(), top(), &mut segs),
                    7 | 8 => push(top(), left(), &mut segs),
                    5 => {
                        push(bottom(), left(), &mut segs);
                        push(top(), right(), &mut segs);
                    }
                    10 => {
                        push(bottom(), right(), &mut segs);
                        push(top(), left(), &mut segs);
                    }
                    _ => {}
                }
            }
        }
        segs
    }

    /// Extract contours at regularly spaced elevations.
    ///
    /// `interval` is the spacing between contour levels.  Returns all segments
    /// across all levels.
    pub fn extract_all(&self, interval: f64) -> Vec<ContourSegment> {
        let (h_min, h_max) = self.field.height_range();
        if interval <= 0.0 || h_min >= h_max {
            return vec![];
        }
        let mut all = Vec::new();
        let mut h = (h_min / interval).ceil() * interval;
        while h <= h_max {
            all.extend(self.extract_contour(h));
            h += interval;
        }
        all
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WaterSurface — flat water plane with shoreline detection
// ─────────────────────────────────────────────────────────────────────────────

/// A flat water surface over a [`Heightfield`].
#[derive(Debug, Clone)]
pub struct WaterSurface<'a> {
    /// The underlying terrain.
    pub terrain: &'a Heightfield,
    /// Water level elevation in world-Y units.
    pub level: f64,
    /// Water color `[r, g, b, a]`.
    pub color: [f32; 4],
}

impl<'a> WaterSurface<'a> {
    /// Create a new water surface at the given level.
    pub fn new(terrain: &'a Heightfield, level: f64, color: [f32; 4]) -> Self {
        Self {
            terrain,
            level,
            color,
        }
    }

    /// Return `true` if grid cell `(i, j)` is underwater.
    pub fn is_underwater(&self, i: usize, j: usize) -> bool {
        self.terrain.height(i, j).is_some_and(|h| h < self.level)
    }

    /// Fraction `[0, 1]` of heightfield cells that are underwater.
    pub fn flood_fraction(&self) -> f64 {
        let total = self.terrain.nx * self.terrain.ny;
        if total == 0 {
            return 0.0;
        }
        let flooded = (0..self.terrain.ny)
            .flat_map(|j| (0..self.terrain.nx).map(move |i| (i, j)))
            .filter(|&(i, j)| self.is_underwater(i, j))
            .count();
        flooded as f64 / total as f64
    }

    /// Extract shoreline contour segments at the water level.
    pub fn shoreline(&self) -> Vec<ContourSegment> {
        ContourExtractor::new(self.terrain).extract_contour(self.level)
    }

    /// Build a water quad mesh covering the entire terrain AABB at `level`.
    ///
    /// Returns `(vertices [x,y,z], indices)` for a single quad rendered at the
    /// water elevation.
    pub fn water_quad(&self) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let x_max = (self.terrain.nx.saturating_sub(1)) as f64 * self.terrain.dx;
        let z_max = (self.terrain.ny.saturating_sub(1)) as f64 * self.terrain.dz;
        let y = self.level;
        let verts = vec![
            [0.0, y, 0.0],
            [x_max, y, 0.0],
            [x_max, y, z_max],
            [0.0, y, z_max],
        ];
        let tris = vec![[0, 1, 2], [0, 2, 3]];
        (verts, tris)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ErosionMap — erosion overlay visualization
// ─────────────────────────────────────────────────────────────────────────────

/// Erosion intensity record for a single terrain cell.
#[derive(Debug, Clone, Copy)]
pub struct ErosionCell {
    /// Total sediment deposited (positive) or removed (negative) in metres.
    pub net_sediment_m: f64,
    /// Cumulative water flux in m²/s.
    pub water_flux: f64,
}

/// An erosion map overlay for a [`Heightfield`], with one record per cell.
#[derive(Debug, Clone)]
pub struct ErosionMap {
    /// Number of columns (matches the terrain).
    pub nx: usize,
    /// Number of rows (matches the terrain).
    pub ny: usize,
    /// Per-cell erosion data, row-major.
    pub cells: Vec<ErosionCell>,
}

impl ErosionMap {
    /// Create an erosion map with all cells zero-initialised.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            cells: vec![
                ErosionCell {
                    net_sediment_m: 0.0,
                    water_flux: 0.0
                };
                nx * ny
            ],
        }
    }

    /// Get a reference to the cell at `(i, j)`.
    pub fn cell(&self, i: usize, j: usize) -> Option<&ErosionCell> {
        if i >= self.nx || j >= self.ny {
            return None;
        }
        Some(&self.cells[j * self.nx + i])
    }

    /// Set the cell at `(i, j)`.
    pub fn set_cell(&mut self, i: usize, j: usize, c: ErosionCell) -> bool {
        if i >= self.nx || j >= self.ny {
            return false;
        }
        self.cells[j * self.nx + i] = c;
        true
    }

    /// Map `net_sediment_m` to an RGBA color.
    ///
    /// Negative (erosion) → red; zero → white; positive (deposition) → blue.
    pub fn color_for_cell(&self, i: usize, j: usize) -> [f32; 4] {
        let c = match self.cell(i, j) {
            Some(c) => c,
            None => return [0.0, 0.0, 0.0, 0.0],
        };
        if c.net_sediment_m < 0.0 {
            let t = (c.net_sediment_m.abs() / 10.0).min(1.0) as f32;
            [t, 1.0 - t, 1.0 - t, 1.0]
        } else {
            let t = (c.net_sediment_m / 10.0).min(1.0) as f32;
            [1.0 - t, 1.0 - t, t, 1.0]
        }
    }

    /// Maximum absolute `net_sediment_m` across the map.
    pub fn max_abs_sediment(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| c.net_sediment_m.abs())
            .fold(0.0_f64, f64::max)
    }

    /// Maximum `water_flux` across the map.
    pub fn max_flux(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| c.water_flux)
            .fold(0.0_f64, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VegetationDensityMap — vegetation coverage visualization
// ─────────────────────────────────────────────────────────────────────────────

/// Vegetation density at a single terrain cell, `[0, 1]`.
#[derive(Debug, Clone, Copy)]
pub struct VegetationCell {
    /// Normalised plant coverage fraction `[0, 1]`.
    pub density: f32,
    /// Dominant plant-type index (user-defined).
    pub plant_type: u8,
}

/// A vegetation density map overlay.
#[derive(Debug, Clone)]
pub struct VegetationDensityMap {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// Per-cell vegetation data.
    pub cells: Vec<VegetationCell>,
}

impl VegetationDensityMap {
    /// Create a map with zero density.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            cells: vec![
                VegetationCell {
                    density: 0.0,
                    plant_type: 0
                };
                nx * ny
            ],
        }
    }

    /// Get the cell at `(i, j)`.
    pub fn cell(&self, i: usize, j: usize) -> Option<&VegetationCell> {
        if i >= self.nx || j >= self.ny {
            return None;
        }
        Some(&self.cells[j * self.nx + i])
    }

    /// Set the cell at `(i, j)`.
    pub fn set_cell(&mut self, i: usize, j: usize, c: VegetationCell) -> bool {
        if i >= self.nx || j >= self.ny {
            return false;
        }
        self.cells[j * self.nx + i] = c;
        true
    }

    /// RGBA color for a vegetation cell (green gradient by density).
    pub fn color_for_cell(&self, i: usize, j: usize) -> [f32; 4] {
        let d = self.cell(i, j).map_or(0.0, |c| c.density);
        [0.0, d, 0.0, 1.0]
    }

    /// Mean vegetation density across the whole map.
    pub fn mean_density(&self) -> f32 {
        if self.cells.is_empty() {
            return 0.0;
        }
        self.cells.iter().map(|c| c.density).sum::<f32>() / self.cells.len() as f32
    }

    /// Populate density based on slope and altitude from a terrain.
    ///
    /// Simple rule: vegetation grows where slope < `max_slope_rad` and
    /// altitude is in `[alt_min, alt_max]`.
    pub fn populate_from_terrain(
        &mut self,
        terrain: &Heightfield,
        max_slope_rad: f64,
        alt_min: f64,
        alt_max: f64,
    ) {
        for j in 0..self.ny.min(terrain.ny) {
            for i in 0..self.nx.min(terrain.nx) {
                let h = terrain.heights[j * terrain.nx + i];
                let slope = terrain.slope_rad(i, j);
                let density = if slope < max_slope_rad && h >= alt_min && h <= alt_max {
                    let slope_factor = (1.0 - slope / max_slope_rad).clamp(0.0, 1.0);
                    let alt_factor = {
                        let t = ((h - alt_min) / (alt_max - alt_min)).clamp(0.0, 1.0);
                        (1.0 - (2.0 * t - 1.0).abs()).max(0.0)
                    };
                    (slope_factor * alt_factor) as f32
                } else {
                    0.0
                };
                self.cells[j * self.nx + i].density = density;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SlopeAspectMap — slope angle and cardinal aspect
// ─────────────────────────────────────────────────────────────────────────────

/// Cardinal or inter-cardinal compass direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Aspect {
    /// North (azimuth ~0°).
    North,
    /// North-east (~45°).
    NorthEast,
    /// East (~90°).
    East,
    /// South-east (~135°).
    SouthEast,
    /// South (~180°).
    South,
    /// South-west (~225°).
    SouthWest,
    /// West (~270°).
    West,
    /// North-west (~315°).
    NorthWest,
    /// Flat (slope below threshold; no aspect).
    Flat,
}

/// Per-cell slope and aspect data.
#[derive(Debug, Clone, Copy)]
pub struct SlopeAspectCell {
    /// Slope angle in radians.
    pub slope_rad: f64,
    /// Aspect azimuth in radians (measured clockwise from North = +Z axis).
    pub aspect_rad: f64,
}

impl SlopeAspectCell {
    /// Cardinal aspect category (flat threshold 0.5°).
    pub fn aspect_cardinal(&self) -> Aspect {
        const FLAT_THRESH: f64 = 0.5f64 * std::f64::consts::PI / 180.0;
        if self.slope_rad < FLAT_THRESH {
            return Aspect::Flat;
        }
        let deg = self.aspect_rad.to_degrees().rem_euclid(360.0);
        match deg as u32 {
            0..=22 | 338..=360 => Aspect::North,
            23..=67 => Aspect::NorthEast,
            68..=112 => Aspect::East,
            113..=157 => Aspect::SouthEast,
            158..=202 => Aspect::South,
            203..=247 => Aspect::SouthWest,
            248..=292 => Aspect::West,
            293..=337 => Aspect::NorthWest,
            _ => Aspect::North,
        }
    }
}

/// A slope-and-aspect map computed from a [`Heightfield`].
#[derive(Debug, Clone)]
pub struct SlopeAspectMap {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// Per-cell data.
    pub cells: Vec<SlopeAspectCell>,
}

impl SlopeAspectMap {
    /// Compute the slope-aspect map for the given heightfield.
    pub fn compute(hf: &Heightfield) -> Self {
        let n = hf.nx * hf.ny;
        let mut cells = Vec::with_capacity(n);
        for j in 0..hf.ny {
            for i in 0..hf.nx {
                let ci = i.clamp(1, hf.nx.saturating_sub(2));
                let cj = j.clamp(1, hf.ny.saturating_sub(2));
                let h = |x: usize, z: usize| hf.heights[z * hf.nx + x];
                let dhdx = (h(ci + 1, cj) - h(ci - 1, cj)) / (2.0 * hf.dx);
                let dhdz = (h(ci, cj + 1) - h(ci, cj - 1)) / (2.0 * hf.dz);
                let slope_rad = (dhdx * dhdx + dhdz * dhdz).sqrt().atan();
                // Aspect: azimuth of the downslope direction (clockwise from +Z = North)
                let aspect_rad = dhdz.atan2(dhdx);
                cells.push(SlopeAspectCell {
                    slope_rad,
                    aspect_rad,
                });
            }
        }
        Self {
            nx: hf.nx,
            ny: hf.ny,
            cells,
        }
    }

    /// Get the cell at `(i, j)`.
    pub fn cell(&self, i: usize, j: usize) -> Option<&SlopeAspectCell> {
        if i >= self.nx || j >= self.ny {
            return None;
        }
        Some(&self.cells[j * self.nx + i])
    }

    /// RGBA colour for slope-aspect visualization (HSV mapped).
    ///
    /// Hue encodes aspect; saturation and value encode slope.
    pub fn color_for_cell(&self, i: usize, j: usize) -> [f32; 4] {
        let cell = match self.cell(i, j) {
            Some(c) => c,
            None => return [0.0, 0.0, 0.0, 0.0],
        };
        // Map aspect azimuth [0, 2π] → hue [0, 1]
        let hue = (cell.aspect_rad + std::f64::consts::PI) / (2.0 * std::f64::consts::PI);
        // Slope → saturation (more slope = more saturated)
        let sat = (cell.slope_rad / (std::f64::consts::PI * 0.5)).clamp(0.0, 1.0);
        let val = 1.0f32;
        // Simple HSV → RGB
        let h6 = hue as f32 * 6.0;
        let hi = h6 as u32 % 6;
        let f = h6 - h6.floor();
        let s = sat as f32;
        let p = val * (1.0 - s);
        let q = val * (1.0 - s * f);
        let t_ = val * (1.0 - s * (1.0 - f));
        let (r, g, b) = match hi {
            0 => (val, t_, p),
            1 => (q, val, p),
            2 => (p, val, t_),
            3 => (p, q, val),
            4 => (t_, p, val),
            _ => (val, p, q),
        };
        [r, g, b, 1.0]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TerrainShadow — ray-marching shadow map
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for a directional sun light for shadow computation.
#[derive(Debug, Clone, Copy)]
pub struct SunLight {
    /// Azimuth angle in radians (clockwise from +Z = North).
    pub azimuth_rad: f64,
    /// Elevation angle above the horizon in radians.
    pub elevation_rad: f64,
}

impl SunLight {
    /// Create a new sun light.
    pub fn new(azimuth_rad: f64, elevation_rad: f64) -> Self {
        Self {
            azimuth_rad,
            elevation_rad,
        }
    }

    /// Unit direction vector toward the sun `[dx, dy, dz]`.
    pub fn direction(&self) -> [f64; 3] {
        let el = self.elevation_rad;
        let az = self.azimuth_rad;
        let cos_el = el.cos();
        [cos_el * az.sin(), el.sin(), cos_el * az.cos()]
    }
}

/// Per-cell shadow flags for a terrain.
#[derive(Debug, Clone)]
pub struct TerrainShadow {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// `true` where the cell is in shadow.
    pub shadow: Vec<bool>,
}

impl TerrainShadow {
    /// Compute shadows by ray-marching from each cell toward the sun.
    ///
    /// For each surface point, a ray is cast in the sun direction and the
    /// height is compared to the terrain along the ray.  If any terrain sample
    /// exceeds the ray height, the origin is in shadow.
    ///
    /// `max_steps` limits the march length (longer = more accurate shadows).
    pub fn compute(hf: &Heightfield, sun: SunLight, max_steps: usize) -> Self {
        let [sdx, sdy, sdz] = sun.direction();
        let mut shadow = vec![false; hf.nx * hf.ny];

        for j in 0..hf.ny {
            for i in 0..hf.nx {
                let origin_y = hf.heights[j * hf.nx + i];
                let ox = i as f64 * hf.dx;
                let oz = j as f64 * hf.dz;

                let step_len = hf.dx.min(hf.dz);
                let in_shadow = (1..=max_steps).any(|s| {
                    let t = s as f64 * step_len;
                    let rx = ox + t * sdx;
                    let rz = oz + t * sdz;
                    let ry = origin_y + t * sdy;

                    // Grid indices for the ray position
                    let ri = (rx / hf.dx) as isize;
                    let rj = (rz / hf.dz) as isize;
                    if ri < 0 || rj < 0 {
                        return false;
                    }
                    let ri = ri as usize;
                    let rj = rj as usize;
                    if ri >= hf.nx || rj >= hf.ny {
                        return false;
                    }
                    hf.heights[rj * hf.nx + ri] > ry
                });
                shadow[j * hf.nx + i] = in_shadow;
            }
        }
        Self {
            nx: hf.nx,
            ny: hf.ny,
            shadow,
        }
    }

    /// Return `true` if cell `(i, j)` is in shadow.
    pub fn is_shadowed(&self, i: usize, j: usize) -> bool {
        if i >= self.nx || j >= self.ny {
            return false;
        }
        self.shadow[j * self.nx + i]
    }

    /// Fraction of the terrain in shadow.
    pub fn shadow_fraction(&self) -> f64 {
        let n = self.shadow.len();
        if n == 0 {
            return 0.0;
        }
        let count = self.shadow.iter().filter(|&&s| s).count();
        count as f64 / n as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CrossSection — terrain profile along a polyline
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D sample along a terrain cross-section profile.
#[derive(Debug, Clone, Copy)]
pub struct ProfileSample {
    /// Cumulative horizontal distance from the start of the profile in world units.
    pub distance: f64,
    /// Terrain elevation (Y) at this sample.
    pub elevation: f64,
}

/// A terrain cross-section profile sampled along a polyline in the XZ plane.
#[derive(Debug, Clone)]
pub struct CrossSection {
    /// Profile samples ordered by increasing distance.
    pub samples: Vec<ProfileSample>,
}

impl CrossSection {
    /// Sample the terrain along the line segment from `(x0, z0)` to `(x1, z1)`.
    ///
    /// `n_samples` controls resolution (inclusive at both ends).
    pub fn along_segment(
        hf: &Heightfield,
        x0: f64,
        z0: f64,
        x1: f64,
        z1: f64,
        n_samples: usize,
    ) -> Self {
        let n = n_samples.max(2);
        let total_dist = ((x1 - x0).powi(2) + (z1 - z0).powi(2)).sqrt();
        let mut samples = Vec::with_capacity(n);
        for k in 0..n {
            let t = k as f64 / (n - 1) as f64;
            let wx = x0 + t * (x1 - x0);
            let wz = z0 + t * (z1 - z0);
            let u = wx / hf.dx;
            let v = wz / hf.dz;
            let elevation = hf.sample(u, v).unwrap_or(0.0);
            samples.push(ProfileSample {
                distance: t * total_dist,
                elevation,
            });
        }
        Self { samples }
    }

    /// Maximum elevation along the profile.
    pub fn max_elevation(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.elevation)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum elevation along the profile.
    pub fn min_elevation(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.elevation)
            .fold(f64::INFINITY, f64::min)
    }

    /// Total horizontal distance of the profile.
    pub fn length(&self) -> f64 {
        self.samples.last().map_or(0.0, |s| s.distance)
    }

    /// Elevation gain (sum of all positive increments).
    pub fn elevation_gain(&self) -> f64 {
        self.samples
            .windows(2)
            .map(|w| (w[1].elevation - w[0].elevation).max(0.0))
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElevationColorRamp — elevation-to-color mapping
// ─────────────────────────────────────────────────────────────────────────────

/// A single stop in an elevation color ramp.
#[derive(Debug, Clone)]
pub struct ColorStop {
    /// Elevation at this stop.
    pub elevation: f64,
    /// RGBA color `[r, g, b, a]` in `[0, 1]`.
    pub color: [f32; 4],
}

/// An elevation-to-color ramp with piecewise-linear interpolation between stops.
///
/// Stops must be sorted by ascending elevation before calling `sample`.
#[derive(Debug, Clone, Default)]
pub struct ElevationColorRamp {
    /// Color stops, should be sorted by `elevation`.
    pub stops: Vec<ColorStop>,
}

impl ElevationColorRamp {
    /// Create an empty ramp.
    pub fn new() -> Self {
        Self { stops: Vec::new() }
    }

    /// Add a color stop and maintain sorted order.
    pub fn add_stop(&mut self, elevation: f64, color: [f32; 4]) {
        self.stops.push(ColorStop { elevation, color });
        self.stops.sort_by(|a, b| {
            a.elevation
                .partial_cmp(&b.elevation)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Sample the ramp at a given elevation.
    ///
    /// Returns the color of the nearest stop if only one stop exists.
    /// Returns `[0.5, 0.5, 0.5, 1.0]` for an empty ramp.
    pub fn sample(&self, elevation: f64) -> [f32; 4] {
        match self.stops.len() {
            0 => [0.5, 0.5, 0.5, 1.0],
            1 => self.stops[0].color,
            _ => {
                if elevation <= self.stops[0].elevation {
                    return self.stops[0].color;
                }
                let last = self.stops.last().expect("collection should not be empty");
                if elevation >= last.elevation {
                    return last.color;
                }
                let idx = self.stops.partition_point(|s| s.elevation <= elevation);
                let lo = &self.stops[idx - 1];
                let hi = &self.stops[idx];
                let t = ((elevation - lo.elevation) / (hi.elevation - lo.elevation)) as f32;
                let mut out = [0.0f32; 4];
                for k in 0..4 {
                    out[k] = lo.color[k] + t * (hi.color[k] - lo.color[k]);
                }
                out
            }
        }
    }

    /// Build a standard hypsometric (elevation-based) ramp.
    ///
    /// Ranges from deep ocean blue at the lowest elevation to snow-white at the
    /// highest, passing through green and brown in between.
    pub fn hypsometric(min_elev: f64, max_elev: f64) -> Self {
        let mut ramp = Self::new();
        let range = max_elev - min_elev;
        if range <= 0.0 {
            return ramp;
        }
        ramp.add_stop(min_elev, [0.05, 0.22, 0.56, 1.0]); // deep blue
        ramp.add_stop(min_elev + 0.2 * range, [0.20, 0.60, 0.25, 1.0]); // lowland green
        ramp.add_stop(min_elev + 0.5 * range, [0.55, 0.45, 0.25, 1.0]); // mid brown
        ramp.add_stop(min_elev + 0.75 * range, [0.70, 0.60, 0.45, 1.0]); // high brown
        ramp.add_stop(max_elev, [0.97, 0.97, 0.97, 1.0]); // snow white
        ramp
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Heightfield generation helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Generate a heightfield from a sinusoidal test function.
///
/// `h(i,j) = amp * sin(freq_x * i * dx) * cos(freq_z * j * dz)`
pub fn sinusoidal_terrain(
    nx: usize,
    ny: usize,
    dx: f64,
    dz: f64,
    amp: f64,
    freq_x: f64,
    freq_z: f64,
) -> Heightfield {
    let mut hf = Heightfield::new(nx, ny, dx, dz);
    for j in 0..ny {
        for i in 0..nx {
            let h = amp * (freq_x * i as f64 * dx).sin() * (freq_z * j as f64 * dz).cos();
            hf.set_height(i, j, h);
        }
    }
    hf
}

/// Generate a simple conical hill centred at the middle of the heightfield.
///
/// Height falls linearly from `peak_height` at the centre to zero at the edges.
pub fn conical_terrain(nx: usize, ny: usize, dx: f64, dz: f64, peak_height: f64) -> Heightfield {
    let mut hf = Heightfield::new(nx, ny, dx, dz);
    let cx = (nx as f64 - 1.0) * 0.5 * dx;
    let cz = (ny as f64 - 1.0) * 0.5 * dz;
    let max_r = (cx * cx + cz * cz).sqrt();
    for j in 0..ny {
        for i in 0..nx {
            let px = i as f64 * dx - cx;
            let pz = j as f64 * dz - cz;
            let r = (px * px + pz * pz).sqrt();
            let h = peak_height * (1.0 - (r / max_r).min(1.0));
            hf.set_height(i, j, h);
        }
    }
    hf
}

// ─────────────────────────────────────────────────────────────────────────────
// Hillshading — analytical hillshading from normal and sun
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the hillshade intensity `[0, 1]` at cell `(i, j)` for the given sun.
///
/// Returns the Lambertian (N·L) factor clamped to `[0, 1]`.
pub fn hillshade(hf: &Heightfield, i: usize, j: usize, sun: SunLight) -> f64 {
    let [nx, ny, nz] = hf.normal(i, j);
    let [lx, ly, lz] = sun.direction();
    let dot = nx * lx + ny * ly + nz * lz;
    dot.clamp(0.0, 1.0)
}

/// Compute a full hillshade image (one intensity per cell).
pub fn hillshade_map(hf: &Heightfield, sun: SunLight) -> Vec<f64> {
    let mut map = Vec::with_capacity(hf.nx * hf.ny);
    for j in 0..hf.ny {
        for i in 0..hf.nx {
            map.push(hillshade(hf, i, j, sun));
        }
    }
    map
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Heightfield ──────────────────────────────────────────────────────────

    #[test]
    fn test_heightfield_basic_set_get() {
        let mut hf = Heightfield::new(4, 4, 1.0, 1.0);
        assert!(hf.set_height(2, 3, 7.5));
        assert!((hf.height(2, 3).unwrap() - 7.5).abs() < 1e-12);
    }

    #[test]
    fn test_heightfield_out_of_bounds() {
        let hf = Heightfield::new(4, 4, 1.0, 1.0);
        assert!(hf.height(4, 0).is_none());
        assert!(hf.height(0, 4).is_none());
    }

    #[test]
    fn test_heightfield_bilinear_sample_corners() {
        let mut hf = Heightfield::new(3, 3, 1.0, 1.0);
        hf.set_height(0, 0, 0.0);
        hf.set_height(1, 0, 1.0);
        hf.set_height(0, 1, 0.0);
        hf.set_height(1, 1, 1.0);
        // At u=0.5, v=0.0 → bilinear gives 0.5
        let h = hf.sample(0.5, 0.0).unwrap();
        assert!((h - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_heightfield_bilinear_sample_out_of_range() {
        let hf = Heightfield::new(3, 3, 1.0, 1.0);
        assert!(hf.sample(-1.0, 0.0).is_none());
        assert!(hf.sample(5.0, 0.0).is_none());
    }

    #[test]
    fn test_heightfield_world_pos() {
        let mut hf = Heightfield::new(5, 5, 2.0, 3.0);
        hf.set_height(2, 1, 4.0);
        let p = hf.world_pos(2, 1).unwrap();
        assert!((p[0] - 4.0).abs() < 1e-9); // x = 2*2
        assert!((p[1] - 4.0).abs() < 1e-9); // y = height
        assert!((p[2] - 3.0).abs() < 1e-9); // z = 1*3
    }

    #[test]
    fn test_heightfield_min_max() {
        let mut hf = Heightfield::new(3, 3, 1.0, 1.0);
        hf.heights = vec![-5.0, 0.0, 10.0, 3.0, 7.0, 1.0, 2.0, 8.0, -1.0];
        assert!((hf.min_height() - (-5.0)).abs() < 1e-9);
        assert!((hf.max_height() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_heightfield_lod_mesh_step1() {
        let hf = Heightfield::new(4, 4, 1.0, 1.0);
        let (verts, tris) = hf.lod_mesh(1);
        assert_eq!(verts.len(), 16); // 4×4
        assert_eq!(tris.len(), 18); // 3×3 quads × 2 triangles
    }

    #[test]
    fn test_heightfield_lod_mesh_step2() {
        let hf = Heightfield::new(5, 5, 1.0, 1.0);
        let (verts_full, _) = hf.lod_mesh(1);
        let (verts_lod, _) = hf.lod_mesh(2);
        assert!(verts_lod.len() < verts_full.len());
    }

    #[test]
    fn test_heightfield_normal_flat() {
        let hf = Heightfield::new(5, 5, 1.0, 1.0); // all zeros → flat
        let n = hf.normal(2, 2);
        // Flat terrain → normal should be [0, 1, 0]
        assert!(n[0].abs() < 1e-6);
        assert!((n[1] - 1.0).abs() < 1e-6);
        assert!(n[2].abs() < 1e-6);
    }

    #[test]
    fn test_heightfield_slope_flat() {
        let hf = Heightfield::new(5, 5, 1.0, 1.0);
        let s = hf.slope_rad(2, 2);
        assert!(s.abs() < 1e-9);
    }

    // ── SplatLayer ────────────────────────────────────────────────────────────

    #[test]
    fn test_splat_layer_weight_in_range() {
        let layer = SplatLayer::new("grass", 0.0, 0.5, 0.0, 100.0, [0.0, 1.0, 0.0, 1.0]);
        let w = layer.weight(0.25, 50.0);
        assert!(w > 0.0 && w <= 1.0);
    }

    #[test]
    fn test_splat_layer_weight_out_of_altitude() {
        let layer = SplatLayer::new("rock", 0.0, 1.0, 200.0, 400.0, [0.5, 0.5, 0.5, 1.0]);
        let w = layer.weight(0.5, 50.0); // altitude 50 < alt_min 200
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_texture_splat_blend_normalised() {
        let mut splat = TextureSplat::new();
        splat.add_layer(SplatLayer::new(
            "a",
            0.0,
            1.0,
            0.0,
            100.0,
            [1.0, 0.0, 0.0, 1.0],
        ));
        splat.add_layer(SplatLayer::new(
            "b",
            0.0,
            1.0,
            0.0,
            100.0,
            [0.0, 1.0, 0.0, 1.0],
        ));
        let c = splat.blend(0.3, 50.0);
        // Alpha channel should be normalized to ≈1
        assert!(
            (c[3] - 1.0).abs() < 0.05,
            "alpha should be ~1, got {}",
            c[3]
        );
    }

    #[test]
    fn test_texture_splat_empty() {
        let splat = TextureSplat::new();
        let c = splat.blend(0.0, 0.0);
        assert_eq!(c, [0.5, 0.5, 0.5, 1.0]);
    }

    // ── ContourExtractor ──────────────────────────────────────────────────────

    #[test]
    fn test_contour_flat_terrain_no_segments() {
        let hf = Heightfield::new(4, 4, 1.0, 1.0); // all zeros
        let ext = ContourExtractor::new(&hf);
        // Contour at 1.0: no terrain reaches it
        let segs = ext.extract_contour(1.0);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_contour_single_peak() {
        let mut hf = Heightfield::new(5, 5, 1.0, 1.0);
        hf.set_height(2, 2, 2.0); // single high point
        let ext = ContourExtractor::new(&hf);
        let segs = ext.extract_contour(1.0);
        assert!(!segs.is_empty(), "should find contour segments around peak");
    }

    #[test]
    fn test_contour_extract_all() {
        let hf = sinusoidal_terrain(8, 8, 1.0, 1.0, 3.0, 0.5, 0.5);
        let ext = ContourExtractor::new(&hf);
        let segs = ext.extract_all(0.5);
        // Should produce some segments for a wave terrain
        assert!(!segs.is_empty());
    }

    #[test]
    fn test_contour_elevation_tagged() {
        let hf = conical_terrain(5, 5, 1.0, 1.0, 10.0);
        let ext = ContourExtractor::new(&hf);
        let segs = ext.extract_contour(3.0);
        for seg in &segs {
            assert!((seg.elevation - 3.0).abs() < 1e-9);
        }
    }

    // ── WaterSurface ──────────────────────────────────────────────────────────

    #[test]
    fn test_water_surface_flood_fraction_all() {
        let hf = Heightfield::new(4, 4, 1.0, 1.0); // all heights = 0
        let water = WaterSurface::new(&hf, 1.0, [0.0, 0.3, 0.8, 0.8]); // level above all
        assert!((water.flood_fraction() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_water_surface_flood_fraction_none() {
        let mut hf = Heightfield::new(4, 4, 1.0, 1.0);
        for v in hf.heights.iter_mut() {
            *v = 10.0;
        }
        let water = WaterSurface::new(&hf, 1.0, [0.0, 0.3, 0.8, 0.8]);
        assert!(water.flood_fraction() < 1e-9);
    }

    #[test]
    fn test_water_surface_shoreline() {
        let hf = conical_terrain(7, 7, 1.0, 1.0, 5.0);
        let water = WaterSurface::new(&hf, 2.0, [0.0, 0.0, 1.0, 1.0]);
        let shore = water.shoreline();
        assert!(
            !shore.is_empty(),
            "conical hill with water level should have shoreline"
        );
    }

    #[test]
    fn test_water_quad_has_two_triangles() {
        let hf = Heightfield::new(4, 4, 2.0, 2.0);
        let water = WaterSurface::new(&hf, 1.0, [0.0, 0.0, 1.0, 1.0]);
        let (verts, tris) = water.water_quad();
        assert_eq!(verts.len(), 4);
        assert_eq!(tris.len(), 2);
    }

    // ── ErosionMap ────────────────────────────────────────────────────────────

    #[test]
    fn test_erosion_map_set_get() {
        let mut em = ErosionMap::new(4, 4);
        assert!(em.set_cell(
            1,
            2,
            ErosionCell {
                net_sediment_m: -3.0,
                water_flux: 0.5
            }
        ));
        let c = em.cell(1, 2).unwrap();
        assert!((c.net_sediment_m - (-3.0)).abs() < 1e-9);
    }

    #[test]
    fn test_erosion_map_color_negative() {
        let mut em = ErosionMap::new(2, 2);
        em.set_cell(
            0,
            0,
            ErosionCell {
                net_sediment_m: -10.0,
                water_flux: 0.0,
            },
        );
        let c = em.color_for_cell(0, 0);
        assert!(c[0] > c[1], "erosion (negative) should be reddish");
    }

    #[test]
    fn test_erosion_map_color_positive() {
        let mut em = ErosionMap::new(2, 2);
        em.set_cell(
            0,
            0,
            ErosionCell {
                net_sediment_m: 10.0,
                water_flux: 0.0,
            },
        );
        let c = em.color_for_cell(0, 0);
        assert!(c[2] > c[0], "deposition (positive) should be bluish");
    }

    #[test]
    fn test_erosion_map_max_abs_sediment() {
        let mut em = ErosionMap::new(3, 3);
        em.set_cell(
            1,
            1,
            ErosionCell {
                net_sediment_m: -7.0,
                water_flux: 0.0,
            },
        );
        em.set_cell(
            0,
            0,
            ErosionCell {
                net_sediment_m: 3.0,
                water_flux: 0.0,
            },
        );
        assert!((em.max_abs_sediment() - 7.0).abs() < 1e-9);
    }

    // ── VegetationDensityMap ──────────────────────────────────────────────────

    #[test]
    fn test_vegetation_map_set_get() {
        let mut vm = VegetationDensityMap::new(3, 3);
        vm.set_cell(
            1,
            1,
            VegetationCell {
                density: 0.8,
                plant_type: 2,
            },
        );
        let c = vm.cell(1, 1).unwrap();
        assert!((c.density - 0.8).abs() < 1e-6);
        assert_eq!(c.plant_type, 2);
    }

    #[test]
    fn test_vegetation_map_color_green() {
        let mut vm = VegetationDensityMap::new(2, 2);
        vm.set_cell(
            0,
            0,
            VegetationCell {
                density: 1.0,
                plant_type: 0,
            },
        );
        let c = vm.color_for_cell(0, 0);
        assert!(c[1] > 0.9, "high-density cell should be green");
        assert!(c[0] < 0.1, "red channel should be near zero");
    }

    #[test]
    fn test_vegetation_map_mean_density() {
        let mut vm = VegetationDensityMap::new(2, 1);
        vm.set_cell(
            0,
            0,
            VegetationCell {
                density: 0.4,
                plant_type: 0,
            },
        );
        vm.set_cell(
            1,
            0,
            VegetationCell {
                density: 0.8,
                plant_type: 0,
            },
        );
        assert!((vm.mean_density() - 0.6).abs() < 1e-6);
    }

    #[test]
    fn test_vegetation_populate_from_terrain() {
        // Use gentle slope (peak_height=1.0) so slope_rad < max_slope_rad=0.5
        let hf = conical_terrain(7, 7, 1.0, 1.0, 1.0);
        let mut vm = VegetationDensityMap::new(7, 7);
        vm.populate_from_terrain(&hf, 0.5, 0.0, 1.0);
        // At least some cells should have non-zero density
        assert!(vm.cells.iter().any(|c| c.density > 0.0));
    }

    // ── SlopeAspectMap ────────────────────────────────────────────────────────

    #[test]
    fn test_slope_aspect_flat() {
        let hf = Heightfield::new(5, 5, 1.0, 1.0);
        let sa = SlopeAspectMap::compute(&hf);
        let cell = sa.cell(2, 2).unwrap();
        assert!(
            cell.slope_rad.abs() < 1e-9,
            "flat terrain should have zero slope"
        );
        assert_eq!(cell.aspect_cardinal(), Aspect::Flat);
    }

    #[test]
    fn test_slope_aspect_tilted() {
        let mut hf = Heightfield::new(5, 5, 1.0, 1.0);
        // Slope increasing in X direction
        for j in 0..5 {
            for i in 0..5 {
                hf.set_height(i, j, i as f64 * 2.0);
            }
        }
        let sa = SlopeAspectMap::compute(&hf);
        let cell = sa.cell(2, 2).unwrap();
        assert!(
            cell.slope_rad > 0.0,
            "tilted terrain should have positive slope"
        );
    }

    #[test]
    fn test_slope_aspect_color_valid() {
        let hf = sinusoidal_terrain(8, 8, 1.0, 1.0, 2.0, 1.0, 1.0);
        let sa = SlopeAspectMap::compute(&hf);
        let c = sa.color_for_cell(3, 3);
        for k in 0..4 {
            assert!(
                c[k] >= 0.0 && c[k] <= 1.0,
                "color component {} out of range: {}",
                k,
                c[k]
            );
        }
    }

    // ── TerrainShadow ─────────────────────────────────────────────────────────

    #[test]
    fn test_shadow_flat_no_self_shadow() {
        let hf = Heightfield::new(6, 6, 1.0, 1.0);
        let sun = SunLight::new(0.0, std::f64::consts::FRAC_PI_4); // 45° elevation
        let shad = TerrainShadow::compute(&hf, sun, 10);
        // A completely flat terrain should cast no self-shadows
        assert!(shad.shadow_fraction() < 1e-9);
    }

    #[test]
    fn test_shadow_hill_casts_shadow() {
        let hf = conical_terrain(11, 11, 1.0, 1.0, 5.0);
        let sun = SunLight::new(0.0, 0.3); // low sun angle
        let shad = TerrainShadow::compute(&hf, sun, 20);
        assert!(shad.shadow_fraction() > 0.0, "hill should cast some shadow");
    }

    #[test]
    fn test_shadow_fraction_range() {
        let hf = conical_terrain(7, 7, 1.0, 1.0, 3.0);
        let sun = SunLight::new(1.0, 0.2);
        let shad = TerrainShadow::compute(&hf, sun, 15);
        let f = shad.shadow_fraction();
        assert!((0.0..=1.0).contains(&f));
    }

    // ── CrossSection ──────────────────────────────────────────────────────────

    #[test]
    fn test_cross_section_flat_terrain() {
        let hf = Heightfield::new(10, 10, 1.0, 1.0);
        let cs = CrossSection::along_segment(&hf, 0.0, 0.0, 5.0, 5.0, 6);
        assert_eq!(cs.samples.len(), 6);
        for s in &cs.samples {
            assert!(
                s.elevation.abs() < 1e-9,
                "flat terrain should have zero elevation"
            );
        }
    }

    #[test]
    fn test_cross_section_length() {
        let hf = Heightfield::new(10, 10, 1.0, 1.0);
        let cs = CrossSection::along_segment(&hf, 0.0, 0.0, 3.0, 4.0, 10);
        // 3-4-5 right triangle
        assert!((cs.length() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_cross_section_elevation_gain() {
        let hf = conical_terrain(11, 11, 1.0, 1.0, 10.0);
        let cs = CrossSection::along_segment(&hf, 0.0, 5.0, 10.0, 5.0, 20);
        let gain = cs.elevation_gain();
        assert!(
            gain > 0.0,
            "cross section over a hill should have positive gain"
        );
    }

    // ── ElevationColorRamp ────────────────────────────────────────────────────

    #[test]
    fn test_elevation_ramp_empty() {
        let ramp = ElevationColorRamp::new();
        let c = ramp.sample(0.0);
        assert_eq!(c, [0.5, 0.5, 0.5, 1.0]);
    }

    #[test]
    fn test_elevation_ramp_single_stop() {
        let mut ramp = ElevationColorRamp::new();
        ramp.add_stop(5.0, [1.0, 0.0, 0.0, 1.0]);
        let c = ramp.sample(3.0);
        assert_eq!(c, [1.0, 0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_elevation_ramp_clamp_below() {
        let mut ramp = ElevationColorRamp::new();
        ramp.add_stop(0.0, [0.0, 0.0, 1.0, 1.0]);
        ramp.add_stop(10.0, [1.0, 0.0, 0.0, 1.0]);
        let c = ramp.sample(-5.0);
        assert_eq!(c, [0.0, 0.0, 1.0, 1.0]);
    }

    #[test]
    fn test_elevation_ramp_clamp_above() {
        let mut ramp = ElevationColorRamp::new();
        ramp.add_stop(0.0, [0.0, 0.0, 1.0, 1.0]);
        ramp.add_stop(10.0, [1.0, 0.0, 0.0, 1.0]);
        let c = ramp.sample(20.0);
        assert_eq!(c, [1.0, 0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_elevation_ramp_midpoint_interpolation() {
        let mut ramp = ElevationColorRamp::new();
        ramp.add_stop(0.0, [0.0, 0.0, 0.0, 1.0]);
        ramp.add_stop(10.0, [1.0, 1.0, 1.0, 1.0]);
        let c = ramp.sample(5.0);
        assert!((c[0] - 0.5).abs() < 1e-5, "midpoint should be grey");
    }

    #[test]
    fn test_elevation_ramp_hypsometric_endpoints() {
        let ramp = ElevationColorRamp::hypsometric(0.0, 100.0);
        let low = ramp.sample(0.0);
        let high = ramp.sample(100.0);
        // Low → bluish, high → whitish
        assert!(low[2] > low[0], "sea level should be bluish (blue > red)");
        assert!(high[0] > 0.8, "peak should be whitish (high red)");
    }

    #[test]
    fn test_elevation_ramp_sorted_after_add() {
        let mut ramp = ElevationColorRamp::new();
        ramp.add_stop(10.0, [1.0, 0.0, 0.0, 1.0]);
        ramp.add_stop(0.0, [0.0, 1.0, 0.0, 1.0]);
        ramp.add_stop(5.0, [0.0, 0.0, 1.0, 1.0]);
        for w in ramp.stops.windows(2) {
            assert!(w[0].elevation <= w[1].elevation, "stops not sorted");
        }
    }

    // ── Terrain generation helpers ────────────────────────────────────────────

    #[test]
    fn test_sinusoidal_terrain_range() {
        let hf = sinusoidal_terrain(10, 10, 0.5, 0.5, 2.0, 1.0, 1.0);
        let (mn, mx) = hf.height_range();
        assert!(mn >= -2.0 - 1e-9);
        assert!(mx <= 2.0 + 1e-9);
    }

    #[test]
    fn test_conical_terrain_peak_at_center() {
        let hf = conical_terrain(11, 11, 1.0, 1.0, 5.0);
        let center_h = hf.height(5, 5).unwrap();
        assert!((center_h - 5.0).abs() < 1e-9, "peak should be at centre");
        let edge_h = hf.height(0, 0).unwrap();
        assert!(edge_h < center_h, "edge should be lower than peak");
    }

    // ── Hillshading ───────────────────────────────────────────────────────────

    #[test]
    fn test_hillshade_flat_overhead_sun() {
        let hf = Heightfield::new(5, 5, 1.0, 1.0);
        // Sun directly overhead (elevation = π/2)
        let sun = SunLight::new(0.0, std::f64::consts::FRAC_PI_2);
        let hs = hillshade(&hf, 2, 2, sun);
        assert!(
            (hs - 1.0).abs() < 1e-6,
            "overhead sun on flat terrain → full intensity"
        );
    }

    #[test]
    fn test_hillshade_map_length() {
        let hf = sinusoidal_terrain(6, 8, 1.0, 1.0, 1.0, 1.0, 1.0);
        let sun = SunLight::new(0.5, 0.5);
        let map = hillshade_map(&hf, sun);
        assert_eq!(map.len(), 6 * 8);
    }

    #[test]
    fn test_hillshade_values_in_range() {
        let hf = conical_terrain(9, 9, 1.0, 1.0, 3.0);
        let sun = SunLight::new(0.7, 0.6);
        let map = hillshade_map(&hf, sun);
        for &v in &map {
            assert!(
                (0.0..=1.0).contains(&v),
                "hillshade value {} out of [0,1]",
                v
            );
        }
    }

    // ── SunLight ──────────────────────────────────────────────────────────────

    #[test]
    fn test_sun_direction_unit_length() {
        let sun = SunLight::new(1.2, 0.4);
        let [dx, dy, dz] = sun.direction();
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-9,
            "sun direction should be unit vector"
        );
    }

    #[test]
    fn test_sun_direction_overhead() {
        let sun = SunLight::new(0.0, std::f64::consts::FRAC_PI_2);
        let [_dx, dy, _dz] = sun.direction();
        assert!((dy - 1.0).abs() < 1e-6, "overhead sun → dy = 1");
    }
}
