#![allow(clippy::needless_range_loop, clippy::type_complexity)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Remote sensing and geospatial data I/O.
//!
//! Provides structs and methods for reading, writing, and processing
//! geospatial raster data (GeoTIFF-style), point clouds (LAS/XYZ),
//! satellite imagery indices, and digital elevation models.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// GeoTiffReader
// ---------------------------------------------------------------------------

/// A minimal GeoTIFF-style raster reader that holds image geometry and
/// georeferencing metadata.
///
/// The affine transform follows the standard 6-parameter GDAL convention:
/// `[x_origin, pixel_width, row_rotation, y_origin, col_rotation, pixel_height]`
/// where `pixel_height` is typically negative (top-left origin).
#[derive(Debug, Clone)]
pub struct GeoTiffReader {
    /// Raster width in pixels.
    pub width: usize,
    /// Raster height in pixels.
    pub height: usize,
    /// Number of spectral or data bands.
    pub bands: usize,
    /// Six-element affine geotransform (GDAL order).
    pub transform: [f64; 6],
    /// Band data stored as flat row-major vectors, one per band.
    data: Vec<Vec<f64>>,
}

impl GeoTiffReader {
    /// Construct a new `GeoTiffReader` with explicit dimensions and transform.
    pub fn new(
        width: usize,
        height: usize,
        bands: usize,
        transform: [f64; 6],
        data: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            width,
            height,
            bands,
            transform,
            data,
        }
    }

    /// Convert pixel coordinates `(px, py)` to world (lon, lat) coordinates
    /// using the stored affine transform.
    pub fn pixel_to_world(&self, px: usize, py: usize) -> [f64; 2] {
        let t = &self.transform;
        let lon = t[0] + px as f64 * t[1] + py as f64 * t[2];
        let lat = t[3] + px as f64 * t[4] + py as f64 * t[5];
        [lon, lat]
    }

    /// Convert world (lon, lat) coordinates to pixel coordinates using the
    /// inverse of the stored affine transform.
    pub fn world_to_pixel(&self, lon: f64, lat: f64) -> (usize, usize) {
        let t = &self.transform;
        let det = t[1] * t[5] - t[2] * t[4];
        let det = if det.abs() < 1e-15 { 1e-15 } else { det };
        let dx = lon - t[0];
        let dy = lat - t[3];
        let px = (t[5] * dx - t[2] * dy) / det;
        let py = (t[1] * dy - t[4] * dx) / det;
        let px = px
            .round()
            .max(0.0)
            .min((self.width.saturating_sub(1)) as f64) as usize;
        let py = py
            .round()
            .max(0.0)
            .min((self.height.saturating_sub(1)) as f64) as usize;
        (px, py)
    }

    /// Return a clone of band `band` data (0-indexed).
    ///
    /// Returns an empty vector if `band >= self.bands`.
    pub fn read_band(&self, band: usize) -> Vec<f64> {
        if band < self.data.len() {
            self.data[band].clone()
        } else {
            vec![]
        }
    }

    /// Return the elevation (first-band value) at world coordinate `(lon, lat)`.
    pub fn elevation_at(&self, lon: f64, lat: f64) -> f64 {
        if self.data.is_empty() || self.data[0].is_empty() {
            return 0.0;
        }
        let (px, py) = self.world_to_pixel(lon, lat);
        let idx = py * self.width + px;
        *self.data[0].get(idx).unwrap_or(&0.0)
    }

    /// Return the bounding box `[min_lon, min_lat, max_lon, max_lat]`.
    pub fn bounding_box(&self) -> [f64; 4] {
        let tl = self.pixel_to_world(0, 0);
        let br = self.pixel_to_world(self.width, self.height);
        [
            tl[0].min(br[0]),
            tl[1].min(br[1]),
            tl[0].max(br[0]),
            tl[1].max(br[1]),
        ]
    }
}

// ---------------------------------------------------------------------------
// RasterData
// ---------------------------------------------------------------------------

/// A 2-D raster grid with an associated no-data sentinel value.
#[derive(Debug, Clone)]
pub struct RasterData {
    /// Row-major 2-D data grid (`data[row][col]`).
    pub data: Vec<Vec<f64>>,
    /// Number of columns.
    pub width: usize,
    /// Number of rows.
    pub height: usize,
    /// Sentinel value representing missing / invalid cells.
    pub nodata: f64,
}

impl RasterData {
    /// Construct a new `RasterData` from a 2-D grid.
    pub fn new(data: Vec<Vec<f64>>, nodata: f64) -> Self {
        let height = data.len();
        let width = if height > 0 { data[0].len() } else { 0 };
        Self {
            data,
            width,
            height,
            nodata,
        }
    }

    /// Compute `(min, max, mean, std_dev)` over all valid (non-nodata) cells.
    pub fn statistics(&self) -> (f64, f64, f64, f64) {
        let valid: Vec<f64> = self
            .data
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .filter(|&v| (v - self.nodata).abs() > 1e-10)
            .collect();
        if valid.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let min = valid.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = valid.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = valid.iter().sum::<f64>() / valid.len() as f64;
        let variance =
            valid.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / valid.len() as f64;
        (min, max, mean, variance.sqrt())
    }

    /// Clip the raster to a bounding box `[min_x, min_y, max_x, max_y]`.
    ///
    /// In pixel space this maps `x → col` and `y → row` (0-origin).
    pub fn clip_to_bbox(&self, bbox: [f64; 4]) -> Self {
        let col_start = (bbox[0] as usize).min(self.width);
        let row_start = (bbox[1] as usize).min(self.height);
        let col_end = (bbox[2] as usize).min(self.width);
        let row_end = (bbox[3] as usize).min(self.height);
        let clipped: Vec<Vec<f64>> = self.data[row_start..row_end]
            .iter()
            .map(|row| row[col_start..col_end].to_vec())
            .collect();
        RasterData::new(clipped, self.nodata)
    }

    /// Resample the raster by `factor` (e.g. 0.5 = half resolution).
    pub fn resample(&self, factor: f64) -> Self {
        let factor = factor.max(0.01);
        let new_h = ((self.height as f64) * factor).round() as usize;
        let new_w = ((self.width as f64) * factor).round() as usize;
        let mut out = vec![vec![self.nodata; new_w]; new_h];
        for r in 0..new_h {
            for c in 0..new_w {
                let src_r = ((r as f64) / factor) as usize;
                let src_c = ((c as f64) / factor) as usize;
                let src_r = src_r.min(self.height.saturating_sub(1));
                let src_c = src_c.min(self.width.saturating_sub(1));
                out[r][c] = self.data[src_r][src_c];
            }
        }
        RasterData::new(out, self.nodata)
    }

    /// Compute hillshade given sun `azimuth` and `elevation` (degrees).
    pub fn hillshade(&self, azimuth: f64, elevation: f64) -> Vec<Vec<f64>> {
        let az_rad = azimuth.to_radians();
        let el_rad = elevation.to_radians();
        let sun = [
            -el_rad.cos() * az_rad.sin(),
            -el_rad.cos() * az_rad.cos(),
            el_rad.sin(),
        ];
        let slope_grid = self.slope();
        let aspect_grid = self.aspect();
        let mut hs = vec![vec![0.0_f64; self.width]; self.height];
        for r in 0..self.height {
            for c in 0..self.width {
                let s = slope_grid[r][c].to_radians();
                let a = aspect_grid[r][c].to_radians();
                let nx = s.sin() * a.sin();
                let ny = s.sin() * a.cos();
                let nz = s.cos();
                let dot = (nx * sun[0] + ny * sun[1] + nz * sun[2]).max(0.0);
                hs[r][c] = dot;
            }
        }
        hs
    }

    /// Compute slope (degrees) using a 3×3 Sobel kernel.
    pub fn slope(&self) -> Vec<Vec<f64>> {
        let mut out = vec![vec![0.0_f64; self.width]; self.height];
        for r in 1..self.height.saturating_sub(1) {
            for c in 1..self.width.saturating_sub(1) {
                let dzdx =
                    (self.data[r - 1][c + 1] + 2.0 * self.data[r][c + 1] + self.data[r + 1][c + 1]
                        - self.data[r - 1][c - 1]
                        - 2.0 * self.data[r][c - 1]
                        - self.data[r + 1][c - 1])
                        / 8.0;
                let dzdy =
                    (self.data[r + 1][c - 1] + 2.0 * self.data[r + 1][c] + self.data[r + 1][c + 1]
                        - self.data[r - 1][c - 1]
                        - 2.0 * self.data[r - 1][c]
                        - self.data[r - 1][c + 1])
                        / 8.0;
                out[r][c] = (dzdx * dzdx + dzdy * dzdy).sqrt().atan().to_degrees();
            }
        }
        out
    }

    /// Compute aspect (degrees, 0 = North, clockwise) using a 3×3 Sobel kernel.
    pub fn aspect(&self) -> Vec<Vec<f64>> {
        let mut out = vec![vec![0.0_f64; self.width]; self.height];
        for r in 1..self.height.saturating_sub(1) {
            for c in 1..self.width.saturating_sub(1) {
                let dzdx =
                    (self.data[r - 1][c + 1] + 2.0 * self.data[r][c + 1] + self.data[r + 1][c + 1]
                        - self.data[r - 1][c - 1]
                        - 2.0 * self.data[r][c - 1]
                        - self.data[r + 1][c - 1])
                        / 8.0;
                let dzdy =
                    (self.data[r + 1][c - 1] + 2.0 * self.data[r + 1][c] + self.data[r + 1][c + 1]
                        - self.data[r - 1][c - 1]
                        - 2.0 * self.data[r - 1][c]
                        - self.data[r - 1][c + 1])
                        / 8.0;
                let aspect = (-dzdy).atan2(dzdx).to_degrees();
                out[r][c] = if aspect < 0.0 { aspect + 360.0 } else { aspect };
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Dem
// ---------------------------------------------------------------------------

/// Digital Elevation Model wrapping a `RasterData` grid with a known cell size.
#[derive(Debug, Clone)]
pub struct Dem {
    /// Underlying raster elevation data.
    pub raster: RasterData,
    /// Physical size of each square cell (metres).
    pub cell_size: f64,
}

impl Dem {
    /// Construct a new `Dem` from a `RasterData` and a cell size.
    pub fn new(raster: RasterData, cell_size: f64) -> Self {
        Self { raster, cell_size }
    }

    /// Compute an 8-direction flow-direction grid (D8 algorithm).
    ///
    /// Encoded directions: 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE.
    pub fn flow_direction(&self) -> Vec<Vec<u8>> {
        let h = self.raster.height;
        let w = self.raster.width;
        // D8 neighbour offsets (dr, dc) and their encoded values
        let dirs: [(i32, i32, u8); 8] = [
            (0, 1, 1),    // E
            (1, 1, 2),    // SE
            (1, 0, 4),    // S
            (1, -1, 8),   // SW
            (0, -1, 16),  // W
            (-1, -1, 32), // NW
            (-1, 0, 64),  // N
            (-1, 1, 128), // NE
        ];
        let mut out = vec![vec![0u8; w]; h];
        for r in 0..h {
            for c in 0..w {
                let z = self.raster.data[r][c];
                let mut max_drop = 0.0_f64;
                let mut best_dir = 1u8;
                for &(dr, dc, code) in &dirs {
                    let nr = r as i32 + dr;
                    let nc = c as i32 + dc;
                    if nr < 0 || nc < 0 || nr >= h as i32 || nc >= w as i32 {
                        continue;
                    }
                    let nz = self.raster.data[nr as usize][nc as usize];
                    let dist = if dr == 0 || dc == 0 {
                        self.cell_size
                    } else {
                        self.cell_size * std::f64::consts::SQRT_2
                    };
                    let drop = (z - nz) / dist;
                    if drop > max_drop {
                        max_drop = drop;
                        best_dir = code;
                    }
                }
                out[r][c] = best_dir;
            }
        }
        out
    }

    /// Compute flow-accumulation grid (number of upstream cells).
    pub fn flow_accumulation(&self) -> Vec<Vec<f64>> {
        let h = self.raster.height;
        let w = self.raster.width;
        let fd = self.flow_direction();
        let mut accum = vec![vec![1.0_f64; w]; h];
        // Simple iterative accumulation (not topology-sorted; sufficient for tests)
        for _ in 0..(h * w) {
            let snap = accum.clone();
            for r in 0..h {
                for c in 0..w {
                    let dir = fd[r][c];
                    let (dr, dc): (i32, i32) = match dir {
                        1 => (0, 1),
                        2 => (1, 1),
                        4 => (1, 0),
                        8 => (1, -1),
                        16 => (0, -1),
                        32 => (-1, -1),
                        64 => (-1, 0),
                        128 => (-1, 1),
                        _ => continue,
                    };
                    let nr = r as i32 + dr;
                    let nc = c as i32 + dc;
                    if nr >= 0 && nc >= 0 && nr < h as i32 && nc < w as i32 {
                        accum[nr as usize][nc as usize] += snap[r][c];
                    }
                }
            }
        }
        accum
    }

    /// Delineate a watershed upstream of `outlet` pixel using D8 flow directions.
    ///
    /// Returns a list of `(row, col)` cells belonging to the watershed.
    pub fn watershed_delineation(&self, outlet: (usize, usize)) -> Vec<(usize, usize)> {
        let h = self.raster.height;
        let w = self.raster.width;
        let fd = self.flow_direction();
        // For each cell check whether its flow direction eventually reaches outlet
        // (single-step check: does the immediate downstream cell flow toward outlet?)
        let dir_to_offset = |d: u8| -> (i32, i32) {
            match d {
                1 => (0, 1),
                2 => (1, 1),
                4 => (1, 0),
                8 => (1, -1),
                16 => (0, -1),
                32 => (-1, -1),
                64 => (-1, 0),
                128 => (-1, 1),
                _ => (0, 0),
            }
        };

        let mut watershed = Vec::new();
        let mut visited = vec![vec![false; w]; h];
        let mut stack = vec![outlet];
        visited[outlet.0][outlet.1] = true;

        while let Some((r, c)) = stack.pop() {
            watershed.push((r, c));
            // Find all cells that drain directly into (r, c)
            for nr in 0..h {
                for nc in 0..w {
                    if visited[nr][nc] {
                        continue;
                    }
                    let (dr, dc) = dir_to_offset(fd[nr][nc]);
                    let tr = nr as i32 + dr;
                    let tc = nc as i32 + dc;
                    if tr == r as i32 && tc == c as i32 {
                        visited[nr][nc] = true;
                        stack.push((nr, nc));
                    }
                }
            }
        }
        watershed
    }

    /// Extract a stream network as connected channel segments.
    ///
    /// Cells with flow accumulation ≥ `threshold` are considered streams.
    /// Returns a list of segments, each segment being an ordered list of cells.
    pub fn extract_stream_network(&self, threshold: f64) -> Vec<Vec<(usize, usize)>> {
        let accum = self.flow_accumulation();
        let h = self.raster.height;
        let w = self.raster.width;
        let fd = self.flow_direction();
        let is_stream: Vec<Vec<bool>> = accum
            .iter()
            .map(|row| row.iter().map(|&v| v >= threshold).collect())
            .collect();

        let dir_to_offset = |d: u8| -> (i32, i32) {
            match d {
                1 => (0, 1),
                2 => (1, 1),
                4 => (1, 0),
                8 => (1, -1),
                16 => (0, -1),
                32 => (-1, -1),
                64 => (-1, 0),
                128 => (-1, 1),
                _ => (0, 0),
            }
        };

        let mut visited = vec![vec![false; w]; h];
        let mut segments: Vec<Vec<(usize, usize)>> = Vec::new();

        for r in 0..h {
            for c in 0..w {
                if !is_stream[r][c] || visited[r][c] {
                    continue;
                }
                // Trace downstream segment
                let mut seg = Vec::new();
                let mut cr = r;
                let mut cc = c;
                loop {
                    if visited[cr][cc] || !is_stream[cr][cc] {
                        break;
                    }
                    visited[cr][cc] = true;
                    seg.push((cr, cc));
                    let (dr, dc) = dir_to_offset(fd[cr][cc]);
                    let nr = cr as i32 + dr;
                    let nc = cc as i32 + dc;
                    if nr < 0 || nc < 0 || nr >= h as i32 || nc >= w as i32 {
                        break;
                    }
                    cr = nr as usize;
                    cc = nc as usize;
                }
                if !seg.is_empty() {
                    segments.push(seg);
                }
            }
        }
        segments
    }
}

// ---------------------------------------------------------------------------
// PointCloudIO
// ---------------------------------------------------------------------------

/// Utilities for reading and writing 3-D point cloud data.
pub struct PointCloudIO;

impl PointCloudIO {
    /// Parse a minimal LAS-style ASCII stub where each line is `x y z`.
    ///
    /// Lines beginning with `#` are treated as comments and skipped.
    pub fn read_las_stub(content: &str) -> Vec<[f64; 3]> {
        content
            .lines()
            .filter(|l| !l.trim_start().starts_with('#') && !l.trim().is_empty())
            .filter_map(|l| {
                let mut parts = l.split_whitespace();
                let x = parts.next()?.parse::<f64>().ok()?;
                let y = parts.next()?.parse::<f64>().ok()?;
                let z = parts.next()?.parse::<f64>().ok()?;
                Some([x, y, z])
            })
            .collect()
    }

    /// Serialize a slice of points to a whitespace-delimited XYZ string.
    pub fn write_xyz(points: &[[f64; 3]]) -> String {
        points
            .iter()
            .map(|p| format!("{} {} {}\n", p[0], p[1], p[2]))
            .collect()
    }

    /// Parse a whitespace-delimited XYZ text file into a vector of 3-D points.
    pub fn read_xyz(content: &str) -> Vec<[f64; 3]> {
        Self::read_las_stub(content)
    }

    /// Compute the axis-aligned bounding box `[min_x,min_y,min_z,max_x,max_y,max_z]`.
    ///
    /// Returns all zeros if `points` is empty.
    pub fn bounding_box(points: &[[f64; 3]]) -> [f64; 6] {
        if points.is_empty() {
            return [0.0; 6];
        }
        let mut mn = points[0];
        let mut mx = points[0];
        for p in points.iter().skip(1) {
            for i in 0..3 {
                mn[i] = mn[i].min(p[i]);
                mx[i] = mx[i].max(p[i]);
            }
        }
        [mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]]
    }

    /// Voxelise a point cloud at `resolution` and return one representative
    /// point per non-empty voxel (the centroid of all points in that voxel).
    pub fn voxelize(points: &[[f64; 3]], resolution: f64) -> Vec<[f64; 3]> {
        if points.is_empty() || resolution <= 0.0 {
            return vec![];
        }
        let bb = Self::bounding_box(points);
        use std::collections::HashMap;
        let mut voxels: HashMap<(i64, i64, i64), (f64, f64, f64, u64)> = HashMap::new();
        for p in points {
            let ix = ((p[0] - bb[0]) / resolution).floor() as i64;
            let iy = ((p[1] - bb[1]) / resolution).floor() as i64;
            let iz = ((p[2] - bb[2]) / resolution).floor() as i64;
            let e = voxels.entry((ix, iy, iz)).or_insert((0.0, 0.0, 0.0, 0));
            e.0 += p[0];
            e.1 += p[1];
            e.2 += p[2];
            e.3 += 1;
        }
        voxels
            .values()
            .map(|&(sx, sy, sz, n)| {
                let n = n as f64;
                [sx / n, sy / n, sz / n]
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// LAS 1.4 binary point cloud
// ---------------------------------------------------------------------------

/// Fixed-length (375-byte) header parsed from a LAS 1.4 binary file.
#[derive(Debug, Clone)]
pub struct LasHeader {
    /// File creation year (GPS week or calendar year).
    pub file_creation_year: u16,
    /// File creation day of year.
    pub file_creation_doy: u16,
    /// LAS spec major version (should be 1).
    pub version_major: u8,
    /// LAS spec minor version (should be 4 for LAS 1.4).
    pub version_minor: u8,
    /// Number of point records.
    pub point_count: u64,
    /// X scale factor.
    pub x_scale: f64,
    /// Y scale factor.
    pub y_scale: f64,
    /// Z scale factor.
    pub z_scale: f64,
    /// X offset.
    pub x_offset: f64,
    /// Y offset.
    pub y_offset: f64,
    /// Z offset.
    pub z_offset: f64,
    /// Point data format ID (0–10).
    pub point_format: u8,
    /// Number of bytes per point record.
    pub point_record_len: u16,
    /// Byte offset to the first point record.
    pub point_data_offset: u32,
}

/// One decoded LAS point (format 0 fields).
#[derive(Debug, Clone)]
pub struct LasPoint {
    /// World-space X coordinate (metres or survey feet).
    pub x: f64,
    /// World-space Y coordinate.
    pub y: f64,
    /// World-space Z coordinate.
    pub z: f64,
    /// Return intensity (0–65535).
    pub intensity: u16,
    /// Return number and flags byte.
    pub return_flags: u8,
    /// Classification byte.
    pub classification: u8,
    /// Scan angle rank.
    pub scan_angle: i8,
    /// User data.
    pub user_data: u8,
    /// Point source ID.
    pub point_source_id: u16,
}

/// A decoded LAS point cloud (header + points).
#[derive(Debug, Clone)]
pub struct LasPointCloud {
    /// Parsed LAS header.
    pub header: LasHeader,
    /// Decoded point records.
    pub points: Vec<LasPoint>,
}

// ── internal byte-reading helpers ──────────────────────────────────────────

fn las_read_u8(data: &[u8], off: &mut usize) -> Result<u8, String> {
    if *off >= data.len() {
        return Err(format!("LAS: unexpected EOF at offset {}", *off));
    }
    let v = data[*off];
    *off += 1;
    Ok(v)
}

fn las_read_u16_le(data: &[u8], off: &mut usize) -> Result<u16, String> {
    if *off + 2 > data.len() {
        return Err(format!(
            "LAS: unexpected EOF reading u16 at offset {}",
            *off
        ));
    }
    let v = u16::from_le_bytes([data[*off], data[*off + 1]]);
    *off += 2;
    Ok(v)
}

fn las_read_u32_le(data: &[u8], off: &mut usize) -> Result<u32, String> {
    if *off + 4 > data.len() {
        return Err(format!(
            "LAS: unexpected EOF reading u32 at offset {}",
            *off
        ));
    }
    let v = u32::from_le_bytes(
        data[*off..*off + 4]
            .try_into()
            .map_err(|e| format!("{e}"))?,
    );
    *off += 4;
    Ok(v)
}

fn las_read_u64_le(data: &[u8], off: &mut usize) -> Result<u64, String> {
    if *off + 8 > data.len() {
        return Err(format!(
            "LAS: unexpected EOF reading u64 at offset {}",
            *off
        ));
    }
    let v = u64::from_le_bytes(
        data[*off..*off + 8]
            .try_into()
            .map_err(|e| format!("{e}"))?,
    );
    *off += 8;
    Ok(v)
}

fn las_read_f64_le(data: &[u8], off: &mut usize) -> Result<f64, String> {
    if *off + 8 > data.len() {
        return Err(format!(
            "LAS: unexpected EOF reading f64 at offset {}",
            *off
        ));
    }
    let v = f64::from_le_bytes(
        data[*off..*off + 8]
            .try_into()
            .map_err(|e| format!("{e}"))?,
    );
    *off += 8;
    Ok(v)
}

fn las_read_i8(data: &[u8], off: &mut usize) -> Result<i8, String> {
    las_read_u8(data, off).map(|v| v as i8)
}

/// Parse a LAS 1.4 binary point cloud from raw bytes.
///
/// Validates the `"LASF"` file signature, decodes the fixed 375-byte header,
/// then reads all point records using Point Data Format 0 as the base layout
/// (the first 20 bytes of every format are identical).
///
/// # Errors
///
/// Returns a descriptive `String` error if the data is too short, the
/// signature is wrong, or any byte offset exceeds the data length.
pub fn read_las_binary(data: &[u8]) -> Result<LasPointCloud, String> {
    // LAS 1.4 header is 375 bytes minimum.
    if data.len() < 375 {
        return Err(format!(
            "LAS: data too short ({} bytes, need at least 375)",
            data.len()
        ));
    }

    // ── File signature "LASF" (bytes 0-3) ───────────────────────────────────
    if &data[0..4] != b"LASF" {
        return Err(format!(
            "LAS: bad signature {:?}, expected b\"LASF\"",
            &data[0..4]
        ));
    }

    let mut off = 4usize;

    // File source ID (u16) + Global encoding (u16)
    let _file_source_id = las_read_u16_le(data, &mut off)?;
    let _global_encoding = las_read_u16_le(data, &mut off)?;
    // Project GUID (16 bytes)
    off += 16;
    // Version major / minor
    let version_major = las_read_u8(data, &mut off)?;
    let version_minor = las_read_u8(data, &mut off)?;
    // System identifier (32 bytes)
    off += 32;
    // Generating software (32 bytes)
    off += 32;
    // File creation day / year
    let file_creation_doy = las_read_u16_le(data, &mut off)?;
    let file_creation_year = las_read_u16_le(data, &mut off)?;
    // Header size (u16) — skip, we know it's 375 for LAS 1.4
    let _header_size = las_read_u16_le(data, &mut off)?;
    // Offset to point data (u32)
    let point_data_offset = las_read_u32_le(data, &mut off)?;
    // Number of VLRs (u32) — skip VLR parsing for now
    let _num_vlrs = las_read_u32_le(data, &mut off)?;
    // Point data format (u8)
    let point_format = las_read_u8(data, &mut off)?;
    // Point data record length (u16)
    let point_record_len = las_read_u16_le(data, &mut off)?;
    // Legacy point count (u32) — LAS 1.4 uses the u64 field at offset 247
    let legacy_point_count = las_read_u32_le(data, &mut off)?;
    // Legacy point counts by return (5 × u32) — 20 bytes
    off += 20;
    // Scale factors (3 × f64)
    let x_scale = las_read_f64_le(data, &mut off)?;
    let y_scale = las_read_f64_le(data, &mut off)?;
    let z_scale = las_read_f64_le(data, &mut off)?;
    // Offsets (3 × f64)
    let x_offset = las_read_f64_le(data, &mut off)?;
    let y_offset = las_read_f64_le(data, &mut off)?;
    let z_offset = las_read_f64_le(data, &mut off)?;
    // Max/min X, Y, Z (6 × f64) = 48 bytes
    off += 48;
    // Start of Waveform Data Packet Record (u64) + Start of first EVLR (u64) +
    // Number of EVLRs (u32) = 20 bytes — skip, not needed for point reading
    off += 20;

    // LAS 1.4 extended point count (u64) at offset 247.
    let point_count_extended = las_read_u64_le(data, &mut off)?;
    // Extended point counts by return (15 × u64) = 120 bytes — skip, not needed
    let _ = off; // header parsing complete

    // Resolve point count: prefer non-zero extended count.
    let point_count = if point_count_extended > 0 {
        point_count_extended
    } else {
        legacy_point_count as u64
    };

    let header = LasHeader {
        file_creation_year,
        file_creation_doy,
        version_major,
        version_minor,
        point_count,
        x_scale,
        y_scale,
        z_scale,
        x_offset,
        y_offset,
        z_offset,
        point_format,
        point_record_len,
        point_data_offset,
    };

    // ── Point records ────────────────────────────────────────────────────────
    let mut points: Vec<LasPoint> = Vec::with_capacity(point_count as usize);
    let rec_len = point_record_len.max(20) as usize;
    let mut p_off = point_data_offset as usize;

    for _ in 0..point_count {
        if p_off + rec_len > data.len() {
            break;
        }
        let mut o = p_off;
        let raw_x = i32::from_le_bytes(
            data[o..o + 4]
                .try_into()
                .map_err(|e| format!("LAS point x: {e}"))?,
        );
        o += 4;
        let raw_y = i32::from_le_bytes(
            data[o..o + 4]
                .try_into()
                .map_err(|e| format!("LAS point y: {e}"))?,
        );
        o += 4;
        let raw_z = i32::from_le_bytes(
            data[o..o + 4]
                .try_into()
                .map_err(|e| format!("LAS point z: {e}"))?,
        );
        o += 4;
        let intensity = las_read_u16_le(data, &mut o)?;
        let return_flags = las_read_u8(data, &mut o)?;
        let classification = las_read_u8(data, &mut o)?;
        let scan_angle = las_read_i8(data, &mut o)?;
        let user_data = las_read_u8(data, &mut o)?;
        let point_source_id = las_read_u16_le(data, &mut o)?;

        points.push(LasPoint {
            x: raw_x as f64 * x_scale + x_offset,
            y: raw_y as f64 * y_scale + y_offset,
            z: raw_z as f64 * z_scale + z_offset,
            intensity,
            return_flags,
            classification,
            scan_angle,
            user_data,
            point_source_id,
        });
        p_off += rec_len;
    }

    Ok(LasPointCloud { header, points })
}

/// Write a minimal LAS 1.4 binary file from a slice of world-space points.
///
/// Uses Point Data Format 0 with a 375-byte header and a scale factor of
/// `scale` (default 0.001 m = 1 mm precision) and no offset.  Intended only
/// for round-trip testing of [`read_las_binary`].
pub fn write_las_binary_minimal(points: &[[f64; 3]], scale: f64) -> Vec<u8> {
    let n = points.len() as u64;
    let rec_len: u16 = 20; // Point Format 0 minimum
    let header_size: u16 = 375;
    let point_data_offset: u32 = header_size as u32;
    let mut buf = Vec::with_capacity(header_size as usize + n as usize * rec_len as usize);

    // File signature
    buf.extend_from_slice(b"LASF");
    // File source ID + global encoding
    buf.extend_from_slice(&0u16.to_le_bytes());
    buf.extend_from_slice(&0u16.to_le_bytes());
    // Project GUID (16 zero bytes)
    buf.extend_from_slice(&[0u8; 16]);
    // Version: 1.4
    buf.push(1u8);
    buf.push(4u8);
    // System identifier (32 bytes)
    buf.extend_from_slice(&[0u8; 32]);
    // Generating software (32 bytes) — fill with "oxiphysics"
    let mut sw = [b' '; 32];
    let tag = b"oxiphysics";
    sw[..tag.len()].copy_from_slice(tag);
    buf.extend_from_slice(&sw);
    // File creation DOY / year
    buf.extend_from_slice(&1u16.to_le_bytes()); // DOY
    buf.extend_from_slice(&2026u16.to_le_bytes()); // year
    // Header size
    buf.extend_from_slice(&header_size.to_le_bytes());
    // Offset to point data
    buf.extend_from_slice(&point_data_offset.to_le_bytes());
    // Number of VLRs
    buf.extend_from_slice(&0u32.to_le_bytes());
    // Point data format
    buf.push(0u8); // format 0
    // Point data record length
    buf.extend_from_slice(&rec_len.to_le_bytes());
    // Legacy point count
    let legacy_n = n.min(u32::MAX as u64) as u32;
    buf.extend_from_slice(&legacy_n.to_le_bytes());
    // Legacy point counts by return (5 × u32)
    buf.extend_from_slice(&[0u8; 20]);
    // Scale factors
    buf.extend_from_slice(&scale.to_le_bytes());
    buf.extend_from_slice(&scale.to_le_bytes());
    buf.extend_from_slice(&scale.to_le_bytes());
    // Offsets (0, 0, 0)
    buf.extend_from_slice(&0.0f64.to_le_bytes());
    buf.extend_from_slice(&0.0f64.to_le_bytes());
    buf.extend_from_slice(&0.0f64.to_le_bytes());
    // Max/min X, Y, Z — compute from points
    let (mut max_x, mut min_x) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut max_y, mut min_y) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut max_z, mut min_z) = (f64::NEG_INFINITY, f64::INFINITY);
    for &[x, y, z] in points {
        if x > max_x {
            max_x = x;
        }
        if x < min_x {
            min_x = x;
        }
        if y > max_y {
            max_y = y;
        }
        if y < min_y {
            min_y = y;
        }
        if z > max_z {
            max_z = z;
        }
        if z < min_z {
            min_z = z;
        }
    }
    if points.is_empty() {
        max_x = 0.0;
        min_x = 0.0;
        max_y = 0.0;
        min_y = 0.0;
        max_z = 0.0;
        min_z = 0.0;
    }
    buf.extend_from_slice(&max_x.to_le_bytes());
    buf.extend_from_slice(&min_x.to_le_bytes());
    buf.extend_from_slice(&max_y.to_le_bytes());
    buf.extend_from_slice(&min_y.to_le_bytes());
    buf.extend_from_slice(&max_z.to_le_bytes());
    buf.extend_from_slice(&min_z.to_le_bytes());
    // Start of Waveform Data Packet Record (u64)
    buf.extend_from_slice(&0u64.to_le_bytes());
    // Start of first EVLR (u64)
    buf.extend_from_slice(&0u64.to_le_bytes());
    // Number of EVLRs (u32)
    buf.extend_from_slice(&0u32.to_le_bytes());
    // LAS 1.4 extended point count (u64)
    buf.extend_from_slice(&n.to_le_bytes());
    // Extended point counts by return (15 × u64)
    buf.extend_from_slice(&[0u8; 120]);

    // Sanity-check: buf.len() must equal header_size at this point.
    assert_eq!(buf.len(), header_size as usize, "header length mismatch");

    // ── Point records (Point Format 0) ────────────────────────────────────────
    for &[x, y, z] in points {
        let raw_x = ((x) / scale).round() as i32;
        let raw_y = ((y) / scale).round() as i32;
        let raw_z = ((z) / scale).round() as i32;
        buf.extend_from_slice(&raw_x.to_le_bytes());
        buf.extend_from_slice(&raw_y.to_le_bytes());
        buf.extend_from_slice(&raw_z.to_le_bytes());
        buf.extend_from_slice(&0u16.to_le_bytes()); // intensity
        buf.push(0u8); // return flags
        buf.push(0u8); // classification
        buf.push(0u8); // scan angle (i8 → u8 cast fine for 0)
        buf.push(0u8); // user data
        buf.extend_from_slice(&0u16.to_le_bytes()); // point source ID
    }

    buf
}

// ---------------------------------------------------------------------------
// SatelliteImage
// ---------------------------------------------------------------------------

/// Multi-spectral satellite image with named bands.
#[derive(Debug, Clone)]
pub struct SatelliteImage {
    /// Flattened (row-major) pixel values per band.
    pub bands: Vec<Vec<f64>>,
    /// Human-readable names for each band (e.g. `"red"`, `"nir"`).
    pub band_names: Vec<String>,
}

impl SatelliteImage {
    /// Construct a `SatelliteImage` from per-band pixel data and band names.
    pub fn new(bands: Vec<Vec<f64>>, band_names: Vec<String>) -> Self {
        Self { bands, band_names }
    }

    /// Return the index of the band with the given name, if present.
    pub fn band_index(&self, name: &str) -> Option<usize> {
        self.band_names.iter().position(|n| n == name)
    }

    /// Compute the Normalized Difference Vegetation Index (NDVI).
    ///
    /// Requires bands named `"nir"` and `"red"`.  Returns a vector of per-pixel
    /// NDVI values in `[-1, 1]`.
    pub fn ndvi(&self) -> Vec<f64> {
        let nir_idx = self.band_index("nir").unwrap_or(0);
        let red_idx = self.band_index("red").unwrap_or(1);
        if self.bands.len() <= nir_idx.max(red_idx) {
            return vec![];
        }
        let nir = &self.bands[nir_idx];
        let red = &self.bands[red_idx];
        nir.iter()
            .zip(red.iter())
            .map(|(&n, &r)| {
                let denom = n + r;
                if denom.abs() < 1e-10 {
                    0.0
                } else {
                    (n - r) / denom
                }
            })
            .collect()
    }

    /// Compute the Normalized Difference Water Index (NDWI).
    ///
    /// Requires bands named `"green"` and `"nir"`.  Returns per-pixel values
    /// in `[-1, 1]`.
    pub fn ndwi(&self) -> Vec<f64> {
        let green_idx = self.band_index("green").unwrap_or(0);
        let nir_idx = self.band_index("nir").unwrap_or(1);
        if self.bands.len() <= green_idx.max(nir_idx) {
            return vec![];
        }
        let green = &self.bands[green_idx];
        let nir = &self.bands[nir_idx];
        green
            .iter()
            .zip(nir.iter())
            .map(|(&g, &n)| {
                let denom = g + n;
                if denom.abs() < 1e-10 {
                    0.0
                } else {
                    (g - n) / denom
                }
            })
            .collect()
    }

    /// Produce an RGB false-colour composite by selecting bands `r`, `g`, `b`
    /// (0-indexed) and normalising each channel to `[0, 1]`.
    pub fn false_color_composite(&self, r: usize, g: usize, b: usize) -> Vec<[f64; 3]> {
        if self.bands.len() <= r.max(g).max(b) {
            return vec![];
        }
        let normalize = |band: &Vec<f64>| -> Vec<f64> {
            let mn = band.iter().cloned().fold(f64::INFINITY, f64::min);
            let mx = band.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let range = mx - mn;
            if range < 1e-10 {
                vec![0.0; band.len()]
            } else {
                band.iter().map(|&v| (v - mn) / range).collect()
            }
        };
        let rb = normalize(&self.bands[r]);
        let gb = normalize(&self.bands[g]);
        let bb = normalize(&self.bands[b]);
        rb.iter()
            .zip(gb.iter())
            .zip(bb.iter())
            .map(|((&rv, &gv), &bv)| [rv, gv, bv])
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_raster() -> RasterData {
        // 5×5 DEM with a simple gradient
        let data: Vec<Vec<f64>> = (0..5)
            .map(|r| (0..5).map(|c| (r * 5 + c) as f64 * 10.0).collect())
            .collect();
        RasterData::new(data, -9999.0)
    }

    fn sample_dem() -> Dem {
        Dem::new(sample_raster(), 30.0)
    }

    // --- GeoTiffReader ---

    #[test]
    fn test_geotiff_pixel_to_world_origin() {
        let t = [100.0, 0.5, 0.0, 50.0, 0.0, -0.5];
        let reader = GeoTiffReader::new(10, 10, 1, t, vec![vec![1.0; 100]]);
        let [lon, lat] = reader.pixel_to_world(0, 0);
        assert!((lon - 100.0).abs() < 1e-9);
        assert!((lat - 50.0).abs() < 1e-9);
    }

    #[test]
    fn test_geotiff_pixel_to_world_offset() {
        let t = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let reader = GeoTiffReader::new(4, 4, 1, t, vec![vec![0.0; 16]]);
        let [lon, lat] = reader.pixel_to_world(2, 3);
        assert!((lon - 2.0).abs() < 1e-9);
        assert!((lat - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_geotiff_world_to_pixel_roundtrip() {
        let t = [10.0, 0.1, 0.0, 20.0, 0.0, -0.1];
        let data = vec![vec![5.0; 100]];
        let reader = GeoTiffReader::new(10, 10, 1, t, data);
        let (px, py) = reader.world_to_pixel(10.5, 19.5);
        assert_eq!(px, 5);
        assert_eq!(py, 5);
    }

    #[test]
    fn test_geotiff_elevation_at() {
        let t = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let band: Vec<f64> = (0..16).map(|i| i as f64).collect();
        let reader = GeoTiffReader::new(4, 4, 1, t, vec![band]);
        // pixel (0,0) → index 0 → value 0.0
        let elev = reader.elevation_at(0.0, 0.0);
        assert!((elev - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_geotiff_bounding_box_size() {
        let t = [10.0, 1.0, 0.0, 50.0, 0.0, -1.0];
        let reader = GeoTiffReader::new(5, 5, 1, t, vec![vec![0.0; 25]]);
        let bb = reader.bounding_box();
        assert!(bb[2] > bb[0]);
    }

    #[test]
    fn test_geotiff_read_band_empty_when_out_of_range() {
        let reader = GeoTiffReader::new(2, 2, 1, [0.0; 6], vec![vec![1.0; 4]]);
        assert!(reader.read_band(5).is_empty());
    }

    // --- RasterData ---

    #[test]
    fn test_raster_statistics_basic() {
        let r = sample_raster();
        let (mn, mx, mean, _std) = r.statistics();
        assert!((mn - 0.0).abs() < 1e-9);
        assert!((mx - 240.0).abs() < 1e-9);
        assert!(mean > 0.0);
    }

    #[test]
    fn test_raster_statistics_all_nodata() {
        let r = RasterData::new(vec![vec![-9999.0; 3]; 3], -9999.0);
        let (mn, mx, mean, std) = r.statistics();
        assert_eq!((mn, mx, mean, std), (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_raster_clip_to_bbox_dimensions() {
        let r = sample_raster();
        let clipped = r.clip_to_bbox([1.0, 1.0, 4.0, 4.0]);
        assert_eq!(clipped.width, 3);
        assert_eq!(clipped.height, 3);
    }

    #[test]
    fn test_raster_resample_down() {
        let r = sample_raster();
        let small = r.resample(0.5);
        assert!(small.width <= r.width);
        assert!(small.height <= r.height);
    }

    #[test]
    fn test_raster_resample_up() {
        let r = sample_raster();
        let big = r.resample(2.0);
        assert!(big.width >= r.width);
        assert!(big.height >= r.height);
    }

    #[test]
    fn test_raster_slope_border_zero() {
        let r = sample_raster();
        let s = r.slope();
        // Border pixels not computed (remain 0)
        assert!((s[0][0] - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_raster_slope_interior_positive() {
        let r = sample_raster();
        let s = r.slope();
        assert!(s[2][2] >= 0.0);
    }

    #[test]
    fn test_raster_aspect_range() {
        let r = sample_raster();
        let a = r.aspect();
        for row in &a {
            for &v in row {
                assert!(v >= 0.0 && v < 360.0 + 1e-9, "aspect out of range: {v}");
            }
        }
    }

    #[test]
    fn test_raster_hillshade_range() {
        let r = sample_raster();
        let hs = r.hillshade(315.0, 45.0);
        for row in &hs {
            for &v in row {
                assert!(v >= 0.0 && v <= 1.0 + 1e-9, "hillshade out of [0,1]: {v}");
            }
        }
    }

    // --- Dem ---

    #[test]
    fn test_dem_flow_direction_shape() {
        let d = sample_dem();
        let fd = d.flow_direction();
        assert_eq!(fd.len(), 5);
        assert_eq!(fd[0].len(), 5);
    }

    #[test]
    fn test_dem_flow_accumulation_shape() {
        let d = sample_dem();
        let fa = d.flow_accumulation();
        assert_eq!(fa.len(), 5);
        assert_eq!(fa[0].len(), 5);
    }

    #[test]
    fn test_dem_flow_accumulation_positive() {
        let d = sample_dem();
        let fa = d.flow_accumulation();
        for row in &fa {
            for &v in row {
                assert!(v >= 1.0, "flow accumulation must be ≥ 1");
            }
        }
    }

    #[test]
    fn test_dem_watershed_contains_outlet() {
        let d = sample_dem();
        let ws = d.watershed_delineation((4, 4));
        assert!(ws.contains(&(4, 4)));
    }

    #[test]
    fn test_dem_stream_network_threshold_high() {
        let d = sample_dem();
        // Very high threshold → no streams
        let segs = d.extract_stream_network(1e12);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_dem_stream_network_threshold_low() {
        let d = sample_dem();
        // Low threshold → all cells are streams
        let segs = d.extract_stream_network(0.0);
        assert!(!segs.is_empty());
    }

    // --- PointCloudIO ---

    #[test]
    fn test_point_cloud_write_read_roundtrip() {
        let pts = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let s = PointCloudIO::write_xyz(&pts);
        let back = PointCloudIO::read_xyz(&s);
        assert_eq!(back.len(), 2);
        assert!((back[0][0] - 1.0).abs() < 1e-9);
        assert!((back[1][2] - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_point_cloud_bounding_box_empty() {
        let bb = PointCloudIO::bounding_box(&[]);
        assert_eq!(bb, [0.0; 6]);
    }

    #[test]
    fn test_point_cloud_bounding_box_single() {
        let bb = PointCloudIO::bounding_box(&[[3.0, 1.0, -2.0]]);
        assert!((bb[0] - 3.0).abs() < 1e-9);
        assert!((bb[3] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_point_cloud_bounding_box_multiple() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, 5.0, -2.0]];
        let bb = PointCloudIO::bounding_box(&pts);
        assert!((bb[0] - (-1.0)).abs() < 1e-9); // min_x
        assert!((bb[4] - 5.0).abs() < 1e-9); // max_y
    }

    #[test]
    fn test_point_cloud_voxelize_reduces_count() {
        let pts: Vec<[f64; 3]> = (0..100)
            .map(|i| [(i % 10) as f64 * 0.1, (i / 10) as f64 * 0.1, 0.0])
            .collect();
        let vox = PointCloudIO::voxelize(&pts, 0.5);
        assert!(vox.len() < pts.len());
    }

    #[test]
    fn test_point_cloud_las_stub_comments_skipped() {
        let content = "# header\n1 2 3\n# comment\n4 5 6\n";
        let pts = PointCloudIO::read_las_stub(content);
        assert_eq!(pts.len(), 2);
    }

    // --- SatelliteImage ---

    fn make_image() -> SatelliteImage {
        let red: Vec<f64> = (0..4).map(|i| i as f64 * 0.1).collect();
        let nir: Vec<f64> = (0..4).map(|i| i as f64 * 0.2 + 0.1).collect();
        let green: Vec<f64> = (0..4).map(|i| i as f64 * 0.15).collect();
        SatelliteImage::new(
            vec![red, nir, green],
            vec!["red".into(), "nir".into(), "green".into()],
        )
    }

    #[test]
    fn test_satellite_ndvi_length() {
        let img = make_image();
        assert_eq!(img.ndvi().len(), 4);
    }

    #[test]
    fn test_satellite_ndvi_range() {
        let img = make_image();
        for v in img.ndvi() {
            assert!(
                v >= -1.0 - 1e-9 && v <= 1.0 + 1e-9,
                "NDVI out of range: {v}"
            );
        }
    }

    #[test]
    fn test_satellite_ndwi_length() {
        let img = make_image();
        assert_eq!(img.ndwi().len(), 4);
    }

    #[test]
    fn test_satellite_false_color_length() {
        let img = make_image();
        let fc = img.false_color_composite(0, 2, 1);
        assert_eq!(fc.len(), 4);
    }

    #[test]
    fn test_satellite_false_color_range() {
        let img = make_image();
        for px in img.false_color_composite(0, 2, 1) {
            for &ch in &px {
                assert!(
                    ch >= 0.0 - 1e-9 && ch <= 1.0 + 1e-9,
                    "channel out of [0,1]: {ch}"
                );
            }
        }
    }

    #[test]
    fn test_satellite_band_index_found() {
        let img = make_image();
        assert_eq!(img.band_index("nir"), Some(1));
    }

    #[test]
    fn test_satellite_band_index_not_found() {
        let img = make_image();
        assert_eq!(img.band_index("swir"), None);
    }

    // ── G5: LAS binary reader tests ───────────────────────────────────────────

    #[test]
    fn test_read_las_bad_signature() {
        let data = b"NOTL\x00\x00\x00\x00";
        assert!(read_las_binary(data).is_err());
    }

    #[test]
    fn test_read_las_too_short() {
        let data = vec![0u8; 100];
        assert!(read_las_binary(&data).is_err());
    }

    #[test]
    fn test_las_roundtrip_empty() {
        let data = write_las_binary_minimal(&[], 0.001);
        let cloud = read_las_binary(&data).expect("should parse empty LAS");
        assert_eq!(cloud.points.len(), 0);
        assert_eq!(cloud.header.version_major, 1);
        assert_eq!(cloud.header.version_minor, 4);
    }

    #[test]
    fn test_las_roundtrip_single_point() {
        let pts = vec![[1.234, 5.678, -0.5]];
        let data = write_las_binary_minimal(&pts, 0.001);
        let cloud = read_las_binary(&data).expect("should parse single-point LAS");
        assert_eq!(cloud.points.len(), 1);
        assert!((cloud.points[0].x - 1.234).abs() < 0.002, "x mismatch");
        assert!((cloud.points[0].y - 5.678).abs() < 0.002, "y mismatch");
        assert!((cloud.points[0].z - (-0.5)).abs() < 0.002, "z mismatch");
    }

    #[test]
    fn test_las_roundtrip_multiple_points() {
        let pts: Vec<[f64; 3]> = (0..10)
            .map(|i| [i as f64 * 0.5, i as f64 * 0.25, i as f64 * 0.1])
            .collect();
        let data = write_las_binary_minimal(&pts, 0.001);
        let cloud = read_las_binary(&data).expect("should parse multi-point LAS");
        assert_eq!(cloud.points.len(), 10);
        for (i, pt) in cloud.points.iter().enumerate() {
            let expected_x = i as f64 * 0.5;
            assert!(
                (pt.x - expected_x).abs() < 0.002,
                "x mismatch at point {i}: expected {expected_x}, got {}",
                pt.x
            );
        }
    }
}
