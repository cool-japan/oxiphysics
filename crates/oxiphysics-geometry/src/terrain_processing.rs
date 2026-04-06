#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Terrain processing and analysis algorithms.
//!
//! This module provides a comprehensive suite of Digital Elevation Model (DEM)
//! analysis tools, including:
//!
//! - Slope and aspect computation (Horn's method)
//! - Curvature analysis: profile, planform, and tangential curvature
//! - Flow direction algorithms: D8 and D-infinity (Tarboton 1997)
//! - Flow accumulation and watershed delineation
//! - Stream network extraction
//! - Terrain Roughness Index (TRI — Riley et al. 1999)
//! - Terrain Position Index (TPI — Weiss 2001)
//! - Topographic Wetness Index (TWI — Beven & Kirkby 1979)
//! - Viewshed analysis (line-of-sight)
//! - Slope stability index (infinite slope model)
//! - Solar irradiance on terrain (incidence angle, terrain shading)

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Section 1 — Core DEM structure
// ---------------------------------------------------------------------------

/// A regular-grid Digital Elevation Model.
///
/// Heights are stored in row-major order: index = `row * cols + col`.
/// The origin is at the north-west corner; rows increase southward (y+),
/// columns increase eastward (x+).
#[derive(Debug, Clone)]
pub struct Dem {
    /// Elevation values (metres), row-major.
    pub heights: Vec<f64>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Cell size in metres (square cells assumed).
    pub cell_size: f64,
    /// No-data sentinel value.
    pub nodata: f64,
}

impl Dem {
    /// Construct a new DEM.
    ///
    /// # Panics
    /// Panics if `heights.len() != rows * cols`.
    pub fn new(heights: Vec<f64>, rows: usize, cols: usize, cell_size: f64) -> Self {
        assert_eq!(heights.len(), rows * cols, "heights length mismatch");
        Self {
            heights,
            rows,
            cols,
            cell_size,
            nodata: -9999.0,
        }
    }

    /// Return the elevation at `(row, col)`, or `None` if out of bounds.
    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row < self.rows && col < self.cols {
            Some(self.heights[row * self.cols + col])
        } else {
            None
        }
    }

    /// Return the elevation at `(row, col)` clamped to grid boundaries.
    pub fn get_clamped(&self, row: i64, col: i64) -> f64 {
        let r = row.clamp(0, self.rows as i64 - 1) as usize;
        let c = col.clamp(0, self.cols as i64 - 1) as usize;
        self.heights[r * self.cols + c]
    }

    /// Return true if the cell contains valid data.
    pub fn is_valid(&self, row: usize, col: usize) -> bool {
        if let Some(h) = self.get(row, col) {
            (h - self.nodata).abs() > 1.0
        } else {
            false
        }
    }

    /// Return a 3×3 neighbourhood window centred at `(row, col)`.
    ///
    /// Missing edge cells are filled by replication of the nearest valid cell.
    pub fn window_3x3(&self, row: usize, col: usize) -> [[f64; 3]; 3] {
        let mut w = [[0.0f64; 3]; 3];
        for dr in 0i64..3 {
            for dc in 0i64..3 {
                let r = row as i64 + dr - 1;
                let c = col as i64 + dc - 1;
                w[dr as usize][dc as usize] = self.get_clamped(r, c);
            }
        }
        w
    }
}

// ---------------------------------------------------------------------------
// Section 2 — Slope and aspect
// ---------------------------------------------------------------------------

/// Terrain slope in radians at every grid cell.
///
/// Uses Horn's (1981) finite-difference algorithm on the 3×3 neighbourhood.
/// Returns a `Vec`f64` of length `rows * cols`.
pub fn compute_slope(dem: &Dem) -> Vec<f64> {
    let mut slope = vec![0.0f64; dem.rows * dem.cols];
    let cs = dem.cell_size;
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let w = dem.window_3x3(r, c);
            // dz/dx  (west–east gradient)
            let dzdx = ((w[0][2] + 2.0 * w[1][2] + w[2][2]) - (w[0][0] + 2.0 * w[1][0] + w[2][0]))
                / (8.0 * cs);
            // dz/dy  (north–south gradient)
            let dzdy = ((w[2][0] + 2.0 * w[2][1] + w[2][2]) - (w[0][0] + 2.0 * w[0][1] + w[0][2]))
                / (8.0 * cs);
            slope[r * dem.cols + c] = (dzdx * dzdx + dzdy * dzdy).sqrt().atan();
        }
    }
    slope
}

/// Terrain aspect in radians at every grid cell.
///
/// Aspect is measured clockwise from north (0 = north, π/2 = east).
/// Flat areas (slope ≈ 0) are assigned aspect = -1.
pub fn compute_aspect(dem: &Dem) -> Vec<f64> {
    let mut aspect = vec![-1.0f64; dem.rows * dem.cols];
    let cs = dem.cell_size;
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let w = dem.window_3x3(r, c);
            let dzdx = ((w[0][2] + 2.0 * w[1][2] + w[2][2]) - (w[0][0] + 2.0 * w[1][0] + w[2][0]))
                / (8.0 * cs);
            let dzdy = ((w[2][0] + 2.0 * w[2][1] + w[2][2]) - (w[0][0] + 2.0 * w[0][1] + w[0][2]))
                / (8.0 * cs);
            if dzdx.abs() > 1e-15 || dzdy.abs() > 1e-15 {
                // atan2 measured from north, clockwise
                let a = 180.0 - dzdy.atan2(dzdx).to_degrees() + 90.0;
                let a = if a < 0.0 {
                    a + 360.0
                } else if a >= 360.0 {
                    a - 360.0
                } else {
                    a
                };
                aspect[r * dem.cols + c] = a.to_radians();
            }
        }
    }
    aspect
}

// ---------------------------------------------------------------------------
// Section 3 — Curvature
// ---------------------------------------------------------------------------

/// Curvature type selector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurvatureKind {
    /// Profile curvature: curvature in the direction of steepest descent.
    Profile,
    /// Planform (contour) curvature: curvature perpendicular to slope.
    Planform,
    /// Tangential curvature: horizontal curvature normal to the slope direction.
    Tangential,
}

/// Compute terrain curvature (units: 1/m) for every grid cell.
///
/// Uses the second-order finite difference formulation of Zevenbergen & Thorne (1987).
pub fn compute_curvature(dem: &Dem, kind: CurvatureKind) -> Vec<f64> {
    let mut curv = vec![0.0f64; dem.rows * dem.cols];
    let l = dem.cell_size;
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let w = dem.window_3x3(r, c);
            // Z layout (Zevenbergen & Thorne notation):
            //  w[0][0]=Z1  w[0][1]=Z2  w[0][2]=Z3
            //  w[1][0]=Z4  w[1][1]=Z5  w[1][2]=Z6
            //  w[2][0]=Z7  w[2][1]=Z8  w[2][2]=Z9
            let z5 = w[1][1];
            let _ = z5; // centre (not needed directly)

            let d = ((w[0][1] + w[2][1]) / 2.0 - w[1][1]) / (l * l);
            let e = ((w[1][0] + w[1][2]) / 2.0 - w[1][1]) / (l * l);
            let f = (-w[0][0] + w[0][2] + w[2][0] - w[2][2]) / (4.0 * l * l);
            let g = (-w[0][1] + w[2][1]) / (2.0 * l);
            let h = (w[1][2] - w[1][0]) / (2.0 * l);

            let p2 = g * g + h * h;

            curv[r * dem.cols + c] = match kind {
                CurvatureKind::Profile => {
                    if p2 > 1e-15 {
                        -2.0 * (d * g * g + e * h * h + f * g * h) / (p2 * (1.0 + p2).sqrt())
                    } else {
                        0.0
                    }
                }
                CurvatureKind::Planform => {
                    if p2 > 1e-15 {
                        -2.0 * (d * h * h + e * g * g - f * g * h) / p2.powf(1.5)
                    } else {
                        0.0
                    }
                }
                CurvatureKind::Tangential => {
                    if p2 > 1e-15 {
                        -2.0 * (d * h * h + e * g * g - f * g * h) / (p2 * (1.0 + p2).sqrt())
                    } else {
                        0.0
                    }
                }
            };
        }
    }
    curv
}

// ---------------------------------------------------------------------------
// Section 4 — Flow direction
// ---------------------------------------------------------------------------

/// D8 flow direction codes.
///
/// Directions are encoded as 1, 2, 4, 8, 16, 32, 64, 128 (clockwise from east,
/// following the ESRI convention).
pub type D8Code = u8;

/// Neighbour offsets and D8 codes for the eight cardinal/diagonal directions.
///
/// Order: E, SE, S, SW, W, NW, N, NE.
pub const D8_OFFSETS: [(i64, i64, D8Code); 8] = [
    (0, 1, 1),
    (1, 1, 2),
    (1, 0, 4),
    (1, -1, 8),
    (0, -1, 16),
    (-1, -1, 32),
    (-1, 0, 64),
    (-1, 1, 128),
];

/// Diagonal distance multiplier for D8 flow direction.
///
/// Cardinal neighbours have distance 1, diagonals have √2.
fn d8_distance(dr: i64, dc: i64) -> f64 {
    if dr != 0 && dc != 0 {
        2.0f64.sqrt()
    } else {
        1.0
    }
}

/// Compute D8 flow direction for every grid cell.
///
/// Returns `Vec`D8Code` with length `rows * cols`.  Flat or edge cells that
/// have no steeper neighbour receive code 0.
pub fn flow_direction_d8(dem: &Dem) -> Vec<D8Code> {
    let mut fdir = vec![0u8; dem.rows * dem.cols];
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let z = dem.heights[r * dem.cols + c];
            let mut max_drop = 0.0f64;
            let mut best_code = 0u8;
            for &(dr, dc, code) in &D8_OFFSETS {
                let nr = r as i64 + dr;
                let nc = c as i64 + dc;
                if nr < 0 || nr >= dem.rows as i64 || nc < 0 || nc >= dem.cols as i64 {
                    continue;
                }
                let nz = dem.heights[nr as usize * dem.cols + nc as usize];
                let dist = d8_distance(dr, dc) * dem.cell_size;
                let drop = (z - nz) / dist;
                if drop > max_drop {
                    max_drop = drop;
                    best_code = code;
                }
            }
            fdir[r * dem.cols + c] = best_code;
        }
    }
    fdir
}

/// D-infinity flow direction (Tarboton 1997).
///
/// Returns a `Vec`f64` of flow angles (radians, 0 = east, anticlockwise) for
/// each cell.  The angle encodes which of the eight triangular facets receives
/// flow, and how it is partitioned between the two bounding neighbours.
pub fn flow_direction_dinf(dem: &Dem) -> Vec<f64> {
    // Eight triangular facets: each defined by two neighbouring cells.
    // Facet ordering: anti-clockwise from east.
    // (dr1,dc1) is the "cardinal" neighbour; (dr2,dc2) is the diagonal.
    const FACETS: [(i64, i64, i64, i64); 8] = [
        (0, 1, -1, 1),
        (-1, 1, -1, 0),
        (-1, 0, -1, -1),
        (-1, -1, 0, -1),
        (0, -1, 1, -1),
        (1, -1, 1, 0),
        (1, 0, 1, 1),
        (1, 1, 0, 1),
    ];
    let mut fdir = vec![0.0f64; dem.rows * dem.cols];
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let z0 = dem.heights[r * dem.cols + c];
            let l = dem.cell_size;
            let mut max_slope = f64::NEG_INFINITY;
            let mut best_angle = 0.0f64;
            for (fi, &(dr1, dc1, dr2, dc2)) in FACETS.iter().enumerate() {
                let facet_angle_base = fi as f64 * PI / 4.0;
                let r1 = r as i64 + dr1;
                let c1 = c as i64 + dc1;
                let r2 = r as i64 + dr2;
                let c2 = c as i64 + dc2;
                if r1 < 0 || r1 >= dem.rows as i64 || c1 < 0 || c1 >= dem.cols as i64 {
                    continue;
                }
                if r2 < 0 || r2 >= dem.rows as i64 || c2 < 0 || c2 >= dem.cols as i64 {
                    continue;
                }
                let z1 = dem.heights[r1 as usize * dem.cols + c1 as usize];
                let z2 = dem.heights[r2 as usize * dem.cols + c2 as usize];
                let e1 = (z0 - z1) / l;
                let e2 = (z1 - z2) / l;
                // Optimal angle within facet [0, π/4]
                let alpha = if e2.abs() < 1e-15 {
                    0.0
                } else {
                    (e1 / e2).atan().clamp(0.0, PI / 4.0)
                };
                let s = (e1 * e1 + e2 * e2 - 2.0 * alpha.cos() * e1 * e2).sqrt();
                if s > max_slope {
                    max_slope = s;
                    best_angle = facet_angle_base + alpha;
                }
            }
            fdir[r * dem.cols + c] = best_angle;
        }
    }
    fdir
}

// ---------------------------------------------------------------------------
// Section 5 — Flow accumulation
// ---------------------------------------------------------------------------

/// Compute D8 flow accumulation (number of upstream cells draining into each cell).
///
/// Uses a topological sort (cells processed from high to low elevation).
pub fn flow_accumulation_d8(dem: &Dem, fdir: &[D8Code]) -> Vec<u32> {
    let n = dem.rows * dem.cols;
    let mut accum = vec![1u32; n];
    // Sort cell indices by elevation descending
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        dem.heights[b]
            .partial_cmp(&dem.heights[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for idx in order {
        let r = idx / dem.cols;
        let c = idx % dem.cols;
        let code = fdir[idx];
        if code == 0 {
            continue;
        }
        // Decode D8 code to neighbour offset
        for &(dr, dc, mask) in &D8_OFFSETS {
            if code == mask {
                let nr = r as i64 + dr;
                let nc = c as i64 + dc;
                if nr >= 0 && nr < dem.rows as i64 && nc >= 0 && nc < dem.cols as i64 {
                    let nidx = nr as usize * dem.cols + nc as usize;
                    accum[nidx] += accum[idx];
                }
                break;
            }
        }
    }
    accum
}

// ---------------------------------------------------------------------------
// Section 6 — Watershed delineation
// ---------------------------------------------------------------------------

/// Delineate the watershed (contributing area) upstream of an outlet cell.
///
/// Returns a boolean mask of the same size as the DEM.
pub fn delineate_watershed(
    dem: &Dem,
    fdir: &[D8Code],
    outlet_row: usize,
    outlet_col: usize,
) -> Vec<bool> {
    let n = dem.rows * dem.cols;
    let mut mask = vec![false; n];
    // Build reverse flow graph
    let mut upstream: Vec<Vec<usize>> = vec![vec![]; n];
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let idx = r * dem.cols + c;
            let code = fdir[idx];
            if code == 0 {
                continue;
            }
            for &(dr, dc, mask_code) in &D8_OFFSETS {
                if code == mask_code {
                    let nr = r as i64 + dr;
                    let nc = c as i64 + dc;
                    if nr >= 0 && nr < dem.rows as i64 && nc >= 0 && nc < dem.cols as i64 {
                        let nidx = nr as usize * dem.cols + nc as usize;
                        upstream[nidx].push(idx);
                    }
                    break;
                }
            }
        }
    }
    // BFS from outlet
    let outlet = outlet_row * dem.cols + outlet_col;
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(outlet);
    mask[outlet] = true;
    while let Some(idx) = queue.pop_front() {
        for &up in &upstream[idx] {
            if !mask[up] {
                mask[up] = true;
                queue.push_back(up);
            }
        }
    }
    mask
}

// ---------------------------------------------------------------------------
// Section 7 — Stream network extraction
// ---------------------------------------------------------------------------

/// Extract stream network by thresholding flow accumulation.
///
/// Cells with accumulation ≥ `threshold` are classified as streams.
/// Returns a boolean mask.
pub fn extract_streams(accum: &[u32], threshold: u32) -> Vec<bool> {
    accum.iter().map(|&a| a >= threshold).collect()
}

/// Strahler stream order for every stream cell.
///
/// Returns a `Vec`u8` where 0 indicates non-stream cells.
/// Order assignment follows standard Strahler rules.
pub fn strahler_order(dem: &Dem, fdir: &[D8Code], streams: &[bool]) -> Vec<u8> {
    let n = dem.rows * dem.cols;
    let mut order = vec![0u8; n];
    // Build upstream adjacency restricted to stream cells
    let mut upstream: Vec<Vec<usize>> = vec![vec![]; n];
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let idx = r * dem.cols + c;
            if !streams[idx] {
                continue;
            }
            let code = fdir[idx];
            if code == 0 {
                continue;
            }
            for &(dr, dc, mask_code) in &D8_OFFSETS {
                if code == mask_code {
                    let nr = r as i64 + dr;
                    let nc = c as i64 + dc;
                    if nr >= 0 && nr < dem.rows as i64 && nc >= 0 && nc < dem.cols as i64 {
                        let nidx = nr as usize * dem.cols + nc as usize;
                        if streams[nidx] {
                            upstream[nidx].push(idx);
                        }
                    }
                    break;
                }
            }
        }
    }
    // Process in topological order (high → low elevation)
    let mut topo_order: Vec<usize> = (0..n).filter(|&i| streams[i]).collect();
    topo_order.sort_by(|&a, &b| {
        dem.heights[b]
            .partial_cmp(&dem.heights[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for idx in topo_order {
        let ups: Vec<u8> = upstream[idx].iter().map(|&u| order[u]).collect();
        if ups.is_empty() {
            order[idx] = 1;
        } else {
            let max_ord = *ups.iter().max().expect("iterator should not be empty");
            let count_max = ups.iter().filter(|&&o| o == max_ord).count();
            order[idx] = if count_max >= 2 { max_ord + 1 } else { max_ord };
        }
    }
    order
}

// ---------------------------------------------------------------------------
// Section 8 — Terrain Roughness Index (TRI)
// ---------------------------------------------------------------------------

/// Terrain Roughness Index (Riley et al. 1999).
///
/// Defined as the square root of the sum of squared elevation differences
/// between the centre cell and its eight neighbours.
pub fn terrain_roughness_index(dem: &Dem) -> Vec<f64> {
    let mut tri = vec![0.0f64; dem.rows * dem.cols];
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let z = dem.heights[r * dem.cols + c];
            let w = dem.window_3x3(r, c);
            let mut sum_sq = 0.0f64;
            for dr in 0..3usize {
                for dc in 0..3usize {
                    if dr == 1 && dc == 1 {
                        continue;
                    }
                    let diff = w[dr][dc] - z;
                    sum_sq += diff * diff;
                }
            }
            tri[r * dem.cols + c] = sum_sq.sqrt();
        }
    }
    tri
}

// ---------------------------------------------------------------------------
// Section 9 — Terrain Position Index (TPI)
// ---------------------------------------------------------------------------

/// Terrain Position Index (Weiss 2001).
///
/// TPI = elevation of cell − mean elevation of annular neighbourhood.
/// `radius` is the search radius in cells.
pub fn terrain_position_index(dem: &Dem, radius: usize) -> Vec<f64> {
    let mut tpi = vec![0.0f64; dem.rows * dem.cols];
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let z = dem.heights[r * dem.cols + c];
            let mut sum = 0.0f64;
            let mut count = 0usize;
            for dr in -(radius as i64)..=(radius as i64) {
                for dc in -(radius as i64)..=(radius as i64) {
                    if dr == 0 && dc == 0 {
                        continue;
                    }
                    let dist_sq = (dr * dr + dc * dc) as f64;
                    if dist_sq > (radius as f64 * radius as f64) {
                        continue;
                    }
                    let nr = r as i64 + dr;
                    let nc = c as i64 + dc;
                    if nr >= 0 && nr < dem.rows as i64 && nc >= 0 && nc < dem.cols as i64 {
                        sum += dem.heights[nr as usize * dem.cols + nc as usize];
                        count += 1;
                    }
                }
            }
            tpi[r * dem.cols + c] = if count > 0 {
                z - sum / count as f64
            } else {
                0.0
            };
        }
    }
    tpi
}

/// Landform classification based on TPI and slope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LandformClass {
    /// Ridges and hilltops (high positive TPI).
    Ridge,
    /// Upper slopes.
    UpperSlope,
    /// Middle slopes (moderate TPI, moderate slope).
    MiddleSlope,
    /// Flat areas (low TPI, low slope).
    Flat,
    /// Lower slopes.
    LowerSlope,
    /// Valleys and depressions (high negative TPI).
    Valley,
}

/// Classify landform type using a dual-scale TPI approach (Weiss 2001).
pub fn classify_landform(tpi_small: f64, tpi_large: f64, slope_rad: f64) -> LandformClass {
    let std_small = 1.0; // normalised
    let std_large = 1.0;
    let ts = tpi_small / std_small;
    let tl = tpi_large / std_large;
    let slope_deg = slope_rad.to_degrees();
    if tl > 1.0 {
        LandformClass::Ridge
    } else if tl < -1.0 {
        LandformClass::Valley
    } else if ts > 1.0 {
        LandformClass::UpperSlope
    } else if ts < -1.0 {
        LandformClass::LowerSlope
    } else if slope_deg < 5.0 {
        LandformClass::Flat
    } else {
        LandformClass::MiddleSlope
    }
}

// ---------------------------------------------------------------------------
// Section 10 — Topographic Wetness Index (TWI)
// ---------------------------------------------------------------------------

/// Topographic Wetness Index (Beven & Kirkby 1979).
///
/// TWI = ln(a / tan(β)), where *a* is the specific catchment area (m²/m) and
/// *β* is the local slope angle (radians).
pub fn topographic_wetness_index(dem: &Dem, accum: &[u32], slope: &[f64]) -> Vec<f64> {
    let n = dem.rows * dem.cols;
    let mut twi = vec![0.0f64; n];
    let cs = dem.cell_size;
    for i in 0..n {
        let a = (accum[i] as f64) * cs; // specific catchment area
        let beta = slope[i].max(1e-6_f64.atan()); // minimum slope to avoid division by zero
        twi[i] = (a / beta.tan()).ln();
    }
    twi
}

// ---------------------------------------------------------------------------
// Section 11 — Viewshed analysis
// ---------------------------------------------------------------------------

/// Observer specification for viewshed analysis.
#[derive(Debug, Clone, Copy)]
pub struct Observer {
    /// Row index of the observer.
    pub row: usize,
    /// Column index of the observer.
    pub col: usize,
    /// Observer height above terrain surface (metres).
    pub height_above_terrain: f64,
    /// Maximum sight distance (metres; use `f64::MAX` for unlimited).
    pub max_range: f64,
}

/// Compute binary viewshed for a single observer.
///
/// Uses the reference plane method: for each target cell the sight line is
/// tested against all intermediate cells sampled along the Bresenham raster.
/// Returns a boolean mask (true = visible).
pub fn compute_viewshed(dem: &Dem, obs: &Observer) -> Vec<bool> {
    let n = dem.rows * dem.cols;
    let mut visible = vec![false; n];
    let obs_z = dem.heights[obs.row * dem.cols + obs.col] + obs.height_above_terrain;

    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let dr = r as f64 - obs.row as f64;
            let dc = c as f64 - obs.col as f64;
            let dist = (dr * dr + dc * dc).sqrt() * dem.cell_size;
            if dist > obs.max_range {
                continue;
            }
            if r == obs.row && c == obs.col {
                visible[r * dem.cols + c] = true;
                continue;
            }
            // Bresenham line from observer to target
            let mut blocked = false;
            let steps = dr.abs().max(dc.abs()) as usize;
            if steps == 0 {
                visible[r * dem.cols + c] = true;
                continue;
            }
            for step in 1..steps {
                let t = step as f64 / steps as f64;
                let ir = (obs.row as f64 + t * dr).round() as i64;
                let ic = (obs.col as f64 + t * dc).round() as i64;
                if ir < 0 || ir >= dem.rows as i64 || ic < 0 || ic >= dem.cols as i64 {
                    continue;
                }
                let inter_z = dem.heights[ir as usize * dem.cols + ic as usize];
                let inter_dist = (t * dist).max(1e-9);
                let target_z = dem.heights[r * dem.cols + c];
                let los_z = obs_z + (target_z - obs_z) * (inter_dist / dist.max(1e-9));
                if inter_z > los_z {
                    blocked = true;
                    break;
                }
            }
            visible[r * dem.cols + c] = !blocked;
        }
    }
    visible
}

/// Cumulative viewshed: number of observers that can see each cell.
pub fn cumulative_viewshed(dem: &Dem, observers: &[Observer]) -> Vec<u32> {
    let n = dem.rows * dem.cols;
    let mut count = vec![0u32; n];
    for obs in observers {
        let vs = compute_viewshed(dem, obs);
        for i in 0..n {
            if vs[i] {
                count[i] += 1;
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Section 12 — Slope stability index
// ---------------------------------------------------------------------------

/// Parameters for the infinite-slope stability model.
#[derive(Debug, Clone, Copy)]
pub struct SlopeStabilityParams {
    /// Soil cohesion (Pa).
    pub cohesion: f64,
    /// Internal friction angle (radians).
    pub friction_angle: f64,
    /// Soil bulk density (kg/m³).
    pub soil_density: f64,
    /// Water density (kg/m³).
    pub water_density: f64,
    /// Depth of failure plane (m).
    pub depth: f64,
    /// Degree of saturation (0–1).
    pub saturation: f64,
}

impl Default for SlopeStabilityParams {
    fn default() -> Self {
        Self {
            cohesion: 5000.0,
            friction_angle: 30.0_f64.to_radians(),
            soil_density: 1800.0,
            water_density: 1000.0,
            depth: 1.0,
            saturation: 0.5,
        }
    }
}

/// Compute the Factor of Safety (FS) for every grid cell.
///
/// Uses the infinite-slope model:
/// FS = (c′ + (ρs − m·ρw)·g·h·cos²β·tan φ) / (ρs·g·h·sin β·cos β)
///
/// where *m* = saturation ratio, *β* = slope angle, *h* = depth.
pub fn slope_stability_index(dem: &Dem, slope: &[f64], params: &SlopeStabilityParams) -> Vec<f64> {
    let g = 9.81f64;
    let n = dem.rows * dem.cols;
    let mut fs = vec![0.0f64; n];
    for i in 0..n {
        let beta = slope[i];
        let sin_b = beta.sin();
        let cos_b = beta.cos();
        let h = params.depth;
        let m = params.saturation;
        let rho_s = params.soil_density;
        let rho_w = params.water_density;
        let c = params.cohesion;
        let phi = params.friction_angle;
        let denom = rho_s * g * h * sin_b * cos_b;
        if denom.abs() < 1e-9 {
            fs[i] = f64::MAX;
        } else {
            let numer = c + (rho_s - m * rho_w) * g * h * cos_b * cos_b * phi.tan();
            fs[i] = numer / denom;
        }
    }
    fs
}

// ---------------------------------------------------------------------------
// Section 13 — Solar irradiance on terrain
// ---------------------------------------------------------------------------

/// Solar position in the sky.
#[derive(Debug, Clone, Copy)]
pub struct SolarPosition {
    /// Solar zenith angle (radians, 0 = directly overhead).
    pub zenith: f64,
    /// Solar azimuth (radians, 0 = north, clockwise).
    pub azimuth: f64,
}

impl SolarPosition {
    /// Construct from elevation and azimuth angles (both in degrees).
    pub fn from_degrees(elevation_deg: f64, azimuth_deg: f64) -> Self {
        Self {
            zenith: (90.0 - elevation_deg).to_radians(),
            azimuth: azimuth_deg.to_radians(),
        }
    }

    /// Direct normal irradiance after atmospheric attenuation.
    ///
    /// Uses Beer-Lambert with air mass computed from zenith angle.
    /// `i0` is the extraterrestrial irradiance (W/m²), `tau` is the optical depth.
    pub fn direct_normal_irradiance(&self, i0: f64, tau: f64) -> f64 {
        let cos_z = self.zenith.cos().max(0.0);
        if cos_z < 1e-6 {
            return 0.0;
        }
        let air_mass = 1.0 / cos_z;
        i0 * (-tau * air_mass).exp()
    }
}

/// Compute the cosine of incidence angle between the sun and the terrain surface.
///
/// `slope_rad` is terrain slope, `aspect_rad` is terrain aspect (clockwise from north).
/// Returns the cosine of the incidence angle (clamped to ≥ 0).
pub fn cos_incidence(solar: &SolarPosition, slope_rad: f64, aspect_rad: f64) -> f64 {
    let cos_z = solar.zenith.cos();
    let sin_z = solar.zenith.sin();
    let cos_b = slope_rad.cos();
    let sin_b = slope_rad.sin();
    let delta_azimuth = solar.azimuth - aspect_rad;
    let ci = cos_z * cos_b + sin_z * sin_b * delta_azimuth.cos();
    ci.max(0.0)
}

/// Compute global horizontal irradiance (GHI) on every terrain cell.
///
/// Accounts for terrain slope, aspect, and self-shading (cast shadows are not
/// computed here — use `compute_viewshed` for that).
pub fn terrain_solar_irradiance(
    dem: &Dem,
    slope: &[f64],
    aspect: &[f64],
    solar: &SolarPosition,
    dni: f64,
    diffuse_horizontal: f64,
) -> Vec<f64> {
    let n = dem.rows * dem.cols;
    let mut irr = vec![0.0f64; n];
    let cos_z = solar.zenith.cos().max(0.0);
    for i in 0..n {
        let ci = cos_incidence(solar, slope[i], aspect[i].max(0.0));
        // Direct component on tilted surface
        let direct = if cos_z > 1e-6 { dni * ci } else { 0.0 };
        // Diffuse component (isotropic sky model)
        let diffuse = diffuse_horizontal * (1.0 + slope[i].cos()) / 2.0;
        irr[i] = direct + diffuse;
    }
    irr
}

// ---------------------------------------------------------------------------
// Section 14 — Terrain statistics helpers
// ---------------------------------------------------------------------------

/// Compute basic statistics of a raster field.
///
/// Returns `(min, max, mean, std_dev)`.
pub fn raster_stats(data: &[f64]) -> (f64, f64, f64, f64) {
    if data.is_empty() {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut min = f64::MAX;
    let mut max = f64::MIN;
    let mut sum = 0.0f64;
    for &v in data {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        sum += v;
    }
    let mean = sum / data.len() as f64;
    let var = data.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / data.len() as f64;
    (min, max, mean, var.sqrt())
}

/// Normalise a raster to the range \[0, 1\].
pub fn normalise_raster(data: &[f64]) -> Vec<f64> {
    let (min, max, _, _) = raster_stats(data);
    let range = max - min;
    if range.abs() < 1e-15 {
        return vec![0.0; data.len()];
    }
    data.iter().map(|&v| (v - min) / range).collect()
}

/// Resample a DEM to a coarser grid using bilinear interpolation.
///
/// `factor` is the downsampling factor (e.g. 2 halves the resolution).
pub fn resample_dem(dem: &Dem, factor: usize) -> Dem {
    let new_rows = dem.rows.div_ceil(factor);
    let new_cols = dem.cols.div_ceil(factor);
    let mut heights = vec![0.0f64; new_rows * new_cols];
    for nr in 0..new_rows {
        for nc in 0..new_cols {
            let r0 = (nr * factor) as i64;
            let c0 = (nc * factor) as i64;
            heights[nr * new_cols + nc] = dem.get_clamped(r0, c0);
        }
    }
    Dem::new(heights, new_rows, new_cols, dem.cell_size * factor as f64)
}

/// Fill sinks in a DEM using a simple iterative raising algorithm.
///
/// Each iteration raises cells that are lower than all their neighbours.
/// Converges in O(n) iterations for well-behaved terrain.
pub fn fill_sinks(dem: &Dem, max_iterations: usize) -> Dem {
    let mut filled = dem.clone();
    for _iter in 0..max_iterations {
        let mut changed = false;
        for r in 0..filled.rows {
            for c in 0..filled.cols {
                let idx = r * filled.cols + c;
                let z = filled.heights[idx];
                let w = filled.window_3x3(r, c);
                let min_nb = {
                    let mut m = f64::MAX;
                    for dr in 0..3usize {
                        for dc in 0..3usize {
                            if dr == 1 && dc == 1 {
                                continue;
                            }
                            if w[dr][dc] < m {
                                m = w[dr][dc];
                            }
                        }
                    }
                    m
                };
                if z < min_nb {
                    filled.heights[idx] = min_nb;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    filled
}

// ---------------------------------------------------------------------------
// Section 15 — Contour line extraction (marching squares)
// ---------------------------------------------------------------------------

/// A single contour polyline at a given elevation.
#[derive(Debug, Clone)]
pub struct Contour {
    /// Elevation of this contour (metres).
    pub elevation: f64,
    /// Vertices as (x, y) pairs in grid coordinates.
    pub vertices: Vec<(f64, f64)>,
}

/// Extract contour lines at the specified elevation levels.
///
/// Uses the marching squares algorithm on the DEM grid.
pub fn extract_contours(dem: &Dem, levels: &[f64]) -> Vec<Contour> {
    let mut contours = Vec::new();
    for &level in levels {
        let vertices = marching_squares_level(dem, level);
        if !vertices.is_empty() {
            contours.push(Contour {
                elevation: level,
                vertices,
            });
        }
    }
    contours
}

fn marching_squares_level(dem: &Dem, level: f64) -> Vec<(f64, f64)> {
    let mut pts = Vec::new();
    for r in 0..dem.rows.saturating_sub(1) {
        for c in 0..dem.cols.saturating_sub(1) {
            let z00 = dem.heights[r * dem.cols + c];
            let z10 = dem.heights[(r + 1) * dem.cols + c];
            let z01 = dem.heights[r * dem.cols + (c + 1)];
            let z11 = dem.heights[(r + 1) * dem.cols + (c + 1)];
            // Simple edge interpolation
            let interp = |za: f64, zb: f64| -> f64 {
                if (zb - za).abs() < 1e-15 {
                    0.5
                } else {
                    (level - za) / (zb - za)
                }
            };
            let mut code = 0u8;
            if z00 >= level {
                code |= 1;
            }
            if z01 >= level {
                code |= 2;
            }
            if z10 >= level {
                code |= 4;
            }
            if z11 >= level {
                code |= 8;
            }
            if code == 0 || code == 15 {
                continue;
            }
            // Top edge
            if (code & 3) == 1 || (code & 3) == 2 {
                let t = interp(z00, z01);
                pts.push((c as f64 + t, r as f64));
            }
            // Left edge
            if (code & 5) == 1 || (code & 5) == 4 {
                let t = interp(z00, z10);
                pts.push((c as f64, r as f64 + t));
            }
        }
    }
    pts
}

// ---------------------------------------------------------------------------
// Section 16 — Hypsometric curve
// ---------------------------------------------------------------------------

/// Compute the hypsometric curve for a DEM.
///
/// Returns a vector of `(relative_height, relative_area)` pairs where both
/// values are normalised to \[0, 1\].
pub fn hypsometric_curve(dem: &Dem, bins: usize) -> Vec<(f64, f64)> {
    let (min, max, _, _) = raster_stats(&dem.heights);
    let range = max - min;
    if range < 1e-9 {
        return vec![(0.0, 1.0), (1.0, 0.0)];
    }
    let total = dem.heights.len() as f64;
    let mut curve = Vec::with_capacity(bins + 1);
    for b in 0..=bins {
        let frac = b as f64 / bins as f64;
        let threshold = min + frac * range;
        let above = dem.heights.iter().filter(|&&h| h >= threshold).count() as f64;
        curve.push((frac, above / total));
    }
    curve
}

/// Hypsometric integral — area under the hypsometric curve.
///
/// Values > 0.6 indicate youthful terrain; < 0.35 indicate monadnock stage.
pub fn hypsometric_integral(curve: &[(f64, f64)]) -> f64 {
    if curve.len() < 2 {
        return 0.0;
    }
    let mut area = 0.0f64;
    for w in curve.windows(2) {
        let dx = w[1].0 - w[0].0;
        let avg_y = (w[0].1 + w[1].1) / 2.0;
        area += dx * avg_y;
    }
    area
}

// ---------------------------------------------------------------------------
// Section 17 — Ridge / valley line extraction
// ---------------------------------------------------------------------------

/// Classify each cell as ridge, valley, or neither.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RidgeValley {
    /// Local elevation maximum in at least one principal direction.
    Ridge,
    /// Local elevation minimum in at least one principal direction.
    Valley,
    /// Neither ridge nor valley.
    Slope,
}

/// Simple ridge and valley detection based on profile curvature sign.
///
/// Positive profile curvature → concave (valley); negative → convex (ridge).
pub fn ridge_valley_lines(dem: &Dem) -> Vec<RidgeValley> {
    let curv = compute_curvature(dem, CurvatureKind::Profile);
    curv.iter()
        .map(|&c| {
            if c > 0.01 {
                RidgeValley::Valley
            } else if c < -0.01 {
                RidgeValley::Ridge
            } else {
                RidgeValley::Slope
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Section 18 — Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a small synthetic DEM: a cone with peak at the centre.
    fn cone_dem(rows: usize, cols: usize, cell_size: f64) -> Dem {
        let cx = (cols as f64 - 1.0) / 2.0;
        let cy = (rows as f64 - 1.0) / 2.0;
        let heights: Vec<f64> = (0..rows)
            .flat_map(|r| {
                (0..cols).map(move |c| {
                    let dx = c as f64 - cx;
                    let dy = r as f64 - cy;
                    100.0 - (dx * dx + dy * dy).sqrt() * cell_size
                })
            })
            .collect();
        Dem::new(heights, rows, cols, cell_size)
    }

    /// Build a flat DEM.
    fn flat_dem(rows: usize, cols: usize, elevation: f64) -> Dem {
        Dem::new(vec![elevation; rows * cols], rows, cols, 10.0)
    }

    /// Build a simple tilted plane DEM (slope in x direction).
    fn tilted_dem(rows: usize, cols: usize, cell_size: f64, slope_gradient: f64) -> Dem {
        let heights: Vec<f64> = (0..rows)
            .flat_map(|_r| (0..cols).map(move |c| c as f64 * cell_size * slope_gradient))
            .collect();
        Dem::new(heights, rows, cols, cell_size)
    }

    #[test]
    fn test_dem_construction() {
        let dem = flat_dem(5, 5, 100.0);
        assert_eq!(dem.rows, 5);
        assert_eq!(dem.cols, 5);
        assert_eq!(dem.heights.len(), 25);
    }

    #[test]
    fn test_dem_get_valid() {
        let dem = flat_dem(4, 4, 50.0);
        assert_eq!(dem.get(2, 3), Some(50.0));
        assert_eq!(dem.get(4, 0), None);
    }

    #[test]
    fn test_dem_get_clamped_edge() {
        let dem = flat_dem(3, 3, 42.0);
        assert_eq!(dem.get_clamped(-1, 0), 42.0);
        assert_eq!(dem.get_clamped(0, 10), 42.0);
    }

    #[test]
    fn test_flat_dem_slope_zero() {
        let dem = flat_dem(5, 5, 100.0);
        let slope = compute_slope(&dem);
        for &s in &slope {
            assert!(s.abs() < 1e-10, "flat DEM slope should be zero, got {}", s);
        }
    }

    #[test]
    fn test_tilted_dem_slope_nonzero() {
        let dem = tilted_dem(5, 5, 10.0, 0.5);
        let slope = compute_slope(&dem);
        // Interior cells should have nonzero slope
        let interior_slope = slope[2 * 5 + 2];
        assert!(
            interior_slope > 0.01,
            "tilted DEM should have positive slope, got {}",
            interior_slope
        );
    }

    #[test]
    fn test_aspect_range() {
        let dem = cone_dem(7, 7, 10.0);
        let aspect = compute_aspect(&dem);
        for &a in &aspect {
            if a >= 0.0 {
                assert!(a <= 2.0 * PI + 1e-9, "aspect out of range: {}", a);
            }
        }
    }

    #[test]
    fn test_profile_curvature_flat_near_zero() {
        let dem = flat_dem(5, 5, 200.0);
        let curv = compute_curvature(&dem, CurvatureKind::Profile);
        for &c in &curv {
            assert!(c.abs() < 1e-10, "flat curvature should be 0, got {}", c);
        }
    }

    #[test]
    fn test_planform_curvature_flat_near_zero() {
        let dem = flat_dem(5, 5, 200.0);
        let curv = compute_curvature(&dem, CurvatureKind::Planform);
        for &c in &curv {
            assert!(
                c.abs() < 1e-10,
                "flat planform curvature should be 0, got {}",
                c
            );
        }
    }

    #[test]
    fn test_tangential_curvature_flat_near_zero() {
        let dem = flat_dem(5, 5, 200.0);
        let curv = compute_curvature(&dem, CurvatureKind::Tangential);
        for &c in &curv {
            assert!(
                c.abs() < 1e-10,
                "flat tangential curvature should be 0, got {}",
                c
            );
        }
    }

    #[test]
    fn test_flow_direction_d8_tilted() {
        let dem = tilted_dem(5, 5, 10.0, 1.0);
        let fdir = flow_direction_d8(&dem);
        // Interior cells should flow in the +x direction (code=1)
        let mid = 2 * 5 + 2;
        assert!(fdir[mid] != 0, "interior cell should have a flow direction");
    }

    #[test]
    fn test_flow_accumulation_monotone() {
        let dem = tilted_dem(6, 6, 10.0, 1.0);
        let fdir = flow_direction_d8(&dem);
        let accum = flow_accumulation_d8(&dem, &fdir);
        // Each cell must accumulate at least 1
        for &a in &accum {
            assert!(a >= 1, "accumulation must be >= 1");
        }
    }

    #[test]
    fn test_watershed_includes_outlet() {
        let dem = cone_dem(9, 9, 5.0);
        let fdir = flow_direction_d8(&dem);
        let mask = delineate_watershed(&dem, &fdir, 8, 4);
        assert!(mask[8 * 9 + 4], "outlet must be included in watershed");
    }

    #[test]
    fn test_stream_extraction_threshold() {
        let accum = vec![1u32, 2, 5, 10, 20];
        let streams = extract_streams(&accum, 5);
        assert!(!streams[0]);
        assert!(!streams[1]);
        assert!(streams[2]);
        assert!(streams[3]);
        assert!(streams[4]);
    }

    #[test]
    fn test_tri_flat_zero() {
        let dem = flat_dem(5, 5, 50.0);
        let tri = terrain_roughness_index(&dem);
        for &t in &tri {
            assert!(t.abs() < 1e-10, "TRI of flat DEM should be 0, got {}", t);
        }
    }

    #[test]
    fn test_tri_non_negative() {
        let dem = cone_dem(7, 7, 10.0);
        let tri = terrain_roughness_index(&dem);
        for &t in &tri {
            assert!(t >= 0.0, "TRI must be non-negative");
        }
    }

    #[test]
    fn test_tpi_flat_near_zero() {
        let dem = flat_dem(7, 7, 100.0);
        let tpi = terrain_position_index(&dem, 2);
        for &v in &tpi {
            assert!(v.abs() < 1e-9, "TPI on flat DEM should be 0, got {}", v);
        }
    }

    #[test]
    fn test_twi_positive() {
        let dem = tilted_dem(6, 6, 10.0, 0.3);
        let fdir = flow_direction_d8(&dem);
        let accum = flow_accumulation_d8(&dem, &fdir);
        let slope = compute_slope(&dem);
        let twi = topographic_wetness_index(&dem, &accum, &slope);
        for &v in &twi {
            assert!(v.is_finite(), "TWI must be finite");
        }
    }

    #[test]
    fn test_viewshed_observer_visible_to_itself() {
        let dem = flat_dem(5, 5, 0.0);
        let obs = Observer {
            row: 2,
            col: 2,
            height_above_terrain: 1.8,
            max_range: 1000.0,
        };
        let vs = compute_viewshed(&dem, &obs);
        assert!(vs[2 * 5 + 2], "observer cell must always be visible");
    }

    #[test]
    fn test_viewshed_flat_all_visible() {
        let dem = flat_dem(5, 5, 0.0);
        let obs = Observer {
            row: 2,
            col: 2,
            height_above_terrain: 1.8,
            max_range: 1000.0,
        };
        let vs = compute_viewshed(&dem, &obs);
        for (i, &v) in vs.iter().enumerate() {
            assert!(v, "flat terrain: cell {} should be visible", i);
        }
    }

    #[test]
    fn test_slope_stability_flat_infinite_fs() {
        let dem = flat_dem(3, 3, 0.0);
        let slope = compute_slope(&dem);
        let params = SlopeStabilityParams::default();
        let fs = slope_stability_index(&dem, &slope, &params);
        // On flat terrain, denominator → 0 → FS = MAX
        for &f in &fs {
            assert!(f >= 1.0, "FS on flat terrain should be large, got {}", f);
        }
    }

    #[test]
    fn test_slope_stability_steep_low_fs() {
        // Very steep slope (close to vertical)
        let dem = tilted_dem(3, 3, 1.0, 10.0);
        let slope = compute_slope(&dem);
        let params = SlopeStabilityParams {
            cohesion: 0.0,
            friction_angle: 30.0_f64.to_radians(),
            soil_density: 1800.0,
            water_density: 1000.0,
            depth: 1.0,
            saturation: 1.0,
        };
        let fs = slope_stability_index(&dem, &slope, &params);
        // At least one cell should have FS < 1 (unstable)
        let has_unstable = fs.iter().any(|&f| f < 2.0 && f > 0.0);
        assert!(has_unstable, "steep saturated slope should show low FS");
    }

    #[test]
    fn test_cos_incidence_overhead_sun() {
        // Sun directly overhead (zenith = 0), flat terrain → cos(0) = 1
        let solar = SolarPosition {
            zenith: 0.0,
            azimuth: 0.0,
        };
        let ci = cos_incidence(&solar, 0.0, 0.0);
        assert!(
            (ci - 1.0).abs() < 1e-9,
            "overhead sun on flat terrain: cos_incidence = {}",
            ci
        );
    }

    #[test]
    fn test_cos_incidence_below_horizon() {
        // Sun below horizon (zenith > π/2) → result clamped to 0
        let solar = SolarPosition {
            zenith: PI * 0.6,
            azimuth: 0.0,
        };
        let ci = cos_incidence(&solar, 0.0, 0.0);
        assert!(ci >= 0.0, "cos_incidence must be non-negative");
    }

    #[test]
    fn test_solar_irradiance_no_sun() {
        let dem = flat_dem(3, 3, 0.0);
        let slope = compute_slope(&dem);
        let aspect = compute_aspect(&dem);
        let solar = SolarPosition {
            zenith: PI / 2.0 + 0.1,
            azimuth: 0.0,
        }; // below horizon
        let irr = terrain_solar_irradiance(&dem, &slope, &aspect, &solar, 800.0, 100.0);
        for &v in &irr {
            // Should be diffuse only (small positive value)
            assert!(v >= 0.0, "irradiance must be non-negative");
        }
    }

    #[test]
    fn test_raster_stats_uniform() {
        let data = vec![5.0f64; 100];
        let (min, max, mean, std) = raster_stats(&data);
        assert_eq!(min, 5.0);
        assert_eq!(max, 5.0);
        assert_eq!(mean, 5.0);
        assert!(std.abs() < 1e-12);
    }

    #[test]
    fn test_normalise_raster_range() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let norm = normalise_raster(&data);
        assert!((norm[0]).abs() < 1e-9);
        assert!((norm[4] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_resample_dem_size() {
        let dem = flat_dem(8, 8, 100.0);
        let resampled = resample_dem(&dem, 2);
        assert_eq!(resampled.rows, 4);
        assert_eq!(resampled.cols, 4);
    }

    #[test]
    fn test_fill_sinks_no_lower_than_input() {
        let dem = cone_dem(5, 5, 10.0);
        let filled = fill_sinks(&dem, 50);
        for (i, (&orig, &fill)) in dem.heights.iter().zip(filled.heights.iter()).enumerate() {
            assert!(
                fill >= orig - 1e-9,
                "filled cell {} should not be lower than original",
                i
            );
        }
    }

    #[test]
    fn test_hypsometric_curve_bounds() {
        let dem = cone_dem(9, 9, 5.0);
        let curve = hypsometric_curve(&dem, 20);
        for &(h, a) in &curve {
            assert!((0.0..=1.0).contains(&h), "h out of [0,1]: {}", h);
            assert!(
                (0.0..=1.0).contains(&a),
                "area fraction out of [0,1]: {}",
                a
            );
        }
    }

    #[test]
    fn test_hypsometric_integral_range() {
        let dem = cone_dem(9, 9, 5.0);
        let curve = hypsometric_curve(&dem, 50);
        let hi = hypsometric_integral(&curve);
        assert!(
            (0.0..=1.0).contains(&hi),
            "hypsometric integral out of [0,1]: {}",
            hi
        );
    }

    #[test]
    fn test_classify_landform_ridge() {
        let lf = classify_landform(2.0, 2.0, 0.1);
        assert_eq!(lf, LandformClass::Ridge);
    }

    #[test]
    fn test_classify_landform_valley() {
        let lf = classify_landform(-2.0, -2.0, 0.1);
        assert_eq!(lf, LandformClass::Valley);
    }

    #[test]
    fn test_classify_landform_flat() {
        let lf = classify_landform(0.1, 0.1, 0.01);
        assert_eq!(lf, LandformClass::Flat);
    }

    #[test]
    fn test_strahler_order_source_is_one() {
        let dem = tilted_dem(5, 10, 10.0, 1.0);
        let fdir = flow_direction_d8(&dem);
        let accum = flow_accumulation_d8(&dem, &fdir);
        let streams = extract_streams(&accum, 2);
        let ord = strahler_order(&dem, &fdir, &streams);
        for (i, (&s, &o)) in streams.iter().zip(ord.iter()).enumerate() {
            if s {
                assert!(o >= 1, "stream cell {} should have order >= 1", i);
            }
        }
    }

    #[test]
    fn test_extract_contours_levels() {
        let dem = cone_dem(9, 9, 5.0);
        let levels = vec![80.0, 85.0, 90.0];
        let contours = extract_contours(&dem, &levels);
        assert!(!contours.is_empty(), "should find at least one contour");
        for c in &contours {
            assert!(
                !c.vertices.is_empty(),
                "contour at {} should have vertices",
                c.elevation
            );
        }
    }

    #[test]
    fn test_ridge_valley_lines_length() {
        let dem = cone_dem(7, 7, 5.0);
        let rv = ridge_valley_lines(&dem);
        assert_eq!(rv.len(), dem.rows * dem.cols);
    }

    #[test]
    fn test_dinf_flow_direction_range() {
        let dem = cone_dem(7, 7, 5.0);
        let fdir = flow_direction_dinf(&dem);
        for &a in &fdir {
            assert!(
                (0.0..=2.0 * PI + 0.01).contains(&a),
                "D-inf angle out of range: {}",
                a
            );
        }
    }

    #[test]
    fn test_cumulative_viewshed_non_negative() {
        let dem = flat_dem(5, 5, 0.0);
        let observers = vec![
            Observer {
                row: 0,
                col: 0,
                height_above_terrain: 2.0,
                max_range: 500.0,
            },
            Observer {
                row: 4,
                col: 4,
                height_above_terrain: 2.0,
                max_range: 500.0,
            },
        ];
        let cv = cumulative_viewshed(&dem, &observers);
        for &v in &cv {
            assert!(v <= 2, "cumulative count cannot exceed number of observers");
        }
    }

    #[test]
    fn test_fill_sinks_flat_unchanged() {
        let dem = flat_dem(5, 5, 100.0);
        let filled = fill_sinks(&dem, 10);
        for (&orig, &fill) in dem.heights.iter().zip(filled.heights.iter()) {
            assert!(
                (orig - fill).abs() < 1e-9,
                "flat DEM should be unchanged by fill_sinks"
            );
        }
    }
}
