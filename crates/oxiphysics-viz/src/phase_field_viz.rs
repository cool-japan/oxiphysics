#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase field visualization utilities.
//!
//! Provides data structures and coloring algorithms for phase-field simulation output:
//!
//! - [`OrderParameterField`] — 2-D order-parameter grid with coloring
//! - `InterfaceWidthViz` — visualization of the diffuse-interface width
//! - [`GrainStructure`] — grain-ID coloring with random-hue palette
//! - [`SolidificationFront`] — tracking and rendering of the solidification front
//! - [`DendriteTips`] — marking of dendrite tip positions
//! - [`ThermalPhaseOverlay`] — thermal field overlaid on a phase field
//! - [`CompositionFieldViz`] — composition field coloring
//! - `MisorientationViz` — misorientation angle display
//! - [`PhaseFractionPie`] — pie-chart data for phase fractions
//! - [`ConvergenceMonitor`] — phase-field convergence tracking

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use crate::colormap::{Colormap, map_scalar};
use crate::primitives::Color;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Map a normalised `t` in `[0,1]` through a given colormap.
#[inline]
fn cmap(t: f64, cm: Colormap) -> Color {
    map_scalar(t, 0.0, 1.0, cm)
}

/// Clamp `v` to `[lo, hi]`.
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.clamp(lo, hi)
}

/// Linear interpolation between `a` and `b` at `t`.
#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Convert a hue `h` in `[0,1]` to an RGB triple via HSV (S=1, V=1).
fn hue_to_rgb(h: f64) -> [f32; 3] {
    let h6 = h * 6.0;
    let i = h6.floor() as u32 % 6;
    let f = h6 - h6.floor();
    let q = 1.0 - f;
    let t = f;
    let (r, g, b) = match i {
        0 => (1.0, t, 0.0),
        1 => (q, 1.0, 0.0),
        2 => (0.0, 1.0, t),
        3 => (0.0, q, 1.0),
        4 => (t, 0.0, 1.0),
        _ => (1.0, 0.0, q),
    };
    [r as f32, g as f32, b as f32]
}

/// A simple deterministic hash of a `u32` to a float in `[0, 1)`.
fn hash_u32_to_f64(v: u32) -> f64 {
    // Xorshift-based integer hash
    let mut x = v.wrapping_add(0x9e3779b9);
    x = ((x >> 16) ^ x).wrapping_mul(0x45d9f3b);
    x = ((x >> 16) ^ x).wrapping_mul(0x45d9f3b);
    x = (x >> 16) ^ x;
    (x as f64) / (u32::MAX as f64 + 1.0)
}

// ---------------------------------------------------------------------------
// OrderParameterField
// ---------------------------------------------------------------------------

/// A 2-D order-parameter (phase field) grid, with values in `[−1, 1]` or `[0, 1]`.
///
/// The grid is stored row-major: `data[y * nx + x]`.
#[derive(Debug, Clone)]
pub struct OrderParameterField {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// Order-parameter values (typically in `[0, 1]`).
    pub data: Vec<f64>,
    /// Physical domain width.
    pub lx: f64,
    /// Physical domain height.
    pub ly: f64,
    /// Simulation time of this snapshot.
    pub time: f64,
}

impl OrderParameterField {
    /// Create a new field filled with `fill_value`.
    pub fn new(nx: usize, ny: usize, lx: f64, ly: f64, fill_value: f64) -> Self {
        Self {
            nx,
            ny,
            data: vec![fill_value; nx * ny],
            lx,
            ly,
            time: 0.0,
        }
    }

    /// Number of cells.
    pub fn len(&self) -> usize {
        self.nx * self.ny
    }

    /// Return `true` if the grid is empty.
    pub fn is_empty(&self) -> bool {
        self.nx == 0 || self.ny == 0
    }

    /// Read value at `(ix, iy)`.
    pub fn get(&self, ix: usize, iy: usize) -> f64 {
        self.data[iy * self.nx + ix]
    }

    /// Write value at `(ix, iy)`.
    pub fn set(&mut self, ix: usize, iy: usize, value: f64) {
        self.data[iy * self.nx + ix] = value;
    }

    /// Map each cell through a colormap.
    ///
    /// Values are normalised from `[min, max]` of the field (or `[0, 1]` if
    /// `auto_range` is `false`, in which case values are clamped to `[0, 1]`).
    pub fn to_rgba_buffer(&self, colormap: Colormap, auto_range: bool) -> Vec<u8> {
        let (vmin, vmax) = if auto_range {
            let mn = self.data.iter().cloned().fold(f64::INFINITY, f64::min);
            let mx = self.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let range = if (mx - mn).abs() < 1e-15 {
                1.0
            } else {
                mx - mn
            };
            (mn, mn + range)
        } else {
            (0.0, 1.0)
        };
        let mut buf = Vec::with_capacity(self.data.len() * 4);
        for &v in &self.data {
            let t = clamp((v - vmin) / (vmax - vmin), 0.0, 1.0);
            let c = cmap(t, colormap);
            buf.push((c.r * 255.0) as u8);
            buf.push((c.g * 255.0) as u8);
            buf.push((c.b * 255.0) as u8);
            buf.push(255u8);
        }
        buf
    }

    /// Volume fraction of solid (cells where φ ≥ `threshold`).
    pub fn solid_fraction(&self, threshold: f64) -> f64 {
        let n = self.data.iter().filter(|&&v| v >= threshold).count();
        n as f64 / self.data.len() as f64
    }

    /// Compute the gradient magnitude at `(ix, iy)` via central differences.
    pub fn gradient_magnitude(&self, ix: usize, iy: usize) -> f64 {
        let dx = self.lx / self.nx as f64;
        let dy = self.ly / self.ny as f64;
        let get = |i: isize, j: isize| {
            let ii = i.clamp(0, self.nx as isize - 1) as usize;
            let jj = j.clamp(0, self.ny as isize - 1) as usize;
            self.data[jj * self.nx + ii]
        };
        let gx =
            (get(ix as isize + 1, iy as isize) - get(ix as isize - 1, iy as isize)) / (2.0 * dx);
        let gy =
            (get(ix as isize, iy as isize + 1) - get(ix as isize, iy as isize - 1)) / (2.0 * dy);
        (gx * gx + gy * gy).sqrt()
    }

    /// Mean order parameter value.
    pub fn mean(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().sum::<f64>() / self.data.len() as f64
    }

    /// Standard deviation of the order parameter.
    pub fn std_dev(&self) -> f64 {
        if self.data.len() < 2 {
            return 0.0;
        }
        let m = self.mean();
        let var =
            self.data.iter().map(|&v| (v - m) * (v - m)).sum::<f64>() / self.data.len() as f64;
        var.sqrt()
    }
}

// ---------------------------------------------------------------------------
// InterfaceWidthViz
// ---------------------------------------------------------------------------

/// Result of interface-width analysis on an order-parameter field.
#[derive(Debug, Clone)]
pub struct InterfaceWidthResult {
    /// Estimated interface width in physical units (metres or dimensionless).
    pub width: f64,
    /// Number of interface cells detected (cells with `|∇φ|` above a threshold).
    pub interface_cells: usize,
    /// Mean gradient magnitude in interface cells.
    pub mean_gradient: f64,
    /// Maximum gradient magnitude found.
    pub max_gradient: f64,
    /// Colormap-encoded gradient image (RGBA, row-major).
    pub gradient_rgba: Vec<u8>,
}

/// Compute interface-width statistics and produce a gradient visualisation.
pub fn compute_interface_width(
    field: &OrderParameterField,
    gradient_threshold: f64,
    colormap: Colormap,
) -> InterfaceWidthResult {
    let n = field.nx * field.ny;
    let mut grads = vec![0.0f64; n];
    let mut max_g = 0.0f64;
    for iy in 0..field.ny {
        for ix in 0..field.nx {
            let g = field.gradient_magnitude(ix, iy);
            grads[iy * field.nx + ix] = g;
            if g > max_g {
                max_g = g;
            }
        }
    }
    let interface_cells = grads.iter().filter(|&&g| g >= gradient_threshold).count();
    let mean_gradient = if interface_cells > 0 {
        grads
            .iter()
            .filter(|&&g| g >= gradient_threshold)
            .sum::<f64>()
            / interface_cells as f64
    } else {
        0.0
    };
    // Approximate interface width from gradient profile:
    // W ≈ 1 / max_gradient (in normalised units)
    let width = if max_g > 1e-15 {
        1.0 / max_g
    } else {
        f64::INFINITY
    };

    let norm = if max_g > 1e-15 { max_g } else { 1.0 };
    let gradient_rgba: Vec<u8> = grads
        .iter()
        .flat_map(|&g| {
            let t = clamp(g / norm, 0.0, 1.0);
            let c = cmap(t, colormap);
            [
                (c.r * 255.0) as u8,
                (c.g * 255.0) as u8,
                (c.b * 255.0) as u8,
                255u8,
            ]
        })
        .collect();

    InterfaceWidthResult {
        width,
        interface_cells,
        mean_gradient,
        max_gradient: max_g,
        gradient_rgba,
    }
}

// ---------------------------------------------------------------------------
// GrainStructure
// ---------------------------------------------------------------------------

/// A 2-D grain-ID map.
///
/// Each cell stores a grain index.  Grain 0 is conventionally the matrix /
/// liquid phase.
#[derive(Debug, Clone)]
pub struct GrainStructure {
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
    /// Grain IDs, row-major (`grain_id[y * nx + x]`).
    pub grain_id: Vec<u32>,
    /// Maximum grain ID present.
    pub max_grain_id: u32,
}

impl GrainStructure {
    /// Construct from a flat grain-ID array.
    ///
    /// Panics if `grain_id.len() != nx * ny`.
    pub fn new(nx: usize, ny: usize, grain_id: Vec<u32>) -> Self {
        assert_eq!(grain_id.len(), nx * ny);
        let max_id = grain_id.iter().cloned().max().unwrap_or(0);
        Self {
            nx,
            ny,
            grain_id,
            max_grain_id: max_id,
        }
    }

    /// Assign a unique random hue to each grain and render to RGBA.
    ///
    /// Grain 0 is rendered black (background / liquid).
    pub fn to_rgba_by_grain_id(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.nx * self.ny * 4);
        for &id in &self.grain_id {
            if id == 0 {
                buf.extend_from_slice(&[0u8, 0, 0, 255]);
            } else {
                let h = hash_u32_to_f64(id);
                let [r, g, b] = hue_to_rgb(h);
                buf.push((r * 255.0) as u8);
                buf.push((g * 255.0) as u8);
                buf.push((b * 255.0) as u8);
                buf.push(255u8);
            }
        }
        buf
    }

    /// Count of unique grain IDs (excluding grain 0).
    pub fn grain_count(&self) -> usize {
        let mut ids: Vec<u32> = self.grain_id.iter().cloned().filter(|&id| id > 0).collect();
        ids.sort_unstable();
        ids.dedup();
        ids.len()
    }

    /// Area fraction of a specific grain (cells with `id / total cells`).
    pub fn area_fraction(&self, id: u32) -> f64 {
        let count = self.grain_id.iter().filter(|&&g| g == id).count();
        count as f64 / self.grain_id.len() as f64
    }

    /// Largest grain ID by cell count.
    pub fn dominant_grain(&self) -> u32 {
        let mut counts: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for &id in &self.grain_id {
            if id > 0 {
                *counts.entry(id).or_insert(0) += 1;
            }
        }
        counts
            .into_iter()
            .max_by_key(|(_, c)| *c)
            .map(|(id, _)| id)
            .unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// SolidificationFront
// ---------------------------------------------------------------------------

/// A detected solidification front: a contour at φ = `iso_level`.
#[derive(Debug, Clone)]
pub struct SolidificationFront {
    /// Isovalue used to extract the front.
    pub iso_level: f64,
    /// Positions of the front contour cells in `(ix, iy)` grid coordinates.
    pub contour_cells: Vec<[usize; 2]>,
    /// Mean position of the front in physical coordinates `[x, y]`.
    pub mean_position: [f64; 2],
    /// Maximum x-coordinate of the front (approximation of front advance).
    pub max_x: f64,
    /// Simulation time at which the front was detected.
    pub time: f64,
}

/// Extract the solidification front at `iso_level` from an order-parameter field.
///
/// A cell belongs to the front if its order parameter crosses `iso_level`
/// relative to at least one of its four neighbours.
pub fn extract_solidification_front(
    field: &OrderParameterField,
    iso_level: f64,
) -> SolidificationFront {
    let mut contour_cells = Vec::new();
    let dx = field.lx / field.nx as f64;
    let dy = field.ly / field.ny as f64;

    for iy in 0..field.ny {
        for ix in 0..field.nx {
            let phi = field.get(ix, iy);
            // Check 4-connected neighbours
            let on_front = [[1i32, 0], [-1, 0], [0, 1], [0, -1]]
                .iter()
                .any(|&[dx2, dy2]| {
                    let nx2 = ix as i32 + dx2;
                    let ny2 = iy as i32 + dy2;
                    if nx2 < 0 || ny2 < 0 || nx2 >= field.nx as i32 || ny2 >= field.ny as i32 {
                        return false;
                    }
                    let phi_n = field.get(nx2 as usize, ny2 as usize);
                    (phi - iso_level) * (phi_n - iso_level) < 0.0
                });
            if on_front {
                contour_cells.push([ix, iy]);
            }
        }
    }

    let (mut sum_x, mut sum_y) = (0.0f64, 0.0f64);
    let mut max_x = f64::NEG_INFINITY;
    for &[ix, iy] in &contour_cells {
        let x = (ix as f64 + 0.5) * dx;
        let y = (iy as f64 + 0.5) * dy;
        sum_x += x;
        sum_y += y;
        if x > max_x {
            max_x = x;
        }
    }
    let n = contour_cells.len().max(1) as f64;
    SolidificationFront {
        iso_level,
        contour_cells,
        mean_position: [sum_x / n, sum_y / n],
        max_x: if max_x.is_finite() { max_x } else { 0.0 },
        time: field.time,
    }
}

/// Render the solidification front as an RGBA overlay on a base image.
///
/// `base_rgba` must have length `nx * ny * 4`.  Front cells are painted with `color`.
pub fn render_front_overlay(
    base_rgba: &[u8],
    nx: usize,
    _ny: usize,
    front: &SolidificationFront,
    color: Color,
) -> Vec<u8> {
    let mut buf = base_rgba.to_vec();
    for &[ix, iy] in &front.contour_cells {
        let idx = (iy * nx + ix) * 4;
        if idx + 3 < buf.len() {
            buf[idx] = (color.r * 255.0) as u8;
            buf[idx + 1] = (color.g * 255.0) as u8;
            buf[idx + 2] = (color.b * 255.0) as u8;
            buf[idx + 3] = 255;
        }
    }
    buf
}

// ---------------------------------------------------------------------------
// DendriteTips
// ---------------------------------------------------------------------------

/// A marked dendrite tip in the simulation domain.
#[derive(Debug, Clone, Copy)]
pub struct DendriteTip {
    /// Grid coordinates `(ix, iy)`.
    pub grid_pos: [usize; 2],
    /// Physical coordinates `[x, y]`.
    pub phys_pos: [f64; 2],
    /// Local order-parameter value at the tip.
    pub phi: f64,
    /// Local gradient magnitude at the tip.
    pub grad_mag: f64,
    /// Dendrite arm index (for bookkeeping).
    pub arm_id: u32,
}

/// Collection of detected dendrite tips.
#[derive(Debug, Default, Clone)]
pub struct DendriteTips {
    /// Detected tips.
    pub tips: Vec<DendriteTip>,
}

impl DendriteTips {
    /// Create an empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Detect dendrite tips as local maxima of gradient magnitude in the interface region.
    ///
    /// A cell qualifies as a tip if:
    /// 1. Its gradient magnitude is at least `grad_threshold`.
    /// 2. It has the highest gradient in a `radius`-cell neighbourhood.
    pub fn detect(field: &OrderParameterField, grad_threshold: f64, radius: usize) -> Self {
        let dx = field.lx / field.nx as f64;
        let dy = field.ly / field.ny as f64;

        let grads: Vec<f64> = (0..field.ny)
            .flat_map(|iy| (0..field.nx).map(move |ix| (ix, iy)))
            .map(|(ix, iy)| field.gradient_magnitude(ix, iy))
            .collect();

        let mut tips = Vec::new();
        let mut arm_id = 0u32;

        for iy in 0..field.ny {
            for ix in 0..field.nx {
                let g = grads[iy * field.nx + ix];
                if g < grad_threshold {
                    continue;
                }
                // Check local maximum
                let is_max = (ix.saturating_sub(radius)..=(ix + radius).min(field.nx - 1))
                    .flat_map(|nx2| {
                        let iy_lo = iy.saturating_sub(radius);
                        let iy_hi = (iy + radius).min(field.ny - 1);
                        (iy_lo..=iy_hi).map(move |ny2| (nx2, ny2))
                    })
                    .all(|(nx2, ny2)| grads[ny2 * field.nx + nx2] <= g);
                if is_max {
                    arm_id += 1;
                    tips.push(DendriteTip {
                        grid_pos: [ix, iy],
                        phys_pos: [(ix as f64 + 0.5) * dx, (iy as f64 + 0.5) * dy],
                        phi: field.get(ix, iy),
                        grad_mag: g,
                        arm_id,
                    });
                }
            }
        }
        Self { tips }
    }

    /// Number of detected tips.
    pub fn len(&self) -> usize {
        self.tips.len()
    }

    /// Return `true` if no tips were detected.
    pub fn is_empty(&self) -> bool {
        self.tips.is_empty()
    }

    /// Render tips as marker pixels on an RGBA buffer.
    ///
    /// Each tip is drawn as a `marker_radius`-pixel filled circle.
    pub fn render_markers(
        &self,
        base_rgba: &[u8],
        nx: usize,
        ny: usize,
        marker_radius: usize,
        color: Color,
    ) -> Vec<u8> {
        let mut buf = base_rgba.to_vec();
        let r = marker_radius as isize;
        for tip in &self.tips {
            let [cx, cy] = tip.grid_pos;
            for dy2 in -r..=r {
                for dx2 in -r..=r {
                    if dx2 * dx2 + dy2 * dy2 > r * r {
                        continue;
                    }
                    let px = cx as isize + dx2;
                    let py = cy as isize + dy2;
                    if px < 0 || py < 0 || px >= nx as isize || py >= ny as isize {
                        continue;
                    }
                    let idx = (py as usize * nx + px as usize) * 4;
                    if idx + 3 < buf.len() {
                        buf[idx] = (color.r * 255.0) as u8;
                        buf[idx + 1] = (color.g * 255.0) as u8;
                        buf[idx + 2] = (color.b * 255.0) as u8;
                        buf[idx + 3] = 255;
                    }
                }
            }
        }
        buf
    }
}

// ---------------------------------------------------------------------------
// ThermalPhaseOverlay
// ---------------------------------------------------------------------------

/// A blended visualisation of a thermal field overlaid on an order-parameter field.
///
/// The thermal field is alpha-composited on top of the phase field image.
#[derive(Debug, Clone)]
pub struct ThermalPhaseOverlay {
    /// Width of the domain in cells.
    pub nx: usize,
    /// Height of the domain in cells.
    pub ny: usize,
    /// Blending weight for the thermal field `[0, 1]`; phase weight = `1 - alpha`.
    pub thermal_alpha: f32,
    /// Colormap used for the phase field.
    pub phase_colormap: Colormap,
    /// Colormap used for the thermal field.
    pub thermal_colormap: Colormap,
}

impl ThermalPhaseOverlay {
    /// Construct a new overlay.
    pub fn new(
        nx: usize,
        ny: usize,
        thermal_alpha: f32,
        phase_colormap: Colormap,
        thermal_colormap: Colormap,
    ) -> Self {
        Self {
            nx,
            ny,
            thermal_alpha: thermal_alpha.clamp(0.0, 1.0),
            phase_colormap,
            thermal_colormap,
        }
    }

    /// Blend the phase and thermal fields into a single RGBA buffer.
    ///
    /// Both `phase_field.data` and `thermal_field` must have `nx * ny` elements.
    /// The thermal field is normalised using its min/max.
    pub fn render(&self, phase_field: &OrderParameterField, thermal_field: &[f64]) -> Vec<u8> {
        assert_eq!(phase_field.data.len(), self.nx * self.ny);
        assert_eq!(thermal_field.len(), self.nx * self.ny);

        let t_min = thermal_field.iter().cloned().fold(f64::INFINITY, f64::min);
        let t_max = thermal_field
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let t_range = if (t_max - t_min).abs() < 1e-15 {
            1.0
        } else {
            t_max - t_min
        };

        let p_min = phase_field
            .data
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let p_max = phase_field
            .data
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let p_range = if (p_max - p_min).abs() < 1e-15 {
            1.0
        } else {
            p_max - p_min
        };

        let alpha = self.thermal_alpha as f64;
        let mut buf = Vec::with_capacity(self.nx * self.ny * 4);

        for i in 0..self.nx * self.ny {
            let tp = clamp((phase_field.data[i] - p_min) / p_range, 0.0, 1.0);
            let tt = clamp((thermal_field[i] - t_min) / t_range, 0.0, 1.0);
            let cp = cmap(tp, self.phase_colormap);
            let ct = cmap(tt, self.thermal_colormap);
            let r = lerp(cp.r as f64, ct.r as f64, alpha);
            let g = lerp(cp.g as f64, ct.g as f64, alpha);
            let b = lerp(cp.b as f64, ct.b as f64, alpha);
            buf.push((r * 255.0) as u8);
            buf.push((g * 255.0) as u8);
            buf.push((b * 255.0) as u8);
            buf.push(255u8);
        }
        buf
    }
}

// ---------------------------------------------------------------------------
// CompositionFieldViz
// ---------------------------------------------------------------------------

/// Visualiser for a composition (solute concentration) field.
#[derive(Debug, Clone)]
pub struct CompositionFieldViz {
    /// Width of the domain in cells.
    pub nx: usize,
    /// Height of the domain in cells.
    pub ny: usize,
    /// Physical minimum composition value.
    pub c_min: f64,
    /// Physical maximum composition value.
    pub c_max: f64,
    /// Colormap used for the composition field.
    pub colormap: Colormap,
}

impl CompositionFieldViz {
    /// Construct a new visualiser.
    pub fn new(nx: usize, ny: usize, c_min: f64, c_max: f64, colormap: Colormap) -> Self {
        Self {
            nx,
            ny,
            c_min,
            c_max,
            colormap,
        }
    }

    /// Render the composition field `data` to RGBA.
    ///
    /// Values outside `[c_min, c_max]` are clamped.
    pub fn render(&self, data: &[f64]) -> Vec<u8> {
        let range = if (self.c_max - self.c_min).abs() < 1e-15 {
            1.0
        } else {
            self.c_max - self.c_min
        };
        let mut buf = Vec::with_capacity(data.len() * 4);
        for &v in data {
            let t = clamp((v - self.c_min) / range, 0.0, 1.0);
            let c = cmap(t, self.colormap);
            buf.push((c.r * 255.0) as u8);
            buf.push((c.g * 255.0) as u8);
            buf.push((c.b * 255.0) as u8);
            buf.push(255u8);
        }
        buf
    }

    /// Mean composition (average of `data`).
    pub fn mean_composition(data: &[f64]) -> f64 {
        if data.is_empty() {
            return 0.0;
        }
        data.iter().sum::<f64>() / data.len() as f64
    }

    /// Histogram of composition values.
    ///
    /// Returns `n_bins` bucket counts normalised to `[0, n]`.
    pub fn histogram(data: &[f64], c_min: f64, c_max: f64, n_bins: usize) -> Vec<usize> {
        if n_bins == 0 || data.is_empty() {
            return vec![0; n_bins.max(1)];
        }
        let range = if (c_max - c_min).abs() < 1e-15 {
            1.0
        } else {
            c_max - c_min
        };
        let mut bins = vec![0usize; n_bins];
        for &v in data {
            let t = clamp((v - c_min) / range, 0.0, 1.0);
            let bin = (t * n_bins as f64) as usize;
            let bin = bin.min(n_bins - 1);
            bins[bin] += 1;
        }
        bins
    }
}

// ---------------------------------------------------------------------------
// MisorientationViz
// ---------------------------------------------------------------------------

/// Misorientation angle between grains, stored per grain pair.
#[derive(Debug, Clone, Copy)]
pub struct GrainPairMisorientation {
    /// First grain ID.
    pub grain_a: u32,
    /// Second grain ID.
    pub grain_b: u32,
    /// Misorientation angle in radians.
    pub angle_rad: f64,
}

impl GrainPairMisorientation {
    /// Construct a new pair.
    pub fn new(grain_a: u32, grain_b: u32, angle_rad: f64) -> Self {
        Self {
            grain_a,
            grain_b,
            angle_rad,
        }
    }

    /// Return the angle in degrees.
    pub fn angle_deg(&self) -> f64 {
        self.angle_rad.to_degrees()
    }
}

/// Render a misorientation map: each boundary cell is coloured by the local
/// misorientation angle.
///
/// A boundary cell is one with a neighbour in a different grain.
/// Non-boundary cells are painted with `background`.
pub fn render_misorientation_map(
    grain_structure: &GrainStructure,
    misorientation_table: &[GrainPairMisorientation],
    max_angle_deg: f64,
    background: Color,
    colormap: Colormap,
) -> Vec<u8> {
    let nx = grain_structure.nx;
    let ny = grain_structure.ny;
    let n = nx * ny;
    let mut buf = Vec::with_capacity(n * 4);

    // Build lookup: (min_id, max_id) → angle
    let mut lookup: std::collections::HashMap<(u32, u32), f64> = std::collections::HashMap::new();
    for m in misorientation_table {
        let key = (m.grain_a.min(m.grain_b), m.grain_a.max(m.grain_b));
        lookup.insert(key, m.angle_rad.to_degrees());
    }

    for iy in 0..ny {
        for ix in 0..nx {
            let id = grain_structure.grain_id[iy * nx + ix];
            // Find max misorientation with 4-connected neighbours
            let mut max_mis = 0.0f64;
            let mut is_boundary = false;
            for &[ddx, ddy] in &[[1i32, 0], [-1, 0], [0, 1i32], [0, -1]] {
                let nx2 = ix as i32 + ddx;
                let ny2 = iy as i32 + ddy;
                if nx2 < 0 || ny2 < 0 || nx2 >= nx as i32 || ny2 >= ny as i32 {
                    continue;
                }
                let nb_id = grain_structure.grain_id[ny2 as usize * nx + nx2 as usize];
                if nb_id != id {
                    is_boundary = true;
                    let key = (id.min(nb_id), id.max(nb_id));
                    let angle = lookup.get(&key).cloned().unwrap_or(0.0);
                    if angle > max_mis {
                        max_mis = angle;
                    }
                }
            }
            if is_boundary {
                let t = clamp(max_mis / max_angle_deg.max(1.0), 0.0, 1.0);
                let c = cmap(t, colormap);
                buf.push((c.r * 255.0) as u8);
                buf.push((c.g * 255.0) as u8);
                buf.push((c.b * 255.0) as u8);
                buf.push(255u8);
            } else {
                buf.push((background.r * 255.0) as u8);
                buf.push((background.g * 255.0) as u8);
                buf.push((background.b * 255.0) as u8);
                buf.push((background.a * 255.0) as u8);
            }
        }
    }
    buf
}

// ---------------------------------------------------------------------------
// PhaseFractionPie
// ---------------------------------------------------------------------------

/// Data for a pie chart showing phase fractions.
#[derive(Debug, Clone)]
pub struct PhaseFractionPie {
    /// Phase names.
    pub names: Vec<String>,
    /// Volume fractions (must sum to 1 after normalisation).
    pub fractions: Vec<f64>,
    /// Colors assigned to each phase.
    pub colors: Vec<Color>,
}

impl PhaseFractionPie {
    /// Compute phase fractions from an order-parameter field using user-defined
    /// threshold bands.
    ///
    /// `bands`: list of `(name, lo, hi, color)` — a cell belongs to a phase if
    /// its value lies in `[lo, hi)`.  Overlapping bands are allowed; a cell
    /// goes to the first matching band.
    pub fn from_order_parameter(
        field: &OrderParameterField,
        bands: &[(&str, f64, f64, Color)],
    ) -> Self {
        let n = field.data.len();
        let mut counts = vec![0usize; bands.len()];
        for &v in &field.data {
            for (k, &(_, lo, hi, _)) in bands.iter().enumerate() {
                if v >= lo && v < hi {
                    counts[k] += 1;
                    break;
                }
            }
        }
        let n_f = n as f64;
        let fractions: Vec<f64> = counts.iter().map(|&c| c as f64 / n_f).collect();
        let names: Vec<String> = bands.iter().map(|(s, _, _, _)| s.to_string()).collect();
        let colors: Vec<Color> = bands.iter().map(|(_, _, _, c)| *c).collect();
        Self {
            names,
            fractions,
            colors,
        }
    }

    /// Generate pie-segment arc angles `[(start_rad, end_rad, color)]`.
    pub fn arc_segments(&self) -> Vec<(f64, f64, Color)> {
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut segments = Vec::with_capacity(self.fractions.len());
        let mut start = 0.0f64;
        for (i, &frac) in self.fractions.iter().enumerate() {
            let end = start + frac * two_pi;
            segments.push((start, end, self.colors[i]));
            start = end;
        }
        segments
    }

    /// Return the index of the dominant (largest fraction) phase.
    pub fn dominant_phase_index(&self) -> Option<usize> {
        self.fractions
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }

    /// Total fraction covered by the defined bands.
    pub fn total_fraction(&self) -> f64 {
        self.fractions.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// ConvergenceMonitor
// ---------------------------------------------------------------------------

/// A single convergence sample in a phase-field iteration.
#[derive(Debug, Clone, Copy)]
pub struct ConvergenceSample {
    /// Iteration number.
    pub iteration: u64,
    /// Simulation time.
    pub time: f64,
    /// Residual (e.g. max ‖δφ‖).
    pub residual: f64,
    /// Free energy (optional, `f64::NAN` if not available).
    pub free_energy: f64,
    /// Solid fraction at this step.
    pub solid_fraction: f64,
}

/// Monitor and visualise the convergence of a phase-field simulation.
#[derive(Debug, Default, Clone)]
pub struct ConvergenceMonitor {
    /// Recorded samples.
    pub samples: Vec<ConvergenceSample>,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Whether convergence was declared.
    pub converged: bool,
}

impl ConvergenceMonitor {
    /// Create a new monitor with the given tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            samples: Vec::new(),
            tolerance,
            converged: false,
        }
    }

    /// Record a new sample.  Marks `converged = true` if `residual < tolerance`.
    pub fn record(
        &mut self,
        iteration: u64,
        time: f64,
        residual: f64,
        free_energy: f64,
        solid_fraction: f64,
    ) {
        if residual < self.tolerance {
            self.converged = true;
        }
        self.samples.push(ConvergenceSample {
            iteration,
            time,
            residual,
            free_energy,
            solid_fraction,
        });
    }

    /// Number of recorded samples.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Return `true` if no samples have been recorded.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Latest residual, or `f64::NAN` if empty.
    pub fn latest_residual(&self) -> f64 {
        self.samples.last().map(|s| s.residual).unwrap_or(f64::NAN)
    }

    /// Minimum residual across all samples.
    pub fn min_residual(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.residual)
            .fold(f64::INFINITY, f64::min)
    }

    /// Rate of convergence: slope of log(residual) vs iteration (least-squares).
    ///
    /// Returns `f64::NAN` for fewer than 2 samples.
    pub fn convergence_rate(&self) -> f64 {
        let n = self.samples.len();
        if n < 2 {
            return f64::NAN;
        }
        let xs: Vec<f64> = self.samples.iter().map(|s| s.iteration as f64).collect();
        let ys: Vec<f64> = self
            .samples
            .iter()
            .map(|s| {
                if s.residual > 0.0 {
                    s.residual.ln()
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect();
        let finite_pairs: Vec<(f64, f64)> = xs
            .iter()
            .zip(ys.iter())
            .filter(|(_, y)| y.is_finite())
            .map(|(&x, &y)| (x, y))
            .collect();
        let m = finite_pairs.len();
        if m < 2 {
            return f64::NAN;
        }
        let sum_x: f64 = finite_pairs.iter().map(|(x, _)| x).sum();
        let sum_y: f64 = finite_pairs.iter().map(|(_, y)| y).sum();
        let sum_xx: f64 = finite_pairs.iter().map(|(x, _)| x * x).sum();
        let sum_xy: f64 = finite_pairs.iter().map(|(x, y)| x * y).sum();
        let mf = m as f64;
        let denom = mf * sum_xx - sum_x * sum_x;
        if denom.abs() < 1e-15 {
            return f64::NAN;
        }
        (mf * sum_xy - sum_x * sum_y) / denom
    }

    /// Render the residual history to an RGBA plot buffer of size `width × height`.
    ///
    /// Uses a log-scale y-axis.  Returns an empty vec if no samples are present.
    pub fn render_residual_plot(
        &self,
        width: usize,
        height: usize,
        bg: Color,
        line_color: Color,
    ) -> Vec<u8> {
        if self.samples.is_empty() || width == 0 || height == 0 {
            return Vec::new();
        }
        let n = self.samples.len();
        // Fill background
        let mut buf = [
            (bg.r * 255.0) as u8,
            (bg.g * 255.0) as u8,
            (bg.b * 255.0) as u8,
            255u8,
        ]
        .repeat(width * height);

        let log_res: Vec<f64> = self
            .samples
            .iter()
            .map(|s| {
                if s.residual > 0.0 {
                    s.residual.ln()
                } else {
                    -50.0
                }
            })
            .collect();
        let log_min = log_res.iter().cloned().fold(f64::INFINITY, f64::min);
        let log_max = log_res.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let log_range = if (log_max - log_min).abs() < 1e-12 {
            1.0
        } else {
            log_max - log_min
        };

        let lc_r = (line_color.r * 255.0) as u8;
        let lc_g = (line_color.g * 255.0) as u8;
        let lc_b = (line_color.b * 255.0) as u8;

        for i in 1..n {
            let x0 = ((i - 1) * (width - 1)) / (n - 1);
            let x1 = (i * (width - 1)) / (n - 1);
            let y0_t = clamp((log_res[i - 1] - log_min) / log_range, 0.0, 1.0);
            let y1_t = clamp((log_res[i] - log_min) / log_range, 0.0, 1.0);
            // Invert: high residual at top → low y pixel
            let y0 = ((1.0 - y0_t) * (height - 1) as f64) as usize;
            let y1 = ((1.0 - y1_t) * (height - 1) as f64) as usize;

            // Draw line segment using Bresenham
            let steps = (x1 as isize - x0 as isize)
                .abs()
                .max((y1 as isize - y0 as isize).abs())
                .max(1);
            for step in 0..=steps {
                let t = step as f64 / steps as f64;
                let px = (lerp(x0 as f64, x1 as f64, t)) as usize;
                let py = (lerp(y0 as f64, y1 as f64, t)) as usize;
                if px < width && py < height {
                    let idx = (py * width + px) * 4;
                    buf[idx] = lc_r;
                    buf[idx + 1] = lc_g;
                    buf[idx + 2] = lc_b;
                    buf[idx + 3] = 255;
                }
            }
        }
        buf
    }

    /// Export convergence history to CSV.
    pub fn to_csv(&self) -> String {
        let mut lines = vec!["iteration,time,residual,free_energy,solid_fraction".to_string()];
        for s in &self.samples {
            lines.push(format!(
                "{},{:.8},{:.6e},{:.8},{:.6}",
                s.iteration, s.time, s.residual, s.free_energy, s.solid_fraction
            ));
        }
        lines.join("\n")
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- helpers ----

    fn make_field(nx: usize, ny: usize) -> OrderParameterField {
        let mut f = OrderParameterField::new(nx, ny, 1.0, 1.0, 0.0);
        // Fill left half with 1.0 (solid), right half with 0.0 (liquid)
        for iy in 0..ny {
            for ix in 0..nx / 2 {
                f.set(ix, iy, 1.0);
            }
        }
        f
    }

    // ---- OrderParameterField ----

    #[test]
    fn test_order_param_field_solid_fraction() {
        let f = make_field(10, 10);
        let sf = f.solid_fraction(0.5);
        assert!((sf - 0.5).abs() < 1e-9, "expected 0.5, got {sf}");
    }

    #[test]
    fn test_order_param_field_mean() {
        let f = make_field(10, 10);
        assert!((f.mean() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_order_param_field_std_dev_uniform() {
        let f = OrderParameterField::new(4, 4, 1.0, 1.0, 0.5);
        assert!(f.std_dev().abs() < 1e-12);
    }

    #[test]
    fn test_order_param_field_to_rgba_length() {
        let f = make_field(8, 8);
        let rgba = f.to_rgba_buffer(Colormap::Viridis, true);
        assert_eq!(rgba.len(), 8 * 8 * 4);
    }

    #[test]
    fn test_order_param_field_gradient_zero_uniform() {
        let f = OrderParameterField::new(5, 5, 1.0, 1.0, 0.7);
        let g = f.gradient_magnitude(2, 2);
        assert!(
            g.abs() < 1e-10,
            "uniform field gradient should be 0, got {g}"
        );
    }

    #[test]
    fn test_order_param_field_get_set() {
        let mut f = OrderParameterField::new(4, 4, 1.0, 1.0, 0.0);
        f.set(2, 3, 0.75);
        assert!((f.get(2, 3) - 0.75).abs() < 1e-12);
    }

    // ---- InterfaceWidthViz ----

    #[test]
    fn test_interface_width_uniform_field() {
        let f = OrderParameterField::new(5, 5, 1.0, 1.0, 0.5);
        let result = compute_interface_width(&f, 0.01, Colormap::Jet);
        // No gradients → interface cells = 0, width = infinity
        assert_eq!(result.interface_cells, 0);
        assert!(result.width.is_infinite());
    }

    #[test]
    fn test_interface_width_step_field() {
        let f = make_field(10, 10);
        let result = compute_interface_width(&f, 0.0, Colormap::Viridis);
        // Step field has gradient along x = 5 column
        assert!(result.interface_cells > 0, "should have interface cells");
        assert!(result.max_gradient > 0.0);
    }

    #[test]
    fn test_interface_width_rgba_buffer_size() {
        let f = make_field(6, 6);
        let result = compute_interface_width(&f, 0.0, Colormap::Jet);
        assert_eq!(result.gradient_rgba.len(), 6 * 6 * 4);
    }

    // ---- GrainStructure ----

    #[test]
    fn test_grain_structure_grain_count() {
        let ids = vec![0u32, 1, 2, 1, 2, 3];
        let gs = GrainStructure::new(6, 1, ids);
        assert_eq!(gs.grain_count(), 3);
    }

    #[test]
    fn test_grain_structure_area_fraction() {
        let ids = vec![1u32; 4];
        let gs = GrainStructure::new(4, 1, ids);
        assert!((gs.area_fraction(1) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_grain_structure_dominant_grain() {
        let ids = vec![1u32, 1, 2, 1, 2];
        let gs = GrainStructure::new(5, 1, ids);
        assert_eq!(gs.dominant_grain(), 1);
    }

    #[test]
    fn test_grain_rgba_background_is_black() {
        let ids = vec![0u32; 4];
        let gs = GrainStructure::new(4, 1, ids);
        let rgba = gs.to_rgba_by_grain_id();
        // Grain 0 → black (0, 0, 0, 255)
        assert_eq!(rgba[0], 0);
        assert_eq!(rgba[1], 0);
        assert_eq!(rgba[2], 0);
        assert_eq!(rgba[3], 255);
    }

    #[test]
    fn test_grain_rgba_nonzero_grain_nonblack() {
        let ids = vec![42u32];
        let gs = GrainStructure::new(1, 1, ids);
        let rgba = gs.to_rgba_by_grain_id();
        // Grain 42 should get a non-black hue-based colour
        let is_nonblack = rgba[0] > 0 || rgba[1] > 0 || rgba[2] > 0;
        assert!(is_nonblack, "non-zero grain should be coloured");
    }

    // ---- SolidificationFront ----

    #[test]
    fn test_extract_front_step_field() {
        let f = make_field(10, 10);
        let front = extract_solidification_front(&f, 0.5);
        assert!(
            !front.contour_cells.is_empty(),
            "step field must yield a front"
        );
    }

    #[test]
    fn test_extract_front_uniform_field_empty() {
        let f = OrderParameterField::new(5, 5, 1.0, 1.0, 1.0);
        let front = extract_solidification_front(&f, 0.5);
        assert!(front.contour_cells.is_empty(), "uniform field has no front");
    }

    #[test]
    fn test_render_front_overlay_correct_size() {
        let f = make_field(8, 8);
        let front = extract_solidification_front(&f, 0.5);
        let base = vec![0u8; 8 * 8 * 4];
        let overlay = render_front_overlay(&base, 8, 8, &front, Color::red());
        assert_eq!(overlay.len(), 8 * 8 * 4);
    }

    // ---- DendriteTips ----

    #[test]
    fn test_dendrite_tips_detect_empty_uniform() {
        let f = OrderParameterField::new(8, 8, 1.0, 1.0, 0.5);
        let tips = DendriteTips::detect(&f, 0.1, 1);
        assert!(tips.is_empty(), "uniform field should yield no tips");
    }

    #[test]
    fn test_dendrite_tips_detect_step_field_finds_tips() {
        let f = make_field(10, 10);
        let tips = DendriteTips::detect(&f, 0.0, 1);
        // Step field interface should produce some tips (local gradient maxima)
        // We only verify the render doesn't panic
        let base = vec![0u8; 10 * 10 * 4];
        let out = tips.render_markers(&base, 10, 10, 1, Color::green());
        assert_eq!(out.len(), 10 * 10 * 4);
    }

    // ---- ThermalPhaseOverlay ----

    #[test]
    fn test_thermal_overlay_output_size() {
        let f = make_field(6, 6);
        let thermal = vec![300.0f64; 36];
        let overlay = ThermalPhaseOverlay::new(6, 6, 0.5, Colormap::Viridis, Colormap::Jet);
        let buf = overlay.render(&f, &thermal);
        assert_eq!(buf.len(), 6 * 6 * 4);
    }

    #[test]
    fn test_thermal_overlay_alpha_zero_equals_phase() {
        let f = make_field(4, 4);
        let thermal = vec![500.0f64; 16];
        let ov = ThermalPhaseOverlay::new(4, 4, 0.0, Colormap::Jet, Colormap::Viridis);
        let buf = ov.render(&f, &thermal);
        let phase_buf = f.to_rgba_buffer(Colormap::Jet, true);
        // With alpha=0, result should equal phase image
        assert_eq!(buf, phase_buf);
    }

    // ---- CompositionFieldViz ----

    #[test]
    fn test_composition_viz_render_size() {
        let viz = CompositionFieldViz::new(5, 5, 0.0, 1.0, Colormap::Viridis);
        let data = vec![0.5f64; 25];
        let buf = viz.render(&data);
        assert_eq!(buf.len(), 25 * 4);
    }

    #[test]
    fn test_composition_mean() {
        let data = vec![0.2f64, 0.4, 0.6, 0.8];
        let m = CompositionFieldViz::mean_composition(&data);
        assert!((m - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_composition_histogram_bucket_sum() {
        let data: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let hist = CompositionFieldViz::histogram(&data, 0.0, 1.0, 10);
        let total: usize = hist.iter().sum();
        assert_eq!(total, 100);
    }

    #[test]
    fn test_composition_histogram_uniform_bins() {
        // 100 evenly spaced values → each of 10 bins gets 10 counts
        let data: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let hist = CompositionFieldViz::histogram(&data, 0.0, 1.0, 10);
        for (i, &count) in hist.iter().enumerate() {
            assert_eq!(count, 10, "bin {i} should have 10 counts, got {count}");
        }
    }

    // ---- MisorientationViz ----

    #[test]
    fn test_misorientation_angle_deg() {
        let m = GrainPairMisorientation::new(1, 2, std::f64::consts::PI / 6.0);
        assert!((m.angle_deg() - 30.0).abs() < 1e-9);
    }

    #[test]
    fn test_misorientation_map_size() {
        let ids = vec![1u32, 1, 2, 2];
        let gs = GrainStructure::new(4, 1, ids);
        let table = vec![GrainPairMisorientation::new(1, 2, 0.5)];
        let buf = render_misorientation_map(&gs, &table, 60.0, Color::black(), Colormap::Jet);
        assert_eq!(buf.len(), 4 * 4);
    }

    // ---- PhaseFractionPie ----

    #[test]
    fn test_phase_fraction_pie_total() {
        let f = make_field(10, 10);
        let bands = [
            ("solid", 0.5, 1.01, Color::red()),
            ("liquid", -0.01, 0.5, Color::blue()),
        ];
        let pie = PhaseFractionPie::from_order_parameter(&f, &bands);
        let total = pie.total_fraction();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "fractions should sum to 1, got {total}"
        );
    }

    #[test]
    fn test_phase_fraction_pie_dominant() {
        let f = OrderParameterField::new(10, 10, 1.0, 1.0, 0.9);
        let bands = [
            ("solid", 0.5, 1.01, Color::red()),
            ("liquid", -0.01, 0.5, Color::blue()),
        ];
        let pie = PhaseFractionPie::from_order_parameter(&f, &bands);
        assert_eq!(pie.dominant_phase_index(), Some(0));
    }

    #[test]
    fn test_phase_fraction_arc_angles_span_2pi() {
        let f = make_field(10, 10);
        let bands = [
            ("solid", 0.5, 1.01, Color::red()),
            ("liquid", -0.01, 0.5, Color::blue()),
        ];
        let pie = PhaseFractionPie::from_order_parameter(&f, &bands);
        let arcs = pie.arc_segments();
        let span = arcs.last().map(|(_, e, _)| *e).unwrap_or(0.0);
        assert!((span - 2.0 * std::f64::consts::PI).abs() < 1e-9);
    }

    // ---- ConvergenceMonitor ----

    #[test]
    fn test_convergence_monitor_record_and_converge() {
        let mut mon = ConvergenceMonitor::new(1e-4);
        mon.record(0, 0.0, 1.0, -10.0, 0.1);
        mon.record(1, 0.01, 1e-5, -10.5, 0.15);
        assert!(mon.converged);
        assert_eq!(mon.len(), 2);
    }

    #[test]
    fn test_convergence_monitor_no_convergence() {
        let mut mon = ConvergenceMonitor::new(1e-10);
        mon.record(0, 0.0, 1.0, 0.0, 0.0);
        assert!(!mon.converged);
    }

    #[test]
    fn test_convergence_monitor_min_residual() {
        let mut mon = ConvergenceMonitor::new(1e-8);
        mon.record(0, 0.0, 1.0, 0.0, 0.0);
        mon.record(1, 1.0, 0.1, 0.0, 0.0);
        mon.record(2, 2.0, 0.01, 0.0, 0.0);
        assert!((mon.min_residual() - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_convergence_monitor_convergence_rate_negative() {
        let mut mon = ConvergenceMonitor::new(1e-12);
        for i in 0u64..10 {
            let res = (-0.5 * i as f64).exp();
            mon.record(i, i as f64, res, 0.0, 0.0);
        }
        let rate = mon.convergence_rate();
        assert!(
            rate < 0.0,
            "converging sequence should have negative rate, got {rate}"
        );
    }

    #[test]
    fn test_convergence_monitor_csv_header() {
        let mon = ConvergenceMonitor::new(1e-6);
        let csv = mon.to_csv();
        assert!(csv.starts_with("iteration,time,residual,free_energy,solid_fraction"));
    }

    #[test]
    fn test_convergence_monitor_render_plot_size() {
        let mut mon = ConvergenceMonitor::new(1e-6);
        for i in 0u64..5 {
            mon.record(i, i as f64, 1.0 / (i + 1) as f64, 0.0, 0.0);
        }
        let buf = mon.render_residual_plot(64, 32, Color::black(), Color::white());
        assert_eq!(buf.len(), 64 * 32 * 4);
    }

    #[test]
    fn test_convergence_monitor_render_empty_is_empty_vec() {
        let mon = ConvergenceMonitor::new(1e-6);
        let buf = mon.render_residual_plot(64, 32, Color::black(), Color::white());
        assert!(buf.is_empty());
    }

    // ---- hue_to_rgb helpers ----

    #[test]
    fn test_hue_to_rgb_red() {
        let [r, g, b] = hue_to_rgb(0.0);
        assert!((r - 1.0).abs() < 1e-6);
        assert!(g.abs() < 1e-6);
        assert!(b.abs() < 1e-6);
    }

    #[test]
    fn test_hue_to_rgb_green() {
        // hue = 1/3 → green
        let [r, g, b] = hue_to_rgb(1.0 / 3.0);
        assert!(g > r, "green channel should dominate");
        assert!(g > b, "green channel should dominate");
    }

    // ---- hash_u32_to_f64 ----

    #[test]
    fn test_hash_u32_in_unit_interval() {
        for v in [0u32, 1, 42, u32::MAX / 2, u32::MAX] {
            let h = hash_u32_to_f64(v);
            assert!((0.0..1.0).contains(&h), "hash({v}) = {h} out of [0,1)");
        }
    }

    #[test]
    fn test_hash_u32_different_seeds_differ() {
        let h1 = hash_u32_to_f64(1);
        let h2 = hash_u32_to_f64(2);
        assert!(
            (h1 - h2).abs() > 1e-6,
            "different seeds should produce different hashes"
        );
    }

    // ---- lerp / clamp ----

    #[test]
    fn test_lerp_endpoints() {
        assert!((lerp(0.0, 10.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((lerp(0.0, 10.0, 1.0) - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_midpoint() {
        assert!((lerp(0.0, 4.0, 0.5) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_clamp_below() {
        assert!((clamp(-5.0, 0.0, 1.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_clamp_above() {
        assert!((clamp(5.0, 0.0, 1.0) - 1.0).abs() < 1e-12);
    }
}
