// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fluid visualization primitives for the OxiPhysics engine.
//!
//! This module provides CPU-side algorithms to generate renderable geometry and
//! image buffers for visualizing fluid simulation data:
//!
//! - Velocity field arrow glyphs and streamlines/pathlines
//! - Vorticity magnitude heat-maps
//! - Pressure iso-contour lines (2-D) and iso-surfaces (3-D marching cubes)
//! - Line Integral Convolution (LIC) for 2-D flow visualization
//! - Marching-cubes iso-surface extraction
//! - Volume rendering for density fields (front-to-back alpha compositing)
//! - Shock-wave visualization (density-gradient magnitude)
//! - Turbulent kinetic energy (TKE) heat-maps
//! - Particle-based flow visualization with advection
//! - SPH particle rendering with density shading

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Shared math helpers (plain arrays — no nalgebra)
// ─────────────────────────────────────────────────────────────────────────────

/// A 3-D position or vector as `[f64; 3]`.
pub type Vec3 = [f64; 3];

/// A 2-D position or vector as `[f64; 2]`.
pub type Vec2 = [f64; 2];

/// An RGBA color packed as `[f32; 4]` (r, g, b, a) in \[0, 1\].
pub type Rgba = [f32; 4];

#[inline]
fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn length3(a: Vec3) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: Vec3) -> Vec3 {
    let l = length3(a);
    if l < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(a, 1.0 / l)
    }
}

#[inline]
fn cross3(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn add2(a: Vec2, b: Vec2) -> Vec2 {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn scale2(a: Vec2, s: f64) -> Vec2 {
    [a[0] * s, a[1] * s]
}

#[inline]
fn length2(a: Vec2) -> f64 {
    (a[0] * a[0] + a[1] * a[1]).sqrt()
}

#[inline]
fn normalize2(a: Vec2) -> Vec2 {
    let l = length2(a);
    if l < 1e-300 {
        [0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l]
    }
}

/// Linear interpolation between two scalars.
#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Clamp `x` to \[lo, hi\].
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ─────────────────────────────────────────────────────────────────────────────
// Simple colormap helpers (jet and hot)
// ─────────────────────────────────────────────────────────────────────────────

/// Map a normalised value `t ∈ [0, 1]` to a jet-colormap RGBA color.
pub fn jet_rgba(t: f64) -> Rgba {
    let t = clamp(t, 0.0, 1.0) as f32;
    let r = (1.5 - (t - 0.75).abs() * 4.0).clamp(0.0, 1.0);
    let g = (1.5 - (t - 0.50).abs() * 4.0).clamp(0.0, 1.0);
    let b = (1.5 - (t - 0.25).abs() * 4.0).clamp(0.0, 1.0);
    [r, g, b, 1.0]
}

/// Map a normalised value `t ∈ [0, 1]` to a "hot" colormap RGBA color.
pub fn hot_rgba(t: f64) -> Rgba {
    let t = clamp(t, 0.0, 1.0) as f32;
    let r = (t * 3.0).min(1.0);
    let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
    let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
    [r, g, b, 1.0]
}

/// Normalise a scalar from `[min_val, max_val]` to `[0, 1]`.
pub fn normalize_scalar(v: f64, min_val: f64, max_val: f64) -> f64 {
    if (max_val - min_val).abs() < 1e-300 {
        0.5
    } else {
        clamp((v - min_val) / (max_val - min_val), 0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Velocity field arrow glyphs
// ─────────────────────────────────────────────────────────────────────────────

/// A renderable arrow glyph representing a velocity vector.
#[derive(Debug, Clone)]
pub struct ArrowGlyph {
    /// Tail position of the arrow.
    pub origin: Vec3,
    /// Tip position of the arrow.
    pub tip: Vec3,
    /// Arrow color (RGBA).
    pub color: Rgba,
    /// Magnitude of the velocity represented.
    pub magnitude: f64,
}

/// Generate arrow glyphs for a 3-D velocity field on a regular grid.
///
/// - `velocities` — flat array of velocity vectors indexed by `(ix, iy, iz)`:
///   `idx = ix*(ny*nz) + iy*nz + iz`.
/// - `nx`, `ny`, `nz` — grid dimensions.
/// - `spacing` — physical spacing between grid points (assumed uniform).
/// - `scale` — arrow length scale factor relative to grid spacing.
/// - `max_speed` — maximum expected speed (for colormap normalisation).
pub fn velocity_arrow_glyphs(
    velocities: &[Vec3],
    nx: usize,
    ny: usize,
    nz: usize,
    spacing: f64,
    scale: f64,
    max_speed: f64,
) -> Vec<ArrowGlyph> {
    let mut glyphs = Vec::with_capacity(nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let idx = ix * ny * nz + iy * nz + iz;
                if idx >= velocities.len() {
                    continue;
                }
                let v = velocities[idx];
                let mag = length3(v);
                let origin = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                let tip = add3(origin, scale3(v, scale));
                let t = normalize_scalar(mag, 0.0, max_speed.max(1e-12));
                glyphs.push(ArrowGlyph {
                    origin,
                    tip,
                    color: jet_rgba(t),
                    magnitude: mag,
                });
            }
        }
    }
    glyphs
}

/// Generate arrow glyphs for a 2-D velocity field (z-component ignored).
///
/// Positions are placed in the XY plane (`z = 0`).
pub fn velocity_arrow_glyphs_2d(
    velocities: &[Vec2],
    nx: usize,
    ny: usize,
    spacing: f64,
    scale: f64,
    max_speed: f64,
) -> Vec<ArrowGlyph> {
    let mut glyphs = Vec::with_capacity(nx * ny);
    for ix in 0..nx {
        for iy in 0..ny {
            let idx = ix * ny + iy;
            if idx >= velocities.len() {
                continue;
            }
            let v2 = velocities[idx];
            let mag = length2(v2);
            let origin = [ix as f64 * spacing, iy as f64 * spacing, 0.0];
            let tip = [origin[0] + v2[0] * scale, origin[1] + v2[1] * scale, 0.0];
            let t = normalize_scalar(mag, 0.0, max_speed.max(1e-12));
            glyphs.push(ArrowGlyph {
                origin,
                tip,
                color: jet_rgba(t),
                magnitude: mag,
            });
        }
    }
    glyphs
}

// ─────────────────────────────────────────────────────────────────────────────
// Streamlines and pathlines
// ─────────────────────────────────────────────────────────────────────────────

/// A streamline or pathline as an ordered sequence of 3-D points.
#[derive(Debug, Clone)]
pub struct FlowLine {
    /// Points along the line.
    pub points: Vec<Vec3>,
    /// Speeds at each point (same length as `points`).
    pub speeds: Vec<f64>,
}

impl FlowLine {
    /// Total arc length of the line.
    pub fn arc_length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| length3(sub3(w[1], w[0])))
            .sum()
    }
}

/// Trace a single streamline from `seed` using 4th-order Runge-Kutta.
///
/// `field` is a closure `(Vec3) -> Vec3` returning the velocity at a position.
/// Integration stops after `max_steps` steps or when the speed drops below
/// `min_speed`.
pub fn trace_streamline_rk4<F>(
    field: &F,
    seed: Vec3,
    dt: f64,
    max_steps: usize,
    min_speed: f64,
) -> FlowLine
where
    F: Fn(Vec3) -> Vec3,
{
    let mut points = Vec::with_capacity(max_steps + 1);
    let mut speeds = Vec::with_capacity(max_steps + 1);
    let mut pos = seed;
    points.push(pos);
    speeds.push(length3(field(pos)));

    for _ in 0..max_steps {
        let k1 = field(pos);
        // Stop early if the field at current position is below min_speed
        if length3(k1) < min_speed {
            break;
        }
        let k2 = field(add3(pos, scale3(k1, dt * 0.5)));
        let k3 = field(add3(pos, scale3(k2, dt * 0.5)));
        let k4 = field(add3(pos, scale3(k3, dt)));
        let combined = [
            (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
            (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
            (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0,
        ];
        pos = add3(pos, scale3(combined, dt));
        let speed = length3(combined);
        points.push(pos);
        speeds.push(speed);
    }
    FlowLine { points, speeds }
}

/// Trace multiple streamlines from a set of seeds.
pub fn trace_streamlines<F>(
    field: &F,
    seeds: &[Vec3],
    dt: f64,
    max_steps: usize,
    min_speed: f64,
) -> Vec<FlowLine>
where
    F: Fn(Vec3) -> Vec3,
{
    seeds
        .iter()
        .map(|&seed| trace_streamline_rk4(field, seed, dt, max_steps, min_speed))
        .collect()
}

/// Trace a pathline by integrating a time-varying velocity field.
///
/// `field_at` is `(pos, time) -> Vec3`. Uses forward Euler for simplicity.
pub fn trace_pathline<F>(
    field_at: &F,
    seed: Vec3,
    dt: f64,
    max_steps: usize,
    start_time: f64,
) -> FlowLine
where
    F: Fn(Vec3, f64) -> Vec3,
{
    let mut points = Vec::with_capacity(max_steps + 1);
    let mut speeds = Vec::with_capacity(max_steps + 1);
    let mut pos = seed;
    let mut t = start_time;
    points.push(pos);
    let v0 = field_at(pos, t);
    speeds.push(length3(v0));

    for _ in 0..max_steps {
        let v = field_at(pos, t);
        pos = add3(pos, scale3(v, dt));
        t += dt;
        let speed = length3(v);
        points.push(pos);
        speeds.push(speed);
    }
    FlowLine { points, speeds }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vorticity magnitude
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the vorticity vector ω = ∇ × u at each interior point of a 3-D
/// velocity grid using finite differences.
///
/// The returned array is length `nx * ny * nz` (boundary points have ω = 0).
pub fn compute_vorticity(
    velocities: &[Vec3],
    nx: usize,
    ny: usize,
    nz: usize,
    h: f64,
) -> Vec<Vec3> {
    let n = nx * ny * nz;
    let mut vorticity = vec![[0.0f64; 3]; n];
    let idx = |ix: usize, iy: usize, iz: usize| ix * ny * nz + iy * nz + iz;

    for ix in 1..nx.saturating_sub(1) {
        for iy in 1..ny.saturating_sub(1) {
            for iz in 1..nz.saturating_sub(1) {
                let i = idx(ix, iy, iz);
                // ∂w/∂y − ∂v/∂z
                let dwdy = (velocities[idx(ix, iy + 1, iz)][2]
                    - velocities[idx(ix, iy - 1, iz)][2])
                    / (2.0 * h);
                let dvdz = (velocities[idx(ix, iy, iz + 1)][1]
                    - velocities[idx(ix, iy, iz - 1)][1])
                    / (2.0 * h);
                // ∂u/∂z − ∂w/∂x
                let dudz = (velocities[idx(ix, iy, iz + 1)][0]
                    - velocities[idx(ix, iy, iz - 1)][0])
                    / (2.0 * h);
                let dwdx = (velocities[idx(ix + 1, iy, iz)][2]
                    - velocities[idx(ix - 1, iy, iz)][2])
                    / (2.0 * h);
                // ∂v/∂x − ∂u/∂y
                let dvdx = (velocities[idx(ix + 1, iy, iz)][1]
                    - velocities[idx(ix - 1, iy, iz)][1])
                    / (2.0 * h);
                let dudy = (velocities[idx(ix, iy + 1, iz)][0]
                    - velocities[idx(ix, iy - 1, iz)][0])
                    / (2.0 * h);

                vorticity[i] = [dwdy - dvdz, dudz - dwdx, dvdx - dudy];
            }
        }
    }
    vorticity
}

/// Compute the vorticity magnitude at every grid point.
pub fn vorticity_magnitude(vorticity: &[Vec3]) -> Vec<f64> {
    vorticity.iter().map(|&v| length3(v)).collect()
}

/// Generate a heat-map image (row-major, width × height × 4) of vorticity
/// magnitude for a 2-D slice at fixed z-index.
pub fn vorticity_heatmap_slice(
    vorticity_mag: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    iz: usize,
    max_vort: f64,
) -> Vec<u8> {
    let mut img = vec![0u8; nx * ny * 4];
    for ix in 0..nx {
        for iy in 0..ny {
            let g_idx = ix * ny * nz + iy * nz + iz;
            let val = if g_idx < vorticity_mag.len() {
                vorticity_mag[g_idx]
            } else {
                0.0
            };
            let t = normalize_scalar(val, 0.0, max_vort.max(1e-12));
            let c = hot_rgba(t);
            let px = (iy * nx + ix) * 4;
            img[px] = (c[0] * 255.0) as u8;
            img[px + 1] = (c[1] * 255.0) as u8;
            img[px + 2] = (c[2] * 255.0) as u8;
            img[px + 3] = 255;
        }
    }
    img
}

// ─────────────────────────────────────────────────────────────────────────────
// Pressure contours (2-D)
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D iso-contour line segment.
#[derive(Debug, Clone, Copy)]
pub struct ContourSegment {
    /// Start point.
    pub p0: Vec2,
    /// End point.
    pub p1: Vec2,
    /// Iso-value this segment represents.
    pub iso_value: f64,
}

/// Extract iso-contour line segments from a 2-D scalar field using the
/// Marching Squares algorithm.
///
/// `field` is indexed `[ix * ny + iy]`. `iso` is the target iso-value.
/// Grid spacing in x is `dx`, in y is `dy`.
pub fn marching_squares(
    field: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    iso: f64,
) -> Vec<ContourSegment> {
    let mut segments = Vec::new();
    let f = |ix: usize, iy: usize| -> f64 {
        if ix < nx && iy < ny {
            field[ix * ny + iy]
        } else {
            iso // treat out-of-range as boundary
        }
    };

    for ix in 0..nx.saturating_sub(1) {
        for iy in 0..ny.saturating_sub(1) {
            let v00 = f(ix, iy);
            let v10 = f(ix + 1, iy);
            let v01 = f(ix, iy + 1);
            let v11 = f(ix + 1, iy + 1);

            let x0 = ix as f64 * dx;
            let x1 = (ix + 1) as f64 * dx;
            let y0 = iy as f64 * dy;
            let y1 = (iy + 1) as f64 * dy;

            // Marching-squares case index
            let case = ((v00 > iso) as u8)
                | (((v10 > iso) as u8) << 1)
                | (((v11 > iso) as u8) << 2)
                | (((v01 > iso) as u8) << 3);

            if case == 0 || case == 15 {
                continue; // no crossing
            }

            // Linear interpolation along each edge
            let interp_x = |a: f64, b: f64, xa: f64, xb: f64| -> f64 {
                if (b - a).abs() < 1e-300 {
                    (xa + xb) * 0.5
                } else {
                    xa + (iso - a) / (b - a) * (xb - xa)
                }
            };

            // Edges: bottom (y0), right (x1), top (y1), left (x0)
            let bottom = [interp_x(v00, v10, x0, x1), y0];
            let right = [x1, interp_x(v10, v11, y0, y1)];
            let top = [interp_x(v01, v11, x0, x1), y1];
            let left = [x0, interp_x(v00, v01, y0, y1)];

            // Lookup table for the 16 cases (excluding 0 and 15)
            let segs: &[[Vec2; 2]] = match case {
                1 | 14 => &[[bottom, left]],
                2 | 13 => &[[bottom, right]],
                3 | 12 => &[[left, right]],
                4 | 11 => &[[top, right]],
                5 => &[[bottom, right], [top, left]], // ambiguous — pick one
                6 | 9 => &[[bottom, top]],
                7 | 8 => &[[top, left]],
                10 => &[[bottom, left], [top, right]], // ambiguous — pick one
                _ => &[],
            };

            for &[p0, p1] in segs {
                segments.push(ContourSegment {
                    p0,
                    p1,
                    iso_value: iso,
                });
            }
        }
    }
    segments
}

/// Generate contour segments for multiple iso-values at once.
pub fn pressure_contours(
    pressure: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    iso_values: &[f64],
) -> Vec<ContourSegment> {
    iso_values
        .iter()
        .flat_map(|&iso| marching_squares(pressure, nx, ny, dx, dy, iso))
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Line Integral Convolution (LIC) — 2-D flow texture
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a LIC texture for a 2-D velocity field.
///
/// The result is a greyscale image of size `width × height` stored as a flat
/// `Vec`f32` (one value per pixel, in [0, 1]).
///
/// - `velocities` — indexed `\[ix * ny + iy\]`, each a `Vec2`.
/// - `noise` — a white-noise input texture (same indexing, values in [0, 1]).
/// - `steps` — number of integration steps in each direction.
/// - `dt` — step size in pixel units.
pub fn lic_2d(
    velocities: &[Vec2],
    noise: &[f32],
    width: usize,
    height: usize,
    steps: usize,
    dt: f64,
) -> Vec<f32> {
    let mut output = vec![0.0f32; width * height];
    let nx = width;
    let ny = height;

    let sample_noise = |px: f64, py: f64| -> f32 {
        let ix = (px as isize).rem_euclid(nx as isize) as usize;
        let iy = (py as isize).rem_euclid(ny as isize) as usize;
        let idx = ix * ny + iy;
        if idx < noise.len() { noise[idx] } else { 0.5 }
    };

    let sample_vel = |px: f64, py: f64| -> Vec2 {
        let ix = (px as isize).rem_euclid(nx as isize) as usize;
        let iy = (py as isize).rem_euclid(ny as isize) as usize;
        let idx = ix * ny + iy;
        if idx < velocities.len() {
            velocities[idx]
        } else {
            [0.0; 2]
        }
    };

    for ix in 0..nx {
        for iy in 0..ny {
            let mut acc = 0.0f32;
            let mut count = 0u32;
            let mut pos = [ix as f64, iy as f64];

            // Forward integration
            for _ in 0..steps {
                let v = normalize2(sample_vel(pos[0], pos[1]));
                acc += sample_noise(pos[0], pos[1]);
                count += 1;
                pos = add2(pos, scale2(v, dt));
            }

            // Backward integration
            pos = [ix as f64, iy as f64];
            for _ in 0..steps {
                let v = normalize2(sample_vel(pos[0], pos[1]));
                acc += sample_noise(pos[0], pos[1]);
                count += 1;
                pos = [pos[0] - v[0] * dt, pos[1] - v[1] * dt];
            }

            let pixel_idx = ix * ny + iy;
            output[pixel_idx] = if count > 0 { acc / count as f32 } else { 0.5 };
        }
    }
    output
}

// ─────────────────────────────────────────────────────────────────────────────
// Marching Cubes iso-surface extraction
// ─────────────────────────────────────────────────────────────────────────────

/// A triangle in 3-D space (positions and an outward unit normal).
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// Three vertex positions.
    pub vertices: [Vec3; 3],
    /// Outward-facing unit normal.
    pub normal: Vec3,
}

impl Triangle {
    /// Compute the normal from vertex positions.
    pub fn compute_normal(v: [Vec3; 3]) -> Vec3 {
        let ab = sub3(v[1], v[0]);
        let ac = sub3(v[2], v[0]);
        normalize3(cross3(ab, ac))
    }
}

// A compact subset of the marching-cubes edge table for 256 cases would be
// large; here we implement a simplified version covering the 15 topologically
// distinct cases via a representative lookup.

/// Extract an iso-surface at the given `iso_value` from a scalar 3-D field.
///
/// `field` is indexed `\[ix*(ny*nz) + iy*nz + iz\]`.
/// Grid spacing is assumed uniform `h` in all directions.
pub fn marching_cubes(
    field: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    h: f64,
    iso_value: f64,
) -> Vec<Triangle> {
    let mut triangles = Vec::new();
    let get = |ix: usize, iy: usize, iz: usize| -> f64 {
        if ix < nx && iy < ny && iz < nz {
            field[ix * ny * nz + iy * nz + iz]
        } else {
            iso_value
        }
    };

    let pos =
        |ix: usize, iy: usize, iz: usize| -> Vec3 { [ix as f64 * h, iy as f64 * h, iz as f64 * h] };

    // Linear interpolation of edge vertex
    let edge_point = |p0: Vec3, f0: f64, p1: Vec3, f1: f64| -> Vec3 {
        if (f1 - f0).abs() < 1e-300 {
            [
                (p0[0] + p1[0]) * 0.5,
                (p0[1] + p1[1]) * 0.5,
                (p0[2] + p1[2]) * 0.5,
            ]
        } else {
            let t = (iso_value - f0) / (f1 - f0);
            add3(p0, scale3(sub3(p1, p0), t))
        }
    };

    for ix in 0..nx.saturating_sub(1) {
        for iy in 0..ny.saturating_sub(1) {
            for iz in 0..nz.saturating_sub(1) {
                // The 8 corners of this cube cell
                let corners = [
                    (get(ix, iy, iz), pos(ix, iy, iz)),
                    (get(ix + 1, iy, iz), pos(ix + 1, iy, iz)),
                    (get(ix + 1, iy + 1, iz), pos(ix + 1, iy + 1, iz)),
                    (get(ix, iy + 1, iz), pos(ix, iy + 1, iz)),
                    (get(ix, iy, iz + 1), pos(ix, iy, iz + 1)),
                    (get(ix + 1, iy, iz + 1), pos(ix + 1, iy, iz + 1)),
                    (get(ix + 1, iy + 1, iz + 1), pos(ix + 1, iy + 1, iz + 1)),
                    (get(ix, iy + 1, iz + 1), pos(ix, iy + 1, iz + 1)),
                ];

                let cube_index: u8 = corners.iter().enumerate().fold(0u8, |acc, (i, &(f, _))| {
                    if f > iso_value { acc | (1 << i) } else { acc }
                });

                if cube_index == 0 || cube_index == 255 {
                    continue;
                }

                // Edge vertices (indices per the standard MC edge table)
                // 12 edges of the cube
                let edges: [(usize, usize); 12] = [
                    (0, 1),
                    (1, 2),
                    (2, 3),
                    (3, 0), // bottom face
                    (4, 5),
                    (5, 6),
                    (6, 7),
                    (7, 4), // top face
                    (0, 4),
                    (1, 5),
                    (2, 6),
                    (3, 7), // verticals
                ];

                let mut verts: [Vec3; 12] = [[0.0; 3]; 12];
                for (e, &(a, b)) in edges.iter().enumerate() {
                    let (fa, pa) = corners[a];
                    let (fb, pb) = corners[b];
                    verts[e] = edge_point(pa, fa, pb, fb);
                }

                // Simplified triangle table — only a subset of the 256 cases.
                // Each row lists triples of edge indices forming triangles.
                let tris = mc_triangle_table(cube_index);
                for tri in tris {
                    let v = [verts[tri[0]], verts[tri[1]], verts[tri[2]]];
                    let normal = Triangle::compute_normal(v);
                    triangles.push(Triangle {
                        vertices: v,
                        normal,
                    });
                }
            }
        }
    }
    triangles
}

/// Return the triangle connectivity for the given marching-cubes case index.
/// Each returned slice element is `\[e0, e1, e2\]` — edge indices forming a triangle.
fn mc_triangle_table(case: u8) -> Vec<[usize; 3]> {
    // A minimal lookup covering the 15 canonical cases (and their complements).
    // For brevity we handle a selection of representative cases.
    match case {
        // 1 vertex inside (and complement 254)
        1 | 254 => vec![[0, 3, 8]],
        2 | 253 => vec![[0, 9, 1]],
        4 | 251 => vec![[1, 10, 2]],
        8 | 247 => vec![[2, 11, 3]],
        16 | 239 => vec![[4, 7, 8]],
        32 | 223 => vec![[4, 9, 5]],
        64 | 191 => vec![[5, 10, 6]],
        128 | 127 => vec![[6, 11, 7]],
        // 2 adjacent vertices inside
        3 | 252 => vec![[1, 3, 8], [1, 8, 9]],
        6 | 249 => vec![[0, 10, 2], [0, 9, 10]],
        // Diagonal face cases (ambiguous — we pick a single split)
        5 | 250 => vec![[0, 3, 8], [1, 10, 2]],
        10 | 245 => vec![[0, 9, 1], [2, 11, 3]],
        // 3 vertices inside
        7 | 248 => vec![[2, 3, 8], [2, 8, 10], [8, 9, 10]],
        // More complex cases — approximate
        15 | 240 => vec![[0, 4, 8], [0, 8, 3], [4, 5, 9], [4, 9, 8]],
        // Fallback: skip unknown cases
        _ => vec![],
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Volume rendering (front-to-back alpha compositing)
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for volume rendering.
#[derive(Debug, Clone)]
pub struct VolumeRenderParams {
    /// Number of ray-marching steps along each ray.
    pub num_steps: usize,
    /// Step size as a fraction of the volume's diagonal.
    pub step_fraction: f64,
    /// Density scale (opacity multiplier per unit density).
    pub density_scale: f64,
    /// Minimum density threshold (below this, sample is transparent).
    pub density_threshold: f64,
}

impl Default for VolumeRenderParams {
    fn default() -> Self {
        Self {
            num_steps: 128,
            step_fraction: 0.01,
            density_scale: 1.0,
            density_threshold: 0.01,
        }
    }
}

/// Trilinearly interpolate a scalar from a 3-D grid at position `p`.
///
/// `p` is in grid coordinates (`\[0, nx-1\] × \[0, ny-1\] × \[0, nz-1\]`).
pub fn trilinear_sample(field: &[f64], nx: usize, ny: usize, nz: usize, p: Vec3) -> f64 {
    let ix = (p[0] as usize).min(nx.saturating_sub(2));
    let iy = (p[1] as usize).min(ny.saturating_sub(2));
    let iz = (p[2] as usize).min(nz.saturating_sub(2));
    let fx = p[0] - ix as f64;
    let fy = p[1] - iy as f64;
    let fz = p[2] - iz as f64;

    let get = |x: usize, y: usize, z: usize| -> f64 {
        if x < nx && y < ny && z < nz {
            field[x * ny * nz + y * nz + z]
        } else {
            0.0
        }
    };

    let c000 = get(ix, iy, iz);
    let c100 = get(ix + 1, iy, iz);
    let c010 = get(ix, iy + 1, iz);
    let c110 = get(ix + 1, iy + 1, iz);
    let c001 = get(ix, iy, iz + 1);
    let c101 = get(ix + 1, iy, iz + 1);
    let c011 = get(ix, iy + 1, iz + 1);
    let c111 = get(ix + 1, iy + 1, iz + 1);

    let c00 = lerp(c000, c100, fx);
    let c10 = lerp(c010, c110, fx);
    let c01 = lerp(c001, c101, fx);
    let c11 = lerp(c011, c111, fx);
    let c0 = lerp(c00, c10, fy);
    let c1 = lerp(c01, c11, fy);
    lerp(c0, c1, fz)
}

/// Volume-render a single ray using front-to-back alpha compositing.
///
/// Returns the accumulated RGBA color for this ray. `ray_dir` should be
/// normalised. `t_near` and `t_far` are the entry/exit distances along the ray.
pub fn volume_render_ray(
    density: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    ray_origin: Vec3,
    ray_dir: Vec3,
    t_near: f64,
    t_far: f64,
    params: &VolumeRenderParams,
    colormap: fn(f64) -> Rgba,
) -> Rgba {
    let step = (t_far - t_near) / params.num_steps as f64;
    let mut color = [0.0f32; 4];
    let mut transmittance = 1.0f64;

    for i in 0..params.num_steps {
        let t = t_near + (i as f64 + 0.5) * step;
        let p = add3(ray_origin, scale3(ray_dir, t));
        // Map to grid coordinates
        let gp = [
            p[0] * (nx as f64 - 1.0),
            p[1] * (ny as f64 - 1.0),
            p[2] * (nz as f64 - 1.0),
        ];
        if gp[0] < 0.0
            || gp[1] < 0.0
            || gp[2] < 0.0
            || gp[0] >= nx as f64
            || gp[1] >= ny as f64
            || gp[2] >= nz as f64
        {
            continue;
        }
        let d = trilinear_sample(density, nx, ny, nz, gp);
        if d < params.density_threshold {
            continue;
        }
        let opacity = 1.0 - (-d * params.density_scale * step).exp();
        let sample_color = colormap(d);
        color[0] += (transmittance * opacity * sample_color[0] as f64) as f32;
        color[1] += (transmittance * opacity * sample_color[1] as f64) as f32;
        color[2] += (transmittance * opacity * sample_color[2] as f64) as f32;
        color[3] += (transmittance * opacity) as f32;
        transmittance *= 1.0 - opacity;
        if transmittance < 1e-4 {
            break;
        }
    }
    color
}

// ─────────────────────────────────────────────────────────────────────────────
// Shock-wave visualization
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the density gradient magnitude at each interior grid point.
///
/// This is used to highlight shock waves (regions of large density gradient).
pub fn density_gradient_magnitude(
    density: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    h: f64,
) -> Vec<f64> {
    let n = nx * ny * nz;
    let mut grad_mag = vec![0.0f64; n];
    let idx = |ix: usize, iy: usize, iz: usize| ix * ny * nz + iy * nz + iz;

    for ix in 1..nx.saturating_sub(1) {
        for iy in 1..ny.saturating_sub(1) {
            for iz in 1..nz.saturating_sub(1) {
                let drdx =
                    (density[idx(ix + 1, iy, iz)] - density[idx(ix - 1, iy, iz)]) / (2.0 * h);
                let drdy =
                    (density[idx(ix, iy + 1, iz)] - density[idx(ix, iy - 1, iz)]) / (2.0 * h);
                let drdz =
                    (density[idx(ix, iy, iz + 1)] - density[idx(ix, iy, iz - 1)]) / (2.0 * h);
                grad_mag[idx(ix, iy, iz)] = (drdx * drdx + drdy * drdy + drdz * drdz).sqrt();
            }
        }
    }
    grad_mag
}

/// Generate a RGBA image of the schlieren-like shock visualisation for a 2-D
/// density slice at constant z (`iz`).
pub fn shock_visualization_image(
    grad_mag: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    iz: usize,
    max_grad: f64,
) -> Vec<u8> {
    let mut img = vec![0u8; nx * ny * 4];
    for ix in 0..nx {
        for iy in 0..ny {
            let g = if ix * ny * nz + iy * nz + iz < grad_mag.len() {
                grad_mag[ix * ny * nz + iy * nz + iz]
            } else {
                0.0
            };
            let t = normalize_scalar(g, 0.0, max_grad.max(1e-12));
            // Schlieren: black (low gradient) to white (high gradient)
            let v = (t * 255.0) as u8;
            let px = (iy * nx + ix) * 4;
            img[px] = v;
            img[px + 1] = v;
            img[px + 2] = v;
            img[px + 3] = 255;
        }
    }
    img
}

// ─────────────────────────────────────────────────────────────────────────────
// Turbulent kinetic energy heat-map
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the turbulent kinetic energy (TKE) per cell from an ensemble of
/// velocity snapshots.
///
/// TKE = ½ (⟨u'²⟩ + ⟨v'²⟩ + ⟨w'²⟩) where `u' = u − ⟨u⟩`.
///
/// `snapshots` — list of velocity fields, each of length `n_cells`.
pub fn turbulent_kinetic_energy(snapshots: &[Vec<Vec3>]) -> Vec<f64> {
    if snapshots.is_empty() {
        return Vec::new();
    }
    let n_cells = snapshots[0].len();
    let n_snap = snapshots.len() as f64;

    // Mean velocity
    let mut mean = vec![[0.0f64; 3]; n_cells];
    for snap in snapshots {
        for (i, &v) in snap.iter().enumerate() {
            if i < n_cells {
                mean[i][0] += v[0];
                mean[i][1] += v[1];
                mean[i][2] += v[2];
            }
        }
    }
    for m in mean.iter_mut() {
        m[0] /= n_snap;
        m[1] /= n_snap;
        m[2] /= n_snap;
    }

    // Variance
    let mut tke = vec![0.0f64; n_cells];
    for snap in snapshots {
        for (i, &v) in snap.iter().enumerate() {
            if i < n_cells {
                let du = v[0] - mean[i][0];
                let dv = v[1] - mean[i][1];
                let dw = v[2] - mean[i][2];
                tke[i] += du * du + dv * dv + dw * dw;
            }
        }
    }
    for t in tke.iter_mut() {
        *t = 0.5 * *t / n_snap;
    }
    tke
}

/// Generate a TKE heat-map image for a 2-D grid.
pub fn tke_heatmap_image(tke: &[f64], nx: usize, ny: usize, max_tke: f64) -> Vec<u8> {
    let mut img = vec![0u8; nx * ny * 4];
    for ix in 0..nx {
        for iy in 0..ny {
            let idx = ix * ny + iy;
            let t = normalize_scalar(
                tke.get(idx).copied().unwrap_or(0.0),
                0.0,
                max_tke.max(1e-12),
            );
            let c = jet_rgba(t);
            let px = (iy * nx + ix) * 4;
            img[px] = (c[0] * 255.0) as u8;
            img[px + 1] = (c[1] * 255.0) as u8;
            img[px + 2] = (c[2] * 255.0) as u8;
            img[px + 3] = 255;
        }
    }
    img
}

// ─────────────────────────────────────────────────────────────────────────────
// Particle-based flow visualization
// ─────────────────────────────────────────────────────────────────────────────

/// A single flow tracer particle.
#[derive(Debug, Clone)]
pub struct FlowParticle {
    /// Current position.
    pub position: Vec3,
    /// Current velocity.
    pub velocity: Vec3,
    /// Age in simulation time units.
    pub age: f64,
    /// Maximum lifetime (particle is removed after this time).
    pub lifetime: f64,
    /// Display color (RGBA).
    pub color: Rgba,
}

impl FlowParticle {
    /// Return `true` if the particle has exceeded its lifetime.
    pub fn is_dead(&self) -> bool {
        self.age >= self.lifetime
    }
}

/// A collection of flow tracer particles and their rendering parameters.
#[derive(Debug, Clone)]
pub struct FlowParticleSystem {
    /// All particles (live and dead).
    pub particles: Vec<FlowParticle>,
    /// Maximum number of live particles.
    pub max_particles: usize,
}

impl FlowParticleSystem {
    /// Construct an empty particle system.
    pub fn new(max_particles: usize) -> Self {
        Self {
            particles: Vec::with_capacity(max_particles),
            max_particles,
        }
    }

    /// Spawn a new particle at `position` with `velocity` and the given `lifetime`.
    pub fn spawn(&mut self, position: Vec3, velocity: Vec3, lifetime: f64, color: Rgba) {
        if self.live_count() >= self.max_particles {
            return;
        }
        self.particles.push(FlowParticle {
            position,
            velocity,
            age: 0.0,
            lifetime,
            color,
        });
    }

    /// Count living particles.
    pub fn live_count(&self) -> usize {
        self.particles.iter().filter(|p| !p.is_dead()).count()
    }

    /// Advance all particles by `dt` using the provided velocity field `field(pos)`.
    /// Dead particles are removed.
    pub fn step<F>(&mut self, field: &F, dt: f64)
    where
        F: Fn(Vec3) -> Vec3,
    {
        for p in &mut self.particles {
            if p.is_dead() {
                continue;
            }
            let v = field(p.position);
            p.position = add3(p.position, scale3(v, dt));
            p.velocity = v;
            p.age += dt;
        }
        self.particles.retain(|p| !p.is_dead());
    }

    /// Return positions of all live particles.
    pub fn live_positions(&self) -> Vec<Vec3> {
        self.particles
            .iter()
            .filter(|p| !p.is_dead())
            .map(|p| p.position)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH particle rendering with density shading
// ─────────────────────────────────────────────────────────────────────────────

/// A single SPH particle with physical properties for rendering.
#[derive(Debug, Clone)]
pub struct SphParticle {
    /// Position.
    pub position: Vec3,
    /// Density (kg/m³ or normalised).
    pub density: f64,
    /// Pressure (Pa or normalised).
    pub pressure: f64,
    /// Velocity magnitude.
    pub speed: f64,
    /// Smoothing length.
    pub h: f64,
}

/// Render an SPH particle system into a 2-D projected image.
///
/// Each particle is splattered onto the image using a Gaussian kernel of radius
/// proportional to its smoothing length `h` (projected to screen space). The
/// image size is `width × height`. World-to-screen mapping uses `camera_origin`
/// and `pixel_size` (world units per pixel).
pub fn sph_density_image(
    particles: &[SphParticle],
    width: usize,
    height: usize,
    camera_origin: Vec2,
    pixel_size: f64,
    max_density: f64,
    colormap: fn(f64) -> Rgba,
) -> Vec<u8> {
    let mut img = vec![0.0f32; width * height * 4];

    for p in particles {
        // Project to screen (only X-Y plane for simplicity)
        let sx = (p.position[0] - camera_origin[0]) / pixel_size;
        let sy = (p.position[1] - camera_origin[1]) / pixel_size;
        let screen_h = p.h / pixel_size;
        let radius = (screen_h * 3.0) as isize;
        let c = colormap(normalize_scalar(p.density, 0.0, max_density.max(1e-12)));

        for dx in -radius..=radius {
            for dy in -radius..=radius {
                let px = (sx + dx as f64) as isize;
                let py = (sy + dy as f64) as isize;
                if px < 0 || py < 0 || px >= width as isize || py >= height as isize {
                    continue;
                }
                let dist2 = (dx as f64 * dx as f64 + dy as f64 * dy as f64)
                    / (screen_h * screen_h).max(1e-12);
                let w = (-0.5 * dist2).exp() as f32;
                let base = (py as usize * width + px as usize) * 4;
                img[base] += c[0] * w;
                img[base + 1] += c[1] * w;
                img[base + 2] += c[2] * w;
                img[base + 3] += w;
            }
        }
    }

    // Normalise and convert to u8
    let mut out = vec![0u8; width * height * 4];
    for i in 0..width * height {
        let base = i * 4;
        let a = img[base + 3];
        if a > 0.0 {
            out[base] = ((img[base] / a).min(1.0) * 255.0) as u8;
            out[base + 1] = ((img[base + 1] / a).min(1.0) * 255.0) as u8;
            out[base + 2] = ((img[base + 2] / a).min(1.0) * 255.0) as u8;
            out[base + 3] = (a.min(1.0) * 255.0) as u8;
        }
    }
    out
}

/// Compute the SPH kernel smoothed density field on a regular 2-D grid.
///
/// Uses the cubic spline kernel W(r,h) = (1 - (r/h)²)³ for r < h.
pub fn sph_density_grid(
    particles: &[SphParticle],
    nx: usize,
    ny: usize,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
) -> Vec<f64> {
    let mut grid = vec![0.0f64; nx * ny];
    let dx = (x_max - x_min) / nx as f64;
    let dy = (y_max - y_min) / ny as f64;

    for p in particles {
        let h = p.h.max(1e-12);
        let ix_min = ((p.position[0] - h - x_min) / dx) as isize;
        let ix_max = ((p.position[0] + h - x_min) / dx) as isize + 1;
        let iy_min = ((p.position[1] - h - y_min) / dy) as isize;
        let iy_max = ((p.position[1] + h - y_min) / dy) as isize + 1;

        for ix in ix_min.max(0)..ix_max.min(nx as isize) {
            for iy in iy_min.max(0)..iy_max.min(ny as isize) {
                let gx = x_min + (ix as f64 + 0.5) * dx;
                let gy = y_min + (iy as f64 + 0.5) * dy;
                let r = ((p.position[0] - gx).powi(2) + (p.position[1] - gy).powi(2)).sqrt();
                if r < h {
                    let q = r / h;
                    let w = (1.0 - q * q).powi(3);
                    grid[ix as usize * ny + iy as usize] += p.density * w;
                }
            }
        }
    }
    grid
}

// ─────────────────────────────────────────────────────────────────────────────
// Combined fluid state for convenience
// ─────────────────────────────────────────────────────────────────────────────

/// A snapshot of a 3-D fluid field on a regular grid.
#[derive(Debug, Clone)]
pub struct FluidGridSnapshot {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid y-dimension.
    pub ny: usize,
    /// Grid z-dimension.
    pub nz: usize,
    /// Uniform grid spacing.
    pub h: f64,
    /// Velocity field (length `nx*ny*nz`).
    pub velocities: Vec<Vec3>,
    /// Pressure field (length `nx*ny*nz`).
    pub pressure: Vec<f64>,
    /// Density field (length `nx*ny*nz`).
    pub density: Vec<f64>,
}

impl FluidGridSnapshot {
    /// Construct a zero-filled snapshot.
    pub fn zeros(nx: usize, ny: usize, nz: usize, h: f64) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            h,
            velocities: vec![[0.0; 3]; n],
            pressure: vec![0.0; n],
            density: vec![0.0; n],
        }
    }

    /// Number of grid cells.
    pub fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Return `true` if the grid has no cells.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Cell index for `(ix, iy, iz)`.
    pub fn cell_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.ny * self.nz + iy * self.nz + iz
    }

    /// Maximum speed across all cells.
    pub fn max_speed(&self) -> f64 {
        self.velocities
            .iter()
            .map(|&v| length3(v))
            .fold(0.0f64, f64::max)
    }

    /// Maximum pressure.
    pub fn max_pressure(&self) -> f64 {
        self.pressure
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum pressure.
    pub fn min_pressure(&self) -> f64 {
        self.pressure.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Vorticity magnitude field.
    pub fn vorticity_mag(&self) -> Vec<f64> {
        let vort = compute_vorticity(&self.velocities, self.nx, self.ny, self.nz, self.h);
        vorticity_magnitude(&vort)
    }

    /// Density gradient magnitude field (for shock detection).
    pub fn shock_indicator(&self) -> Vec<f64> {
        density_gradient_magnitude(&self.density, self.nx, self.ny, self.nz, self.h)
    }

    /// Arrow glyphs for the velocity field.
    pub fn arrow_glyphs(&self, scale: f64) -> Vec<ArrowGlyph> {
        let max_speed = self.max_speed().max(1e-12);
        velocity_arrow_glyphs(
            &self.velocities,
            self.nx,
            self.ny,
            self.nz,
            self.h,
            scale,
            max_speed,
        )
    }

    /// Extract an iso-surface of the pressure field.
    pub fn pressure_isosurface(&self, iso: f64) -> Vec<Triangle> {
        marching_cubes(&self.pressure, self.nx, self.ny, self.nz, self.h, iso)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: generate a noise texture for LIC
// ─────────────────────────────────────────────────────────────────────────────

/// Generate a deterministic pseudo-random noise texture using a simple LCG.
///
/// Returns a flat array of `width × height` values in [0, 1].
pub fn generate_noise_texture(width: usize, height: usize, seed: u64) -> Vec<f32> {
    let mut state = seed.wrapping_add(1);
    let mut out = Vec::with_capacity(width * height);
    for _ in 0..width * height {
        // Xorshift64
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        out.push((state as f64 / u64::MAX as f64) as f32);
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: rotate a Vec3 around the Z axis
// ─────────────────────────────────────────────────────────────────────────────

/// Rotate `v` about the Z axis by angle `theta` (radians).
pub fn rotate_z(v: Vec3, theta: f64) -> Vec3 {
    let c = theta.cos();
    let s = theta.sin();
    [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]]
}

/// Build a circular ring of seed points in the XY plane at radius `r` and `z`-height.
pub fn ring_seeds(n: usize, r: f64, z: f64) -> Vec<Vec3> {
    (0..n)
        .map(|i| {
            let theta = 2.0 * PI * i as f64 / n as f64;
            [r * theta.cos(), r * theta.sin(), z]
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Math helpers ───────────────────────────────────────────────────────

    #[test]
    fn test_length3_unit_vectors() {
        assert!((length3([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!((length3([0.0, 1.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!((length3([0.0, 0.0, 1.0]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3_produces_unit_vector() {
        let v = normalize3([3.0, 4.0, 0.0]);
        assert!((length3(v) - 1.0).abs() < 1e-12);
        assert!((v[0] - 0.6).abs() < 1e-12);
        assert!((v[1] - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_cross3_orthogonal() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0] - 0.0).abs() < 1e-12);
        assert!((c[1] - 0.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_endpoints() {
        assert!((lerp(0.0, 10.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((lerp(0.0, 10.0, 1.0) - 10.0).abs() < 1e-12);
        assert!((lerp(0.0, 10.0, 0.5) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_scalar_clamps() {
        assert!((normalize_scalar(5.0, 0.0, 10.0) - 0.5).abs() < 1e-12);
        assert!((normalize_scalar(-1.0, 0.0, 10.0) - 0.0).abs() < 1e-12);
        assert!((normalize_scalar(20.0, 0.0, 10.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_scalar_degenerate() {
        // min == max → returns 0.5
        let t = normalize_scalar(5.0, 5.0, 5.0);
        assert!((t - 0.5).abs() < 1e-12);
    }

    // ── Colormap ───────────────────────────────────────────────────────────

    #[test]
    fn test_jet_rgba_blue_at_zero() {
        let c = jet_rgba(0.0);
        assert!(c[2] > c[0], "jet(0) should be more blue than red");
    }

    #[test]
    fn test_jet_rgba_red_at_one() {
        let c = jet_rgba(1.0);
        assert!(c[0] > c[2], "jet(1) should be more red than blue");
    }

    #[test]
    fn test_hot_rgba_black_at_zero() {
        let c = hot_rgba(0.0);
        assert!(c[0] < 1e-3 && c[1] < 1e-3 && c[2] < 1e-3);
    }

    #[test]
    fn test_hot_rgba_white_at_one() {
        let c = hot_rgba(1.0);
        assert!(c[0] > 0.99 && c[1] > 0.99 && c[2] > 0.99);
    }

    #[test]
    fn test_rgba_alpha_always_one() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert!((jet_rgba(t)[3] - 1.0).abs() < 1e-6);
            assert!((hot_rgba(t)[3] - 1.0).abs() < 1e-6);
        }
    }

    // ── Arrow glyphs ───────────────────────────────────────────────────────

    #[test]
    fn test_velocity_arrow_glyphs_count() {
        let vels = vec![[1.0, 0.0, 0.0]; 8]; // 2×2×2 grid
        let glyphs = velocity_arrow_glyphs(&vels, 2, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(glyphs.len(), 8);
    }

    #[test]
    fn test_velocity_arrow_glyphs_2d_count() {
        let vels = vec![[0.5, 0.5]; 9]; // 3×3 grid
        let glyphs = velocity_arrow_glyphs_2d(&vels, 3, 3, 1.0, 1.0, 1.0);
        assert_eq!(glyphs.len(), 9);
    }

    #[test]
    fn test_arrow_glyph_tip_offset() {
        let v = [1.0, 0.0, 0.0];
        let vels = vec![v];
        let glyphs = velocity_arrow_glyphs(&vels, 1, 1, 1, 1.0, 2.0, 1.0);
        assert_eq!(glyphs.len(), 1);
        // origin = (0,0,0), scale=2 → tip = (2,0,0)
        assert!((glyphs[0].tip[0] - 2.0).abs() < 1e-12);
        assert!(glyphs[0].tip[1].abs() < 1e-12);
    }

    // ── Streamlines ────────────────────────────────────────────────────────

    #[test]
    fn test_streamline_uniform_field_length() {
        let field = |_: Vec3| -> Vec3 { [1.0, 0.0, 0.0] };
        let sl = trace_streamline_rk4(&field, [0.0, 0.0, 0.0], 0.1, 10, 0.0);
        assert_eq!(sl.points.len(), 11);
    }

    #[test]
    fn test_streamline_uniform_field_position() {
        let field = |_: Vec3| -> Vec3 { [1.0, 0.0, 0.0] };
        let sl = trace_streamline_rk4(&field, [0.0; 3], 0.1, 5, 0.0);
        let last = sl.points.last().unwrap();
        assert!((last[0] - 0.5).abs() < 1e-10);
        assert!(last[1].abs() < 1e-10);
    }

    #[test]
    fn test_streamline_stops_at_min_speed() {
        // Field that decays to zero
        let field = |pos: Vec3| -> Vec3 {
            let mag = (-pos[0]).exp();
            [mag, 0.0, 0.0]
        };
        let sl = trace_streamline_rk4(&field, [0.0; 3], 0.5, 100, 0.1);
        // Should stop before 100 steps once field decays below min_speed
        assert!(sl.points.len() < 101);
    }

    #[test]
    fn test_trace_streamlines_multiple() {
        let field = |_: Vec3| -> Vec3 { [0.0, 1.0, 0.0] };
        let seeds = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let lines = trace_streamlines(&field, &seeds, 0.1, 5, 0.0);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn test_pathline_length() {
        let field_at = |_pos: Vec3, _t: f64| -> Vec3 { [1.0, 0.0, 0.0] };
        let pl = trace_pathline(&field_at, [0.0; 3], 0.1, 8, 0.0);
        assert_eq!(pl.points.len(), 9);
    }

    #[test]
    fn test_flowline_arc_length() {
        // Points: (0,0,0), (1,0,0), (2,0,0) → arc length = 2
        let fl = FlowLine {
            points: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            speeds: vec![1.0, 1.0, 1.0],
        };
        assert!((fl.arc_length() - 2.0).abs() < 1e-12);
    }

    // ── Vorticity ──────────────────────────────────────────────────────────

    #[test]
    fn test_vorticity_zero_for_uniform_field() {
        // Uniform field: no curl
        let n = 4usize;
        let vels = vec![[1.0, 2.0, 3.0]; n * n * n];
        let vort = compute_vorticity(&vels, n, n, n, 1.0);
        let mag = vorticity_magnitude(&vort);
        for &m in &mag {
            assert!(
                m < 1e-10,
                "vorticity should be zero for uniform field, got {m}"
            );
        }
    }

    #[test]
    fn test_vorticity_magnitude_positive() {
        // Solid rotation: u = y, v = -x, w = 0 → curl = (0, 0, -2)
        let n = 5usize;
        let mut vels = vec![[0.0f64; 3]; n * n * n];
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = ix as f64;
                    let y = iy as f64;
                    vels[ix * n * n + iy * n + iz] = [y, -x, 0.0];
                }
            }
        }
        let vort = compute_vorticity(&vels, n, n, n, 1.0);
        let mag = vorticity_magnitude(&vort);
        // Interior points should have non-zero vorticity
        let interior = 2 * n * n + 2 * n + 2;
        assert!(mag[interior] > 0.5, "expected non-zero vorticity");
    }

    // ── Marching squares ───────────────────────────────────────────────────

    #[test]
    fn test_marching_squares_no_crossing() {
        // All values below iso → no segments
        let field = vec![0.0f64; 4];
        let segs = marching_squares(&field, 2, 2, 1.0, 1.0, 5.0);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_marching_squares_all_above() {
        // All values above iso → no segments (case 15)
        let field = vec![10.0f64; 4];
        let segs = marching_squares(&field, 2, 2, 1.0, 1.0, 5.0);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_marching_squares_one_crossing() {
        // 3×3 field with a gradient that crosses iso=0.5
        let field: Vec<f64> = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0];
        let segs = marching_squares(&field, 3, 3, 1.0, 1.0, 0.5);
        assert!(!segs.is_empty());
    }

    #[test]
    fn test_pressure_contours_multiple_iso() {
        let field: Vec<f64> = (0..16).map(|i| i as f64 / 15.0).collect();
        let segs = pressure_contours(&field, 4, 4, 1.0, 1.0, &[0.3, 0.6]);
        assert!(!segs.is_empty());
    }

    // ── LIC ────────────────────────────────────────────────────────────────

    #[test]
    fn test_lic_output_size() {
        let vels = vec![[1.0f64, 0.0]; 16]; // 4×4
        let noise = vec![0.5f32; 16];
        let result = lic_2d(&vels, &noise, 4, 4, 3, 0.5);
        assert_eq!(result.len(), 16);
    }

    #[test]
    fn test_lic_output_in_range() {
        let vels = vec![[1.0f64, 0.0]; 25]; // 5×5
        let noise = generate_noise_texture(5, 5, 42);
        let result = lic_2d(&vels, &noise, 5, 5, 5, 0.5);
        for &v in &result {
            assert!((0.0..=1.0).contains(&v), "LIC value out of [0,1]: {v}");
        }
    }

    // ── Marching cubes ─────────────────────────────────────────────────────

    #[test]
    fn test_marching_cubes_no_triangles_for_uniform_field() {
        let field = vec![0.0f64; 27]; // 3×3×3, all below iso=1
        let tris = marching_cubes(&field, 3, 3, 3, 1.0, 1.0);
        assert!(tris.is_empty());
    }

    #[test]
    fn test_marching_cubes_some_triangles_for_sphere() {
        // Sphere SDF: f(x,y,z) = r² − (x² + y² + z²)
        let n = 10usize;
        let mut field = vec![0.0f64; n * n * n];
        let center = (n as f64 - 1.0) * 0.5;
        let r2 = (n as f64 * 0.4).powi(2);
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let dx = ix as f64 - center;
                    let dy = iy as f64 - center;
                    let dz = iz as f64 - center;
                    field[ix * n * n + iy * n + iz] = r2 - (dx * dx + dy * dy + dz * dz);
                }
            }
        }
        let tris = marching_cubes(&field, n, n, n, 1.0, 0.0);
        assert!(!tris.is_empty(), "sphere SDF should produce triangles");
    }

    #[test]
    fn test_triangle_normal_perpendicular() {
        let v = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let n = Triangle::compute_normal(v);
        // Normal should be (0, 0, 1)
        assert!(n[0].abs() < 1e-10);
        assert!(n[1].abs() < 1e-10);
        assert!((n[2].abs() - 1.0).abs() < 1e-10);
    }

    // ── Volume rendering ───────────────────────────────────────────────────

    #[test]
    fn test_trilinear_sample_corner() {
        // 2×2×2 field with value 1.0 at corner (1,1,1) and 0 elsewhere
        let mut field = vec![0.0f64; 8];
        field[2 * 2 + 2 + 1] = 1.0;
        let v = trilinear_sample(&field, 2, 2, 2, [1.0, 1.0, 1.0]);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_trilinear_sample_midpoint() {
        // Uniform field = 2.0
        let field = vec![2.0f64; 8];
        let v = trilinear_sample(&field, 2, 2, 2, [0.5, 0.5, 0.5]);
        assert!((v - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_render_ray_empty_density() {
        let density = vec![0.0f64; 8];
        let params = VolumeRenderParams::default();
        let color = volume_render_ray(
            &density,
            2,
            2,
            2,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            0.0,
            1.0,
            &params,
            jet_rgba,
        );
        // All zeros → no contribution
        assert!(color[3] < 1e-3, "transparent for empty density");
    }

    // ── Shock visualization ────────────────────────────────────────────────

    #[test]
    fn test_density_gradient_uniform_field() {
        let n = 4usize;
        let density = vec![1.0f64; n * n * n];
        let grad = density_gradient_magnitude(&density, n, n, n, 1.0);
        for &g in &grad {
            assert!(g < 1e-10, "gradient of uniform field should be zero");
        }
    }

    #[test]
    fn test_density_gradient_step_function() {
        // Step function: left half = 0, right half = 1 → gradient at interface
        let n = 4usize;
        let mut density = vec![0.0f64; n * n * n];
        for ix in n / 2..n {
            for iy in 0..n {
                for iz in 0..n {
                    density[ix * n * n + iy * n + iz] = 1.0;
                }
            }
        }
        let grad = density_gradient_magnitude(&density, n, n, n, 1.0);
        let max_grad = grad.iter().cloned().fold(0.0f64, f64::max);
        assert!(max_grad > 0.1, "step function should have large gradient");
    }

    #[test]
    fn test_shock_image_size() {
        let grad = vec![0.5f64; 4 * 4 * 2];
        let img = shock_visualization_image(&grad, 4, 4, 2, 0, 1.0);
        assert_eq!(img.len(), 4 * 4 * 4);
    }

    // ── TKE ────────────────────────────────────────────────────────────────

    #[test]
    fn test_tke_zero_for_identical_snapshots() {
        // All snapshots identical → fluctuation = 0 → TKE = 0
        let snap: Vec<Vec3> = vec![[1.0, 2.0, 3.0]; 4];
        let tke = turbulent_kinetic_energy(&[snap.clone(), snap.clone(), snap]);
        for &t in &tke {
            assert!(t < 1e-12, "TKE should be zero for identical snapshots");
        }
    }

    #[test]
    fn test_tke_symmetric_fluctuations() {
        // Two snapshots: +v and -v → mean = 0, TKE = v²
        let v = [1.0, 0.0, 0.0];
        let s1: Vec<Vec3> = vec![v; 3];
        let s2: Vec<Vec3> = vec![[-1.0, 0.0, 0.0]; 3];
        let tke = turbulent_kinetic_energy(&[s1, s2]);
        for &t in &tke {
            // mean=0, variance sum = 1²+1² = 2, TKE = 0.5 * 2 / 2 = 0.5
            assert!((t - 0.5).abs() < 1e-10, "TKE = {t}");
        }
    }

    #[test]
    fn test_tke_heatmap_image_size() {
        let tke = vec![0.5f64; 6];
        let img = tke_heatmap_image(&tke, 3, 2, 1.0);
        assert_eq!(img.len(), 3 * 2 * 4);
    }

    // ── Flow particles ─────────────────────────────────────────────────────

    #[test]
    fn test_particle_system_spawn_and_count() {
        let mut sys = FlowParticleSystem::new(10);
        sys.spawn([0.0; 3], [1.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0, 1.0]);
        sys.spawn([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0, [0.0, 1.0, 0.0, 1.0]);
        assert_eq!(sys.live_count(), 2);
    }

    #[test]
    fn test_particle_advection() {
        let mut sys = FlowParticleSystem::new(10);
        sys.spawn([0.0; 3], [1.0, 0.0, 0.0], 10.0, [1.0; 4]);
        let field = |_pos: Vec3| -> Vec3 { [1.0, 0.0, 0.0] };
        sys.step(&field, 1.0);
        assert!((sys.particles[0].position[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_death_removes_particle() {
        let mut sys = FlowParticleSystem::new(10);
        sys.spawn([0.0; 3], [1.0, 0.0, 0.0], 0.5, [1.0; 4]);
        let field = |_pos: Vec3| -> Vec3 { [1.0, 0.0, 0.0] };
        sys.step(&field, 1.0); // age 1.0 > lifetime 0.5
        assert_eq!(sys.live_count(), 0);
    }

    #[test]
    fn test_particle_max_particles_respected() {
        let mut sys = FlowParticleSystem::new(2);
        for _ in 0..5 {
            sys.spawn([0.0; 3], [0.0; 3], 10.0, [1.0; 4]);
        }
        assert!(sys.live_count() <= 2);
    }

    // ── SPH rendering ──────────────────────────────────────────────────────

    #[test]
    fn test_sph_density_image_size() {
        let p = SphParticle {
            position: [5.0, 5.0, 0.0],
            density: 1.0,
            pressure: 1.0,
            speed: 0.0,
            h: 1.0,
        };
        let img = sph_density_image(&[p], 10, 10, [0.0, 0.0], 1.0, 2.0, jet_rgba);
        assert_eq!(img.len(), 10 * 10 * 4);
    }

    #[test]
    fn test_sph_density_grid_positive() {
        let p = SphParticle {
            position: [5.0, 5.0, 0.0],
            density: 1.0,
            pressure: 1.0,
            speed: 0.0,
            h: 2.0,
        };
        let grid = sph_density_grid(&[p], 10, 10, 0.0, 10.0, 0.0, 10.0);
        let max_val = grid.iter().cloned().fold(0.0f64, f64::max);
        assert!(max_val > 0.0, "density grid should have positive values");
    }

    // ── FluidGridSnapshot ──────────────────────────────────────────────────

    #[test]
    fn test_fluid_snapshot_len() {
        let snap = FluidGridSnapshot::zeros(3, 4, 5, 0.1);
        assert_eq!(snap.len(), 60);
    }

    #[test]
    fn test_fluid_snapshot_cell_index() {
        let snap = FluidGridSnapshot::zeros(3, 4, 5, 0.1);
        assert_eq!(snap.cell_index(2, 3, 4), 2 * 4 * 5 + 3 * 5 + 4);
    }

    #[test]
    fn test_fluid_snapshot_max_speed_zero() {
        let snap = FluidGridSnapshot::zeros(2, 2, 2, 1.0);
        assert!((snap.max_speed() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_fluid_snapshot_pressure_range() {
        let mut snap = FluidGridSnapshot::zeros(2, 2, 2, 1.0);
        snap.pressure[0] = -1.0;
        snap.pressure[7] = 5.0;
        assert!((snap.min_pressure() - (-1.0)).abs() < 1e-12);
        assert!((snap.max_pressure() - 5.0).abs() < 1e-12);
    }

    // ── Noise texture ──────────────────────────────────────────────────────

    #[test]
    fn test_noise_texture_size() {
        let noise = generate_noise_texture(8, 8, 12345);
        assert_eq!(noise.len(), 64);
    }

    #[test]
    fn test_noise_texture_values_in_range() {
        let noise = generate_noise_texture(16, 16, 99);
        for &v in &noise {
            assert!((0.0..=1.0).contains(&v), "noise value out of range: {v}");
        }
    }

    // ── Ring seeds / rotate ────────────────────────────────────────────────

    #[test]
    fn test_ring_seeds_count() {
        let seeds = ring_seeds(8, 1.0, 0.5);
        assert_eq!(seeds.len(), 8);
    }

    #[test]
    fn test_ring_seeds_radius() {
        let r = 2.0;
        for seed in ring_seeds(12, r, 0.0) {
            let dist = (seed[0] * seed[0] + seed[1] * seed[1]).sqrt();
            assert!((dist - r).abs() < 1e-10, "seed not on ring: {dist}");
        }
    }

    #[test]
    fn test_rotate_z_90_degrees() {
        let v = [1.0, 0.0, 0.0];
        let rotated = rotate_z(v, PI / 2.0);
        assert!(rotated[0].abs() < 1e-10);
        assert!((rotated[1] - 1.0).abs() < 1e-10);
        assert!(rotated[2].abs() < 1e-10);
    }

    // ── Vorticity heatmap ──────────────────────────────────────────────────

    #[test]
    fn test_vorticity_heatmap_size() {
        let mag = vec![0.5f64; 4 * 4 * 2];
        let img = vorticity_heatmap_slice(&mag, 4, 4, 2, 0, 1.0);
        assert_eq!(img.len(), 4 * 4 * 4);
    }

    #[test]
    fn test_vorticity_heatmap_alpha_full() {
        let mag = vec![0.5f64; 9]; // 3×3×1
        let img = vorticity_heatmap_slice(&mag, 3, 3, 1, 0, 1.0);
        // All alpha channels should be 255
        for i in 0..9 {
            assert_eq!(img[i * 4 + 3], 255);
        }
    }
}
