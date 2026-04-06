#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Flow visualization primitives for the OxiPhysics engine.
//!
//! This module provides:
//! - [`VectorField`] — 2D/3D vector field, interpolation, divergence, curl.
//! - [`LicRenderer`] — Line Integral Convolution for 2D flow visualization.
//! - [`PathlineTracer`] — particle pathlines, backward tracing, bundles.
//! - [`VorticityVisualization`] — vorticity magnitude, Q-criterion, λ2, swirling strength.
//! - [`PressureField`] — pressure contours, iso-surfaces, pressure coefficient Cp.
//! - [`FlowStatistics`] — time-averaged fields, RMS, energy spectra.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vec2 / Vec3 helpers (local, no external dependency)
// ---------------------------------------------------------------------------

/// A 2D vector of `f64`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec2 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
}

impl Vec2 {
    /// Construct from components.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    /// Euclidean length.
    pub fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Return a unit vector, or zero if length is negligible.
    pub fn normalize(self) -> Self {
        let l = self.length();
        if l < 1e-15 {
            Self::zero()
        } else {
            Self::new(self.x / l, self.y / l)
        }
    }

    /// Dot product.
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y
    }

    /// 2D "cross" (scalar z-component).
    pub fn cross(self, rhs: Self) -> f64 {
        self.x * rhs.y - self.y * rhs.x
    }
}

impl std::ops::Add for Vec2 {
    type Output = Vec2;
    fn add(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl std::ops::Sub for Vec2 {
    type Output = Vec2;
    fn sub(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl std::ops::Mul<f64> for Vec2 {
    type Output = Vec2;
    fn mul(self, s: f64) -> Vec2 {
        Vec2::new(self.x * s, self.y * s)
    }
}

/// A 3D vector of `f64`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// Construct from components.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Euclidean length.
    pub fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Return a unit vector, or zero if length is negligible.
    pub fn normalize(self) -> Self {
        let l = self.length();
        if l < 1e-15 {
            Self::zero()
        } else {
            Self::new(self.x / l, self.y / l, self.z / l)
        }
    }

    /// Dot product.
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Cross product.
    pub fn cross(self, rhs: Self) -> Self {
        Self::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }
}

impl std::ops::Add for Vec3 {
    type Output = Vec3;
    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, rhs: Vec3) -> Vec3 {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl std::ops::Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
}

// ---------------------------------------------------------------------------
// VectorField
// ---------------------------------------------------------------------------

/// Dimensionality of a vector field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FieldDimension {
    /// 2-dimensional field (x, y vectors on a 2D grid).
    Two,
    /// 3-dimensional field (x, y, z vectors on a 3D grid).
    Three,
}

/// A discrete vector field defined on a regular Cartesian grid.
///
/// For 2D fields, the grid has `nx × ny` cells and `nz = 1`.
/// For 3D fields, the grid has `nx × ny × nz` cells.
///
/// Vectors are stored row-major: index `(ix, iy, iz)` → `iz * nx * ny + iy * nx + ix`.
#[derive(Debug, Clone)]
pub struct VectorField {
    /// Number of cells along X.
    pub nx: usize,
    /// Number of cells along Y.
    pub ny: usize,
    /// Number of cells along Z (1 for 2D).
    pub nz: usize,
    /// Physical spacing in X.
    pub dx: f64,
    /// Physical spacing in Y.
    pub dy: f64,
    /// Physical spacing in Z.
    pub dz: f64,
    /// Grid origin (min corner).
    pub origin: Vec3,
    /// Vector components: `u[i] = (ux, uy, uz)` at flat index `i`.
    pub u: Vec<Vec3>,
    /// Field dimensionality.
    pub dim: FieldDimension,
}

impl VectorField {
    /// Construct a 2D vector field initialized to zero.
    pub fn new_2d(nx: usize, ny: usize, dx: f64, dy: f64, origin: Vec2) -> Self {
        Self {
            nx,
            ny,
            nz: 1,
            dx,
            dy,
            dz: 1.0,
            origin: Vec3::new(origin.x, origin.y, 0.0),
            u: vec![Vec3::zero(); nx * ny],
            dim: FieldDimension::Two,
        }
    }

    /// Construct a 3D vector field initialized to zero.
    pub fn new_3d(
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        dy: f64,
        dz: f64,
        origin: Vec3,
    ) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            origin,
            u: vec![Vec3::zero(); nx * ny * nz],
            dim: FieldDimension::Three,
        }
    }

    /// Flat index for grid cell `(ix, iy, iz)`.
    #[inline]
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.nx * self.ny + iy * self.nx + ix
    }

    /// Set the vector at grid cell `(ix, iy, iz)`.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, v: Vec3) {
        let idx = self.index(ix, iy, iz);
        self.u[idx] = v;
    }

    /// Get the vector at grid cell `(ix, iy, iz)`.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> Vec3 {
        self.u[self.index(ix, iy, iz)]
    }

    /// Bilinear interpolation of the 2D vector field at physical position `(px, py)`.
    ///
    /// Returns the zero vector for out-of-domain queries.
    pub fn interpolate_2d(&self, px: f64, py: f64) -> Vec3 {
        let lx = px - self.origin.x;
        let ly = py - self.origin.y;
        let fi = lx / self.dx;
        let fj = ly / self.dy;
        let i0 = fi.floor() as isize;
        let j0 = fj.floor() as isize;
        let ti = fi - fi.floor();
        let tj = fj - fj.floor();

        let clamp = |v: isize, max: usize| -> usize { v.clamp(0, max as isize - 1) as usize };
        let i0c = clamp(i0, self.nx);
        let i1c = clamp(i0 + 1, self.nx);
        let j0c = clamp(j0, self.ny);
        let j1c = clamp(j0 + 1, self.ny);

        // Out-of-domain check
        if i0 < 0 || j0 < 0 || i0 + 1 > self.nx as isize || j0 + 1 > self.ny as isize {
            // Allow clamped lookup; field naturally extrapolates at boundaries
        }

        let v00 = self.get(i0c, j0c, 0);
        let v10 = self.get(i1c, j0c, 0);
        let v01 = self.get(i0c, j1c, 0);
        let v11 = self.get(i1c, j1c, 0);

        let lerp = |a: Vec3, b: Vec3, t: f64| a * (1.0 - t) + b * t;
        lerp(lerp(v00, v10, ti), lerp(v01, v11, ti), tj)
    }

    /// Trilinear interpolation of the 3D vector field at physical position `(px, py, pz)`.
    pub fn interpolate_3d(&self, px: f64, py: f64, pz: f64) -> Vec3 {
        let lx = px - self.origin.x;
        let ly = py - self.origin.y;
        let lz = pz - self.origin.z;
        let fi = lx / self.dx;
        let fj = ly / self.dy;
        let fk = lz / self.dz;
        let ti = fi - fi.floor();
        let tj = fj - fj.floor();
        let tk = fk - fk.floor();

        let clamp = |v: isize, max: usize| -> usize { v.clamp(0, max as isize - 1) as usize };
        let i0 = clamp(fi.floor() as isize, self.nx);
        let i1 = clamp(fi.floor() as isize + 1, self.nx);
        let j0 = clamp(fj.floor() as isize, self.ny);
        let j1 = clamp(fj.floor() as isize + 1, self.ny);
        let k0 = clamp(fk.floor() as isize, self.nz);
        let k1 = clamp(fk.floor() as isize + 1, self.nz);

        let lerp = |a: Vec3, b: Vec3, t: f64| a * (1.0 - t) + b * t;

        let v000 = self.get(i0, j0, k0);
        let v100 = self.get(i1, j0, k0);
        let v010 = self.get(i0, j1, k0);
        let v110 = self.get(i1, j1, k0);
        let v001 = self.get(i0, j0, k1);
        let v101 = self.get(i1, j0, k1);
        let v011 = self.get(i0, j1, k1);
        let v111 = self.get(i1, j1, k1);

        let c00 = lerp(v000, v100, ti);
        let c10 = lerp(v010, v110, ti);
        let c01 = lerp(v001, v101, ti);
        let c11 = lerp(v011, v111, ti);
        let c0 = lerp(c00, c10, tj);
        let c1 = lerp(c01, c11, tj);
        lerp(c0, c1, tk)
    }

    /// Compute the 2D divergence field `∂u/∂x + ∂v/∂y` using central differences.
    pub fn divergence_2d(&self) -> Vec<f64> {
        let n = self.nx * self.ny;
        let mut div = vec![0.0f64; n];
        for j in 0..self.ny {
            for i in 0..self.nx {
                let ip = (i + 1).min(self.nx - 1);
                let im = i.saturating_sub(1);
                let jp = (j + 1).min(self.ny - 1);
                let jm = j.saturating_sub(1);
                let du_dx = (self.get(ip, j, 0).x - self.get(im, j, 0).x)
                    / (2.0 * self.dx * (ip - im).max(1) as f64);
                let dv_dy = (self.get(i, jp, 0).y - self.get(i, jm, 0).y)
                    / (2.0 * self.dy * (jp - jm).max(1) as f64);
                div[self.index(i, j, 0)] = du_dx + dv_dy;
            }
        }
        div
    }

    /// Compute the 2D curl (z-component) `∂v/∂x - ∂u/∂y` using central differences.
    pub fn curl_2d(&self) -> Vec<f64> {
        let n = self.nx * self.ny;
        let mut curl = vec![0.0f64; n];
        for j in 0..self.ny {
            for i in 0..self.nx {
                let ip = (i + 1).min(self.nx - 1);
                let im = i.saturating_sub(1);
                let jp = (j + 1).min(self.ny - 1);
                let jm = j.saturating_sub(1);
                let dv_dx = (self.get(ip, j, 0).y - self.get(im, j, 0).y)
                    / (2.0 * self.dx * (ip - im).max(1) as f64);
                let du_dy = (self.get(i, jp, 0).x - self.get(i, jm, 0).x)
                    / (2.0 * self.dy * (jp - jm).max(1) as f64);
                curl[self.index(i, j, 0)] = dv_dx - du_dy;
            }
        }
        curl
    }

    /// Compute the 3D vorticity vector `∇ × u` at each grid point.
    pub fn vorticity_3d(&self) -> Vec<Vec3> {
        let n = self.nx * self.ny * self.nz;
        let mut omega = vec![Vec3::zero(); n];
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let ip = (i + 1).min(self.nx - 1);
                    let im = i.saturating_sub(1);
                    let jp = (j + 1).min(self.ny - 1);
                    let jm = j.saturating_sub(1);
                    let kp = (k + 1).min(self.nz - 1);
                    let km = k.saturating_sub(1);
                    let two_dx = 2.0 * self.dx * (ip - im).max(1) as f64;
                    let two_dy = 2.0 * self.dy * (jp - jm).max(1) as f64;
                    let two_dz = 2.0 * self.dz * (kp - km).max(1) as f64;
                    let dw_dy = (self.get(i, jp, k).z - self.get(i, jm, k).z) / two_dy;
                    let dv_dz = (self.get(i, j, kp).y - self.get(i, j, km).y) / two_dz;
                    let du_dz = (self.get(i, j, kp).x - self.get(i, j, km).x) / two_dz;
                    let dw_dx = (self.get(ip, j, k).z - self.get(im, j, k).z) / two_dx;
                    let dv_dx = (self.get(ip, j, k).y - self.get(im, j, k).y) / two_dx;
                    let du_dy = (self.get(i, jp, k).x - self.get(i, jm, k).x) / two_dy;
                    omega[self.index(i, j, k)] =
                        Vec3::new(dw_dy - dv_dz, du_dz - dw_dx, dv_dx - du_dy);
                }
            }
        }
        omega
    }

    /// Return the maximum velocity magnitude in the field.
    pub fn max_velocity(&self) -> f64 {
        self.u.iter().map(|v| v.length()).fold(0.0f64, f64::max)
    }

    /// Fill the field with a uniform velocity.
    pub fn fill_uniform(&mut self, vel: Vec3) {
        for v in &mut self.u {
            *v = vel;
        }
    }

    /// Fill with a 2D circular (vortex) flow: u=(−y, x) / r² centred at `centre`.
    pub fn fill_vortex_2d(&mut self, centre: Vec2, strength: f64) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let px = self.origin.x + i as f64 * self.dx;
                let py = self.origin.y + j as f64 * self.dy;
                let rx = px - centre.x;
                let ry = py - centre.y;
                let r2 = rx * rx + ry * ry;
                if r2 > 1e-15 {
                    let idx = self.index(i, j, 0);
                    self.u[idx] = Vec3::new(-ry * strength / r2, rx * strength / r2, 0.0);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// LicRenderer
// ---------------------------------------------------------------------------

/// Configuration for the Line Integral Convolution algorithm.
#[derive(Debug, Clone)]
pub struct LicConfig {
    /// Length (in pixels) of the integration kernel.
    pub kernel_length: usize,
    /// Step size for advection (in grid cells).
    pub step_size: f64,
    /// Contrast enhancement factor applied after convolution.
    pub contrast: f64,
}

impl Default for LicConfig {
    fn default() -> Self {
        Self {
            kernel_length: 20,
            step_size: 0.5,
            contrast: 1.5,
        }
    }
}

/// Performs Line Integral Convolution (LIC) on a 2D vector field.
///
/// The output is a grey-scale image of size `nx × ny`, stored row-major.
/// Each pixel value is in `[0.0, 1.0]`.
pub struct LicRenderer {
    /// LIC configuration.
    pub config: LicConfig,
}

impl LicRenderer {
    /// Create a renderer with default settings.
    pub fn new() -> Self {
        Self {
            config: LicConfig::default(),
        }
    }

    /// Create a renderer with the given configuration.
    pub fn with_config(config: LicConfig) -> Self {
        Self { config }
    }

    /// Render the LIC image from a vector field and a white-noise texture.
    ///
    /// - `field` — 2D vector field (only the XY components are used).
    /// - `noise` — white-noise texture of size `nx × ny` (values in `[0.0, 1.0]`).
    ///   If `None`, a deterministic pseudo-noise is generated internally.
    ///
    /// Returns the LIC output pixel values as a `Vec`f64` of length `nx × ny`.
    pub fn render(&self, field: &VectorField, noise: Option<&[f64]>) -> Vec<f64> {
        let nx = field.nx;
        let ny = field.ny;
        let n = nx * ny;
        let half = self.config.kernel_length / 2;

        // Generate or use provided noise texture
        let noise_buf: Vec<f64>;
        let tex: &[f64] = if let Some(t) = noise {
            t
        } else {
            noise_buf = Self::generate_pseudo_noise(nx, ny);
            &noise_buf
        };

        let mut output = vec![0.0f64; n];

        for j in 0..ny {
            for i in 0..nx {
                let mut sum = 0.0;
                let mut count = 0;

                // Forward integration
                let mut px = i as f64 + 0.5;
                let mut py = j as f64 + 0.5;
                for _step in 0..half {
                    let v = field.interpolate_2d(
                        field.origin.x + px * field.dx,
                        field.origin.y + py * field.dy,
                    );
                    let speed = (v.x * v.x + v.y * v.y).sqrt();
                    if speed < 1e-15 {
                        break;
                    }
                    px += v.x / speed * self.config.step_size;
                    py += v.y / speed * self.config.step_size;
                    let ti = px.clamp(0.0, (nx - 1) as f64) as usize;
                    let tj = py.clamp(0.0, (ny - 1) as f64) as usize;
                    sum += tex[tj * nx + ti];
                    count += 1;
                }

                // Backward integration
                px = i as f64 + 0.5;
                py = j as f64 + 0.5;
                for _step in 0..half {
                    let v = field.interpolate_2d(
                        field.origin.x + px * field.dx,
                        field.origin.y + py * field.dy,
                    );
                    let speed = (v.x * v.x + v.y * v.y).sqrt();
                    if speed < 1e-15 {
                        break;
                    }
                    px -= v.x / speed * self.config.step_size;
                    py -= v.y / speed * self.config.step_size;
                    let ti = px.clamp(0.0, (nx - 1) as f64) as usize;
                    let tj = py.clamp(0.0, (ny - 1) as f64) as usize;
                    sum += tex[tj * nx + ti];
                    count += 1;
                }

                // Add centre pixel
                sum += tex[j * nx + i];
                count += 1;

                output[j * nx + i] = if count > 0 { sum / count as f64 } else { 0.0 };
            }
        }

        // Apply contrast enhancement
        let min_v = output.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_v = output.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if (max_v - min_v) > 1e-15 {
            let range = max_v - min_v;
            for v in &mut output {
                *v = ((*v - min_v) / range * self.config.contrast).clamp(0.0, 1.0);
            }
        }

        output
    }

    /// Generate a deterministic pseudo-noise texture using a simple hash.
    fn generate_pseudo_noise(nx: usize, ny: usize) -> Vec<f64> {
        let n = nx * ny;
        let mut tex = Vec::with_capacity(n);
        for i in 0..n {
            // Simple pseudo-random via xorshift
            let mut v = i as u64 ^ 0xdeadbeef_cafebabe;
            v ^= v << 13;
            v ^= v >> 7;
            v ^= v << 17;
            tex.push((v & 0xFFFF) as f64 / 65535.0);
        }
        tex
    }
}

impl Default for LicRenderer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PathlineTracer
// ---------------------------------------------------------------------------

/// A single pathline through a (possibly time-varying) flow field.
#[derive(Debug, Clone)]
pub struct Pathline {
    /// Sequence of 3D positions along the pathline.
    pub positions: Vec<Vec3>,
    /// Velocity at each position.
    pub velocities: Vec<Vec3>,
    /// Time values corresponding to each position.
    pub times: Vec<f64>,
    /// Whether the integration terminated early (e.g. out of domain).
    pub terminated_early: bool,
}

impl Pathline {
    /// Create a new, empty pathline.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            times: Vec::new(),
            terminated_early: false,
        }
    }

    /// Total arc length of the pathline.
    pub fn arc_length(&self) -> f64 {
        self.positions
            .windows(2)
            .map(|w| {
                let diff = Vec3::new(w[1].x - w[0].x, w[1].y - w[0].y, w[1].z - w[0].z);
                diff.length()
            })
            .sum()
    }

    /// Number of points in the pathline.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Return `true` if the pathline has no points.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

impl Default for Pathline {
    fn default() -> Self {
        Self::new()
    }
}

/// Tracer that integrates particle pathlines through a vector field.
///
/// Uses 4th-order Runge-Kutta (RK4) for forward and backward integration.
pub struct PathlineTracer {
    /// Integration time step.
    pub dt: f64,
    /// Maximum number of integration steps.
    pub max_steps: usize,
    /// Domain bounding box `(min, max)`.
    pub domain_min: Vec3,
    /// Domain bounding box `(min, max)`.
    pub domain_max: Vec3,
}

impl PathlineTracer {
    /// Create a new tracer.
    pub fn new(dt: f64, max_steps: usize, domain_min: Vec3, domain_max: Vec3) -> Self {
        Self {
            dt,
            max_steps,
            domain_min,
            domain_max,
        }
    }

    /// Check whether a position is inside the domain.
    #[inline]
    fn in_domain(&self, p: Vec3) -> bool {
        p.x >= self.domain_min.x
            && p.x <= self.domain_max.x
            && p.y >= self.domain_min.y
            && p.y <= self.domain_max.y
            && p.z >= self.domain_min.z
            && p.z <= self.domain_max.z
    }

    /// Forward-integrate a pathline from `seed` using the given velocity function.
    pub fn trace_forward<F>(&self, seed: Vec3, velocity_fn: &F) -> Pathline
    where
        F: Fn(Vec3, f64) -> Vec3,
    {
        self.trace_impl(seed, velocity_fn, 1.0)
    }

    /// Backward-integrate a pathline from `seed` (reverse time).
    pub fn trace_backward<F>(&self, seed: Vec3, velocity_fn: &F) -> Pathline
    where
        F: Fn(Vec3, f64) -> Vec3,
    {
        self.trace_impl(seed, velocity_fn, -1.0)
    }

    /// Internal RK4 integration.
    fn trace_impl<F>(&self, seed: Vec3, velocity_fn: &F, sign: f64) -> Pathline
    where
        F: Fn(Vec3, f64) -> Vec3,
    {
        let mut pl = Pathline::new();
        let mut pos = seed;
        let mut t = 0.0;
        let h = self.dt * sign;

        pl.positions.push(pos);
        let v0 = velocity_fn(pos, t);
        pl.velocities.push(v0);
        pl.times.push(t);

        for _ in 0..self.max_steps {
            let k1 = velocity_fn(pos, t);
            let k2 = velocity_fn(pos + k1 * (h * 0.5), t + h * 0.5);
            let k3 = velocity_fn(pos + k2 * (h * 0.5), t + h * 0.5);
            let k4 = velocity_fn(pos + k3 * h, t + h);
            let dp = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
            pos = pos + dp;
            t += h;

            if !self.in_domain(pos) {
                pl.terminated_early = true;
                break;
            }
            pl.positions.push(pos);
            let vel = velocity_fn(pos, t);
            pl.velocities.push(vel);
            pl.times.push(t);
        }
        pl
    }

    /// Trace a bundle of pathlines from multiple seed points.
    pub fn trace_bundle<F>(&self, seeds: &[Vec3], velocity_fn: &F) -> Vec<Pathline>
    where
        F: Fn(Vec3, f64) -> Vec3,
    {
        seeds
            .iter()
            .map(|&seed| self.trace_forward(seed, velocity_fn))
            .collect()
    }

    /// Trace pathlines from a uniform seed grid in the XY plane at `z = z_plane`.
    pub fn trace_grid_seeds<F>(
        &self,
        nx_seeds: usize,
        ny_seeds: usize,
        z_plane: f64,
        velocity_fn: &F,
    ) -> Vec<Pathline>
    where
        F: Fn(Vec3, f64) -> Vec3,
    {
        let mut pathlines = Vec::with_capacity(nx_seeds * ny_seeds);
        for j in 0..ny_seeds {
            for i in 0..nx_seeds {
                let px = self.domain_min.x
                    + (i as f64 + 0.5) * (self.domain_max.x - self.domain_min.x) / nx_seeds as f64;
                let py = self.domain_min.y
                    + (j as f64 + 0.5) * (self.domain_max.y - self.domain_min.y) / ny_seeds as f64;
                let seed = Vec3::new(px, py, z_plane);
                pathlines.push(self.trace_forward(seed, velocity_fn));
            }
        }
        pathlines
    }
}

// ---------------------------------------------------------------------------
// VorticityVisualization
// ---------------------------------------------------------------------------

/// Vorticity-based flow identification quantities.
pub struct VorticityVisualization;

impl VorticityVisualization {
    /// Compute vorticity magnitude `|ω|` at each cell of a 3D vector field.
    pub fn vorticity_magnitude(field: &VectorField) -> Vec<f64> {
        field.vorticity_3d().iter().map(|v| v.length()).collect()
    }

    /// Compute the Q-criterion at each cell: `Q = 0.5 * (|Ω|² − |S|²)` where
    /// `Ω` is the vorticity tensor and `S` is the strain-rate tensor.
    ///
    /// Positive Q identifies vortex cores.
    pub fn q_criterion(field: &VectorField) -> Vec<f64> {
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        let n = nx * ny * nz;
        let mut q = vec![0.0f64; n];

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let ip = (i + 1).min(nx - 1);
                    let im = i.saturating_sub(1);
                    let jp = (j + 1).min(ny - 1);
                    let jm = j.saturating_sub(1);
                    let kp = (k + 1).min(nz - 1);
                    let km = k.saturating_sub(1);

                    let two_dx = 2.0 * field.dx * (ip - im).max(1) as f64;
                    let two_dy = 2.0 * field.dy * (jp - jm).max(1) as f64;
                    let two_dz = 2.0 * field.dz * (kp - km).max(1) as f64;

                    // Velocity gradient tensor components
                    let du_dx = (field.get(ip, j, k).x - field.get(im, j, k).x) / two_dx;
                    let du_dy = (field.get(i, jp, k).x - field.get(i, jm, k).x) / two_dy;
                    let du_dz = (field.get(i, j, kp).x - field.get(i, j, km).x) / two_dz;
                    let dv_dx = (field.get(ip, j, k).y - field.get(im, j, k).y) / two_dx;
                    let dv_dy = (field.get(i, jp, k).y - field.get(i, jm, k).y) / two_dy;
                    let dv_dz = (field.get(i, j, kp).y - field.get(i, j, km).y) / two_dz;
                    let dw_dx = (field.get(ip, j, k).z - field.get(im, j, k).z) / two_dx;
                    let dw_dy = (field.get(i, jp, k).z - field.get(i, jm, k).z) / two_dy;
                    let dw_dz = (field.get(i, j, kp).z - field.get(i, j, km).z) / two_dz;

                    // Symmetric strain-rate tensor S
                    let s11 = du_dx;
                    let s22 = dv_dy;
                    let s33 = dw_dz;
                    let s12 = 0.5 * (du_dy + dv_dx);
                    let s13 = 0.5 * (du_dz + dw_dx);
                    let s23 = 0.5 * (dv_dz + dw_dy);

                    // Anti-symmetric vorticity tensor Ω
                    let o12 = 0.5 * (du_dy - dv_dx);
                    let o13 = 0.5 * (du_dz - dw_dx);
                    let o23 = 0.5 * (dv_dz - dw_dy);

                    let s_sq = s11 * s11
                        + s22 * s22
                        + s33 * s33
                        + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                    let o_sq = 2.0 * (o12 * o12 + o13 * o13 + o23 * o23);

                    q[field.index(i, j, k)] = 0.5 * (o_sq - s_sq);
                }
            }
        }
        q
    }

    /// Compute the λ2-criterion eigenvalue at each cell.
    ///
    /// The λ2-criterion identifies vortex cores as regions where the second
    /// eigenvalue of `S² + Ω²` is negative. Here we approximate using the
    /// Q-criterion as `−Q` (a common simplified indicator).
    pub fn lambda2_criterion(field: &VectorField) -> Vec<f64> {
        // Simplified: λ2 ≈ −Q (for incompressible flow)
        let q = Self::q_criterion(field);
        q.into_iter().map(|v| -v).collect()
    }

    /// Compute swirling strength (imaginary part of the complex eigenvalue of ∇u).
    ///
    /// Returns the swirling strength at each grid cell. Non-zero only where the
    /// velocity gradient has complex eigenvalues (i.e., in swirling regions).
    ///
    /// This 2D version uses the 2×2 sub-matrix `\[\[du_dx, du_dy\\], \[dv_dx, dv_dy\]]`.
    pub fn swirling_strength_2d(field: &VectorField) -> Vec<f64> {
        let nx = field.nx;
        let ny = field.ny;
        let n = nx * ny;
        let mut sw = vec![0.0f64; n];

        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1).min(nx - 1);
                let im = i.saturating_sub(1);
                let jp = (j + 1).min(ny - 1);
                let jm = j.saturating_sub(1);
                let two_dx = 2.0 * field.dx * (ip - im).max(1) as f64;
                let two_dy = 2.0 * field.dy * (jp - jm).max(1) as f64;

                let du_dx = (field.get(ip, j, 0).x - field.get(im, j, 0).x) / two_dx;
                let du_dy = (field.get(i, jp, 0).x - field.get(i, jm, 0).x) / two_dy;
                let dv_dx = (field.get(ip, j, 0).y - field.get(im, j, 0).y) / two_dx;
                let dv_dy = (field.get(i, jp, 0).y - field.get(i, jm, 0).y) / two_dy;

                // Eigenvalues of 2x2 matrix [[a,b],[c,d]]
                // λ = (tr ± √(tr²-4det)) / 2
                let tr = du_dx + dv_dy;
                let det = du_dx * dv_dy - du_dy * dv_dx;
                let disc = tr * tr - 4.0 * det;
                if disc < 0.0 {
                    sw[field.index(i, j, 0)] = (-disc).sqrt() * 0.5;
                }
            }
        }
        sw
    }

    /// Identify vortex cores as cells where Q > `threshold`.
    pub fn q_vortex_cores(q_field: &[f64], threshold: f64) -> Vec<usize> {
        q_field
            .iter()
            .enumerate()
            .filter(|&(_, q)| *q > threshold)
            .map(|(i, _)| i)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// PressureField
// ---------------------------------------------------------------------------

/// A scalar pressure field on the same regular grid as a `VectorField`.
#[derive(Debug, Clone)]
pub struct PressureField {
    /// Number of cells along X.
    pub nx: usize,
    /// Number of cells along Y.
    pub ny: usize,
    /// Number of cells along Z.
    pub nz: usize,
    /// Physical spacing in X.
    pub dx: f64,
    /// Physical spacing in Y.
    pub dy: f64,
    /// Physical spacing in Z.
    pub dz: f64,
    /// Grid origin.
    pub origin: Vec3,
    /// Pressure values at each cell.
    pub p: Vec<f64>,
    /// Free-stream pressure `p∞` for Cp calculation.
    pub p_inf: f64,
    /// Free-stream dynamic pressure `q∞ = 0.5 ρ U∞²` for Cp calculation.
    pub q_inf: f64,
}

impl PressureField {
    /// Create a pressure field initialized to zero.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64, origin: Vec3) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            origin,
            p: vec![0.0; nx * ny * nz],
            p_inf: 0.0,
            q_inf: 1.0,
        }
    }

    /// Flat index for `(ix, iy, iz)`.
    #[inline]
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.nx * self.ny + iy * self.nx + ix
    }

    /// Set pressure at cell `(ix, iy, iz)`.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, val: f64) {
        let idx = self.index(ix, iy, iz);
        self.p[idx] = val;
    }

    /// Get pressure at cell `(ix, iy, iz)`.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.p[self.index(ix, iy, iz)]
    }

    /// Compute the pressure coefficient `Cp = (p − p∞) / q∞` at each cell.
    pub fn pressure_coefficient(&self) -> Vec<f64> {
        if self.q_inf.abs() < 1e-15 {
            return vec![0.0; self.p.len()];
        }
        self.p
            .iter()
            .map(|&pi| (pi - self.p_inf) / self.q_inf)
            .collect()
    }

    /// Extract iso-contour indices: cells where `|p\[i\] - iso_value| < tolerance`.
    pub fn iso_contour_2d(&self, iso_value: f64, tolerance: f64) -> Vec<(usize, usize)> {
        let mut cells = Vec::new();
        for j in 0..self.ny {
            for i in 0..self.nx {
                if (self.get(i, j, 0) - iso_value).abs() <= tolerance {
                    cells.push((i, j));
                }
            }
        }
        cells
    }

    /// Extract iso-surface cells in 3D: cells where `|p\[i\] - iso_value| < tolerance`.
    pub fn iso_surface_3d(&self, iso_value: f64, tolerance: f64) -> Vec<(usize, usize, usize)> {
        let mut cells = Vec::new();
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    if (self.get(i, j, k) - iso_value).abs() <= tolerance {
                        cells.push((i, j, k));
                    }
                }
            }
        }
        cells
    }

    /// Bilinear interpolation of pressure at physical position `(px, py)` (2D).
    pub fn interpolate_2d(&self, px: f64, py: f64) -> f64 {
        let fi = (px - self.origin.x) / self.dx;
        let fj = (py - self.origin.y) / self.dy;
        let ti = fi - fi.floor();
        let tj = fj - fj.floor();
        let clamp = |v: isize, max: usize| v.clamp(0, max as isize - 1) as usize;
        let i0 = clamp(fi.floor() as isize, self.nx);
        let i1 = clamp(fi.floor() as isize + 1, self.nx);
        let j0 = clamp(fj.floor() as isize, self.ny);
        let j1 = clamp(fj.floor() as isize + 1, self.ny);
        let p00 = self.get(i0, j0, 0);
        let p10 = self.get(i1, j0, 0);
        let p01 = self.get(i0, j1, 0);
        let p11 = self.get(i1, j1, 0);
        let p0 = p00 * (1.0 - ti) + p10 * ti;
        let p1 = p01 * (1.0 - ti) + p11 * ti;
        p0 * (1.0 - tj) + p1 * tj
    }

    /// Return the minimum and maximum pressure values.
    pub fn min_max(&self) -> (f64, f64) {
        let min = self.p.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = self.p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (min, max)
    }

    /// Return the average pressure over the entire field.
    pub fn mean_pressure(&self) -> f64 {
        if self.p.is_empty() {
            return 0.0;
        }
        self.p.iter().sum::<f64>() / self.p.len() as f64
    }
}

// ---------------------------------------------------------------------------
// FlowStatistics
// ---------------------------------------------------------------------------

/// Time-averaged statistics accumulated from a sequence of instantaneous fields.
///
/// The statistics are maintained in two passes:
/// 1. **Online mean** — computed incrementally using Welford's algorithm.
/// 2. **RMS fluctuations** — computed from the accumulated sum of squares.
#[derive(Debug, Clone)]
pub struct FlowStatistics {
    /// Number of scalar values per snapshot.
    pub size: usize,
    /// Number of snapshots accumulated.
    pub count: usize,
    /// Running mean (Welford's).
    pub mean: Vec<f64>,
    /// Accumulated sum of squared values (for RMS / variance).
    pub sum_sq: Vec<f64>,
    /// Minimum value seen at each position.
    pub min_val: Vec<f64>,
    /// Maximum value seen at each position.
    pub max_val: Vec<f64>,
}

impl FlowStatistics {
    /// Create a zeroed statistics accumulator for fields of length `size`.
    pub fn new(size: usize) -> Self {
        Self {
            size,
            count: 0,
            mean: vec![0.0; size],
            sum_sq: vec![0.0; size],
            min_val: vec![f64::INFINITY; size],
            max_val: vec![f64::NEG_INFINITY; size],
        }
    }

    /// Accumulate one snapshot (scalar field of length `size`).
    ///
    /// Panics in debug mode if `snapshot.len() != self.size`.
    pub fn accumulate(&mut self, snapshot: &[f64]) {
        debug_assert_eq!(snapshot.len(), self.size, "snapshot length mismatch");
        self.count += 1;
        let n = self.count as f64;
        for i in 0..self.size.min(snapshot.len()) {
            let x = snapshot[i];
            let delta = x - self.mean[i];
            self.mean[i] += delta / n;
            self.sum_sq[i] += x * x;
            if x < self.min_val[i] {
                self.min_val[i] = x;
            }
            if x > self.max_val[i] {
                self.max_val[i] = x;
            }
        }
    }

    /// Compute the RMS fluctuation field: `sqrt(mean(x²) − mean(x)²)`.
    pub fn rms_fluctuations(&self) -> Vec<f64> {
        if self.count == 0 {
            return vec![0.0; self.size];
        }
        let n = self.count as f64;
        (0..self.size)
            .map(|i| {
                let mean_sq = self.sum_sq[i] / n;
                let sq_mean = self.mean[i] * self.mean[i];
                (mean_sq - sq_mean).max(0.0).sqrt()
            })
            .collect()
    }

    /// Compute the sample variance field.
    pub fn variance(&self) -> Vec<f64> {
        let rms = self.rms_fluctuations();
        rms.iter().map(|&r| r * r).collect()
    }

    /// Reset all accumulated statistics.
    pub fn reset(&mut self) {
        self.count = 0;
        for i in 0..self.size {
            self.mean[i] = 0.0;
            self.sum_sq[i] = 0.0;
            self.min_val[i] = f64::INFINITY;
            self.max_val[i] = f64::NEG_INFINITY;
        }
    }

    /// Compute the 1D energy spectrum from a scalar field snapshot using DFT.
    ///
    /// Returns a vector of length `n/2 + 1` where each element is the power
    /// spectral density at the corresponding wavenumber.
    pub fn energy_spectrum_1d(signal: &[f64]) -> Vec<f64> {
        let n = signal.len();
        if n == 0 {
            return Vec::new();
        }
        let out_len = n / 2 + 1;
        let mut psd = vec![0.0f64; out_len];

        // DFT (O(n²), sufficient for moderate n in simulation use)
        for k in 0..out_len {
            let mut re = 0.0f64;
            let mut im = 0.0f64;
            for j in 0..n {
                let angle = -2.0 * PI * k as f64 * j as f64 / n as f64;
                re += signal[j] * angle.cos();
                im += signal[j] * angle.sin();
            }
            psd[k] = (re * re + im * im) / n as f64;
        }

        // Double non-DC, non-Nyquist components (one-sided spectrum)
        for k in 1..out_len.saturating_sub(1) {
            psd[k] *= 2.0;
        }
        psd
    }

    /// Compute the 2D energy spectrum, binned radially, from a 2D scalar field.
    ///
    /// `nx` and `ny` are the dimensions; the field is stored row-major.
    /// Returns a vec of length `min(nx,ny)/2 + 1` with the azimuthally averaged PSD.
    pub fn energy_spectrum_2d(field: &[f64], nx: usize, ny: usize) -> Vec<f64> {
        if field.is_empty() || nx == 0 || ny == 0 {
            return Vec::new();
        }
        let kmax = (nx.min(ny) / 2) + 1;
        let mut bins = vec![0.0f64; kmax];
        let mut counts = vec![0usize; kmax];

        for j in 0..ny {
            for i in 0..nx {
                // Wavenumber indices (centred)
                let ki = if i <= nx / 2 {
                    i as f64
                } else {
                    (i as isize - nx as isize) as f64
                };
                let kj = if j <= ny / 2 {
                    j as f64
                } else {
                    (j as isize - ny as isize) as f64
                };
                let kr = (ki * ki + kj * kj).sqrt();
                let bin = kr.round() as usize;
                if bin < kmax {
                    // DFT at (i,j) using brute-force sum
                    let mut re = 0.0f64;
                    let mut im = 0.0f64;
                    for jj in 0..ny {
                        for ii in 0..nx {
                            let angle = -2.0
                                * PI
                                * (ki * ii as f64 / nx as f64 + kj * jj as f64 / ny as f64);
                            re += field[jj * nx + ii] * angle.cos();
                            im += field[jj * nx + ii] * angle.sin();
                        }
                    }
                    bins[bin] += (re * re + im * im) / (nx * ny) as f64;
                    counts[bin] += 1;
                }
            }
        }

        // Average by bin count
        for (b, c) in bins.iter_mut().zip(counts.iter()) {
            if *c > 0 {
                *b /= *c as f64;
            }
        }
        bins
    }

    /// Compute the turbulent kinetic energy (TKE) at each point from three
    /// velocity component RMS fluctuations.
    pub fn turbulent_kinetic_energy(urms: &[f64], vrms: &[f64], wrms: &[f64]) -> Vec<f64> {
        let n = urms.len().min(vrms.len()).min(wrms.len());
        (0..n)
            .map(|i| 0.5 * (urms[i] * urms[i] + vrms[i] * vrms[i] + wrms[i] * wrms[i]))
            .collect()
    }

    /// Compute the integral length scale from a 1D auto-correlation.
    ///
    /// Integrates the auto-correlation from lag 0 until it first crosses zero.
    pub fn integral_length_scale(signal: &[f64], dx: f64) -> f64 {
        let n = signal.len();
        if n < 2 {
            return 0.0;
        }
        let mean = signal.iter().sum::<f64>() / n as f64;
        let var: f64 = signal.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n as f64;
        if var < 1e-15 {
            return 0.0;
        }
        let mut integral = 0.0;
        for lag in 1..n {
            let mut corr = 0.0;
            let cnt = n - lag;
            for i in 0..cnt {
                corr += (signal[i] - mean) * (signal[i + lag] - mean);
            }
            corr /= cnt as f64 * var;
            if corr <= 0.0 {
                break;
            }
            integral += corr * dx;
        }
        integral
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Vec2 ----

    #[test]
    fn test_vec2_length() {
        let v = Vec2::new(3.0, 4.0);
        assert!((v.length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec2_normalize() {
        let v = Vec2::new(3.0, 4.0).normalize();
        assert!((v.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec2_dot() {
        let a = Vec2::new(1.0, 0.0);
        let b = Vec2::new(0.0, 1.0);
        assert!(a.dot(b).abs() < 1e-10);
    }

    #[test]
    fn test_vec2_cross() {
        let a = Vec2::new(1.0, 0.0);
        let b = Vec2::new(0.0, 1.0);
        assert!((a.cross(b) - 1.0).abs() < 1e-10);
    }

    // ---- Vec3 ----

    #[test]
    fn test_vec3_length() {
        let v = Vec3::new(1.0, 2.0, 2.0);
        assert!((v.length() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_normalize_unit() {
        let v = Vec3::new(0.0, 0.0, 5.0).normalize();
        assert!((v.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_cross_orthogonal() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        let z = x.cross(y);
        assert!((z.z - 1.0).abs() < 1e-10);
        assert!(z.x.abs() < 1e-10 && z.y.abs() < 1e-10);
    }

    #[test]
    fn test_vec3_dot_perpendicular() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        assert!(a.dot(b).abs() < 1e-10);
    }

    // ---- VectorField ----

    #[test]
    fn test_vector_field_2d_index() {
        let field = VectorField::new_2d(4, 3, 1.0, 1.0, Vec2::zero());
        assert_eq!(field.index(0, 0, 0), 0);
        assert_eq!(field.index(3, 0, 0), 3);
        assert_eq!(field.index(0, 1, 0), 4);
        assert_eq!(field.index(3, 2, 0), 11);
    }

    #[test]
    fn test_vector_field_set_get() {
        let mut field = VectorField::new_2d(5, 5, 1.0, 1.0, Vec2::zero());
        field.set(2, 3, 0, Vec3::new(1.0, 2.0, 0.0));
        let v = field.get(2, 3, 0);
        assert!((v.x - 1.0).abs() < 1e-10);
        assert!((v.y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_field_fill_uniform() {
        let mut field = VectorField::new_2d(4, 4, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(2.0, 0.0, 0.0));
        assert!((field.max_velocity() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_field_interpolate_2d_uniform() {
        let mut field = VectorField::new_2d(8, 8, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(1.0, 0.5, 0.0));
        let v = field.interpolate_2d(3.5, 3.5);
        assert!((v.x - 1.0).abs() < 1e-10);
        assert!((v.y - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_vector_field_divergence_uniform() {
        let mut field = VectorField::new_2d(8, 8, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(1.0, 1.0, 0.0));
        let div = field.divergence_2d();
        // Divergence of constant field = 0
        for &d in &div {
            assert!(
                d.abs() < 1e-10,
                "divergence of uniform field should be 0, got {}",
                d
            );
        }
    }

    #[test]
    fn test_vector_field_curl_2d_vortex() {
        // A vortex has non-zero curl everywhere (except at the singularity)
        let mut field = VectorField::new_2d(10, 10, 1.0, 1.0, Vec2::new(5.0, 5.0));
        field.fill_vortex_2d(Vec2::new(5.5, 5.5), 1.0);
        let curl = field.curl_2d();
        // At least some cells should have non-zero curl
        let nonzero = curl.iter().any(|&c| c.abs() > 1e-10);
        assert!(nonzero, "vortex field should have non-zero curl");
    }

    #[test]
    fn test_vector_field_3d_index() {
        let field = VectorField::new_3d(3, 4, 5, 1.0, 1.0, 1.0, Vec3::zero());
        assert_eq!(field.index(0, 0, 0), 0);
        assert_eq!(field.index(2, 3, 4), 2 + 3 * 3 + 4 * 3 * 4);
    }

    #[test]
    fn test_vector_field_interpolate_3d_uniform() {
        let mut field = VectorField::new_3d(6, 6, 6, 1.0, 1.0, 1.0, Vec3::zero());
        field.fill_uniform(Vec3::new(3.0, 2.0, 1.0));
        let v = field.interpolate_3d(2.5, 2.5, 2.5);
        assert!((v.x - 3.0).abs() < 1e-10);
        assert!((v.y - 2.0).abs() < 1e-10);
        assert!((v.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_field_vorticity_3d_uniform() {
        let mut field = VectorField::new_3d(5, 5, 5, 1.0, 1.0, 1.0, Vec3::zero());
        field.fill_uniform(Vec3::new(1.0, 1.0, 1.0));
        let omega = field.vorticity_3d();
        for v in &omega {
            assert!(v.length() < 1e-10, "vorticity of uniform field should be 0");
        }
    }

    // ---- LicRenderer ----

    #[test]
    fn test_lic_output_size() {
        let mut field = VectorField::new_2d(16, 16, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(1.0, 0.0, 0.0));
        let lic = LicRenderer::new();
        let output = lic.render(&field, None);
        assert_eq!(output.len(), 16 * 16);
    }

    #[test]
    fn test_lic_output_range() {
        let mut field = VectorField::new_2d(12, 12, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(0.5, 0.5, 0.0));
        let lic = LicRenderer::new();
        let output = lic.render(&field, None);
        for &v in &output {
            assert!(
                (0.0..=1.0).contains(&v),
                "LIC output must be in [0,1], got {}",
                v
            );
        }
    }

    #[test]
    fn test_lic_custom_noise() {
        let mut field = VectorField::new_2d(8, 8, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(1.0, 0.0, 0.0));
        let noise = vec![0.5f64; 64];
        let lic = LicRenderer::new();
        let output = lic.render(&field, Some(&noise));
        assert_eq!(output.len(), 64);
    }

    #[test]
    fn test_lic_config_default() {
        let cfg = LicConfig::default();
        assert_eq!(cfg.kernel_length, 20);
        assert!(cfg.step_size > 0.0);
    }

    // ---- PathlineTracer ----

    #[test]
    fn test_pathline_forward_uniform() {
        let tracer = PathlineTracer::new(
            0.1,
            10,
            Vec3::new(-5.0, -5.0, -5.0),
            Vec3::new(5.0, 5.0, 5.0),
        );
        let seed = Vec3::new(0.0, 0.0, 0.0);
        let vel_fn = |_p: Vec3, _t: f64| Vec3::new(1.0, 0.0, 0.0);
        let pl = tracer.trace_forward(seed, &vel_fn);
        assert!(pl.len() > 1, "pathline should have more than 1 point");
        // Should move in +x direction
        let last = pl.positions.last().unwrap();
        assert!(last.x > 0.0, "last x should be positive");
    }

    #[test]
    fn test_pathline_backward_uniform() {
        let tracer = PathlineTracer::new(
            0.1,
            10,
            Vec3::new(-10.0, -10.0, -10.0),
            Vec3::new(10.0, 10.0, 10.0),
        );
        let seed = Vec3::new(0.0, 0.0, 0.0);
        let vel_fn = |_p: Vec3, _t: f64| Vec3::new(1.0, 0.0, 0.0);
        let pl = tracer.trace_backward(seed, &vel_fn);
        let last = pl.positions.last().unwrap();
        assert!(last.x < 0.0, "backward trace should go in -x direction");
    }

    #[test]
    fn test_pathline_terminates_at_boundary() {
        let tracer = PathlineTracer::new(
            0.5,
            100,
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        );
        let seed = Vec3::new(0.0, 0.0, 0.0);
        let vel_fn = |_p: Vec3, _t: f64| Vec3::new(1.0, 0.0, 0.0);
        let pl = tracer.trace_forward(seed, &vel_fn);
        assert!(pl.terminated_early, "should terminate at domain boundary");
    }

    #[test]
    fn test_pathline_arc_length() {
        let mut pl = Pathline::new();
        pl.positions.push(Vec3::new(0.0, 0.0, 0.0));
        pl.positions.push(Vec3::new(1.0, 0.0, 0.0));
        pl.positions.push(Vec3::new(2.0, 0.0, 0.0));
        assert!((pl.arc_length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_pathline_bundle_count() {
        let tracer = PathlineTracer::new(
            0.1,
            5,
            Vec3::new(-10.0, -10.0, -10.0),
            Vec3::new(10.0, 10.0, 10.0),
        );
        let seeds = vec![
            Vec3::zero(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let vel_fn = |_p: Vec3, _t: f64| Vec3::new(0.5, 0.0, 0.0);
        let bundle = tracer.trace_bundle(&seeds, &vel_fn);
        assert_eq!(bundle.len(), 3);
    }

    #[test]
    fn test_pathline_grid_seeds_count() {
        let tracer = PathlineTracer::new(
            0.1,
            3,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(10.0, 10.0, 1.0),
        );
        let vel_fn = |_p: Vec3, _t: f64| Vec3::new(0.1, 0.0, 0.0);
        let lines = tracer.trace_grid_seeds(3, 3, 0.0, &vel_fn);
        assert_eq!(lines.len(), 9);
    }

    // ---- VorticityVisualization ----

    #[test]
    fn test_vorticity_magnitude_uniform_zero() {
        let mut field = VectorField::new_3d(4, 4, 4, 1.0, 1.0, 1.0, Vec3::zero());
        field.fill_uniform(Vec3::new(1.0, 1.0, 1.0));
        let omega_mag = VorticityVisualization::vorticity_magnitude(&field);
        for &v in &omega_mag {
            assert!(
                v < 1e-10,
                "vorticity magnitude of uniform field = 0, got {}",
                v
            );
        }
    }

    #[test]
    fn test_q_criterion_uniform_zero() {
        let mut field = VectorField::new_3d(4, 4, 4, 1.0, 1.0, 1.0, Vec3::zero());
        field.fill_uniform(Vec3::new(1.0, 0.0, 0.0));
        let q = VorticityVisualization::q_criterion(&field);
        for &v in &q {
            assert!(v.abs() < 1e-10, "Q of uniform field = 0, got {}", v);
        }
    }

    #[test]
    fn test_q_vortex_cores_threshold() {
        let q_field = vec![-1.0, 0.5, 2.0, 3.0, -0.5];
        let cores = VorticityVisualization::q_vortex_cores(&q_field, 1.0);
        assert_eq!(cores.len(), 2, "should find cells with Q > 1.0");
    }

    #[test]
    fn test_lambda2_criterion_sign_flip() {
        let mut field = VectorField::new_3d(3, 3, 3, 1.0, 1.0, 1.0, Vec3::zero());
        field.fill_uniform(Vec3::zero());
        let q = VorticityVisualization::q_criterion(&field);
        let l2 = VorticityVisualization::lambda2_criterion(&field);
        for (qi, li) in q.iter().zip(l2.iter()) {
            assert!((*qi + *li).abs() < 1e-10, "lambda2 = -Q");
        }
    }

    #[test]
    fn test_swirling_strength_uniform_zero() {
        let mut field = VectorField::new_2d(6, 6, 1.0, 1.0, Vec2::zero());
        field.fill_uniform(Vec3::new(1.0, 0.0, 0.0));
        let sw = VorticityVisualization::swirling_strength_2d(&field);
        for &v in &sw {
            assert!(
                v < 1e-10,
                "swirling strength of uniform field = 0, got {}",
                v
            );
        }
    }

    #[test]
    fn test_swirling_strength_vortex_nonzero() {
        let mut field = VectorField::new_2d(10, 10, 1.0, 1.0, Vec2::new(5.0, 5.0));
        field.fill_vortex_2d(Vec2::new(5.5, 5.5), 2.0);
        let sw = VorticityVisualization::swirling_strength_2d(&field);
        let nonzero = sw.iter().any(|&v| v > 1e-10);
        assert!(nonzero, "vortex should have non-zero swirling strength");
    }

    // ---- PressureField ----

    #[test]
    fn test_pressure_field_set_get() {
        let mut pf = PressureField::new(4, 4, 1, 1.0, 1.0, 1.0, Vec3::zero());
        pf.set(2, 3, 0, 101325.0);
        assert!((pf.get(2, 3, 0) - 101325.0).abs() < 1e-6);
    }

    #[test]
    fn test_pressure_coefficient() {
        let mut pf = PressureField::new(3, 3, 1, 1.0, 1.0, 1.0, Vec3::zero());
        pf.p_inf = 100.0;
        pf.q_inf = 50.0;
        pf.set(1, 1, 0, 200.0);
        let cp = pf.pressure_coefficient();
        let idx = pf.index(1, 1, 0);
        assert!(
            (cp[idx] - 2.0).abs() < 1e-10,
            "Cp = (200-100)/50 = 2, got {}",
            cp[idx]
        );
    }

    #[test]
    fn test_pressure_field_min_max() {
        let mut pf = PressureField::new(3, 3, 1, 1.0, 1.0, 1.0, Vec3::zero());
        for j in 0..3 {
            for i in 0..3 {
                pf.set(i, j, 0, (i + j * 3) as f64);
            }
        }
        let (mn, mx) = pf.min_max();
        assert!((mn - 0.0).abs() < 1e-10);
        assert!((mx - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_pressure_field_mean() {
        let mut pf = PressureField::new(2, 2, 1, 1.0, 1.0, 1.0, Vec3::zero());
        pf.set(0, 0, 0, 1.0);
        pf.set(1, 0, 0, 2.0);
        pf.set(0, 1, 0, 3.0);
        pf.set(1, 1, 0, 4.0);
        assert!((pf.mean_pressure() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_pressure_field_iso_contour_2d() {
        let mut pf = PressureField::new(5, 5, 1, 1.0, 1.0, 1.0, Vec3::zero());
        pf.set(2, 2, 0, 100.0);
        pf.set(3, 3, 0, 100.0);
        let cells = pf.iso_contour_2d(100.0, 0.01);
        assert_eq!(cells.len(), 2);
    }

    #[test]
    fn test_pressure_field_iso_surface_3d() {
        let mut pf = PressureField::new(3, 3, 3, 1.0, 1.0, 1.0, Vec3::zero());
        pf.set(1, 1, 1, 42.0);
        let cells = pf.iso_surface_3d(42.0, 0.1);
        assert_eq!(cells.len(), 1);
    }

    #[test]
    fn test_pressure_field_interpolate_2d_uniform() {
        let mut pf = PressureField::new(6, 6, 1, 1.0, 1.0, 1.0, Vec3::zero());
        for v in &mut pf.p {
            *v = 5.0;
        }
        let val = pf.interpolate_2d(2.5, 2.5);
        assert!((val - 5.0).abs() < 1e-10);
    }

    // ---- FlowStatistics ----

    #[test]
    fn test_flow_stats_mean() {
        let mut stats = FlowStatistics::new(3);
        stats.accumulate(&[1.0, 2.0, 3.0]);
        stats.accumulate(&[3.0, 4.0, 5.0]);
        assert!((stats.mean[0] - 2.0).abs() < 1e-10);
        assert!((stats.mean[1] - 3.0).abs() < 1e-10);
        assert!((stats.mean[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_flow_stats_rms_constant_zero_fluctuation() {
        let mut stats = FlowStatistics::new(2);
        stats.accumulate(&[5.0, 3.0]);
        stats.accumulate(&[5.0, 3.0]);
        let rms = stats.rms_fluctuations();
        assert!(rms[0] < 1e-10, "RMS of constant signal should be 0");
        assert!(rms[1] < 1e-10, "RMS of constant signal should be 0");
    }

    #[test]
    fn test_flow_stats_rms_known_value() {
        let mut stats = FlowStatistics::new(1);
        // Values: -1, +1 → mean=0, var=1, rms=1
        stats.accumulate(&[-1.0]);
        stats.accumulate(&[1.0]);
        let rms = stats.rms_fluctuations();
        assert!(
            (rms[0] - 1.0).abs() < 1e-10,
            "RMS should be 1.0, got {}",
            rms[0]
        );
    }

    #[test]
    fn test_flow_stats_reset() {
        let mut stats = FlowStatistics::new(2);
        stats.accumulate(&[1.0, 2.0]);
        stats.reset();
        assert_eq!(stats.count, 0);
        assert!(stats.mean[0].abs() < 1e-10);
    }

    #[test]
    fn test_flow_stats_min_max() {
        let mut stats = FlowStatistics::new(1);
        stats.accumulate(&[3.0]);
        stats.accumulate(&[7.0]);
        stats.accumulate(&[1.0]);
        assert!((stats.min_val[0] - 1.0).abs() < 1e-10);
        assert!((stats.max_val[0] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_energy_spectrum_1d_length() {
        let n = 32;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let psd = FlowStatistics::energy_spectrum_1d(&signal);
        assert_eq!(psd.len(), n / 2 + 1);
    }

    #[test]
    fn test_energy_spectrum_1d_dc_zero_mean() {
        // Zero-mean pure sine: DC component should be small
        let n = 64;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 4.0 * i as f64 / n as f64).sin())
            .collect();
        let psd = FlowStatistics::energy_spectrum_1d(&signal);
        assert!(psd[0] < 1.0, "DC power of zero-mean signal should be small");
    }

    #[test]
    fn test_turbulent_kinetic_energy() {
        let urms = vec![1.0, 2.0];
        let vrms = vec![1.0, 1.0];
        let wrms = vec![1.0, 1.0];
        let tke = FlowStatistics::turbulent_kinetic_energy(&urms, &vrms, &wrms);
        // TKE[0] = 0.5*(1+1+1) = 1.5
        assert!((tke[0] - 1.5).abs() < 1e-10);
        // TKE[1] = 0.5*(4+1+1) = 3.0
        assert!((tke[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_integral_length_scale_constant() {
        // Constant signal has zero variance → integral scale = 0
        let signal = vec![5.0f64; 50];
        let l = FlowStatistics::integral_length_scale(&signal, 0.1);
        assert!(l.abs() < 1e-10);
    }

    #[test]
    fn test_energy_spectrum_2d_length() {
        let nx = 8;
        let ny = 8;
        let field = vec![1.0f64; nx * ny];
        let psd = FlowStatistics::energy_spectrum_2d(&field, nx, ny);
        assert_eq!(psd.len(), nx.min(ny) / 2 + 1);
    }

    #[test]
    fn test_pathline_rk4_circular() {
        // Circular field v=(−y, x, 0): particle should orbit the origin
        let tracer = PathlineTracer::new(
            0.01,
            628,
            Vec3::new(-3.0, -3.0, -1.0),
            Vec3::new(3.0, 3.0, 1.0),
        );
        let seed = Vec3::new(1.0, 0.0, 0.0);
        let vel_fn = |p: Vec3, _t: f64| Vec3::new(-p.y, p.x, 0.0);
        let pl = tracer.trace_forward(seed, &vel_fn);
        assert!(pl.len() > 10, "should trace many steps");
        // After a full revolution (~2π/dt steps), return close to origin distance 1
        let last = pl.positions.last().unwrap();
        let dist = (last.x * last.x + last.y * last.y).sqrt();
        assert!(
            (dist - 1.0).abs() < 0.1,
            "orbit radius should remain ≈ 1, got {}",
            dist
        );
    }

    #[test]
    fn test_flow_stats_variance() {
        let mut stats = FlowStatistics::new(1);
        stats.accumulate(&[2.0]);
        stats.accumulate(&[4.0]);
        let var = stats.variance();
        // var = mean(x²) - mean(x)² = 10 - 9 = 1
        assert!(
            (var[0] - 1.0).abs() < 1e-10,
            "variance = 1.0, got {}",
            var[0]
        );
    }
}
