#![allow(clippy::needless_range_loop, clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Flow field visualization utilities.
//!
//! Provides a 3-D Eulerian velocity field, streamline integration (RK4),
//! pathline integration (Euler), streakline construction, vorticity,
//! divergence, Q-criterion, lambda-2 vortex criterion, vorticity magnitude,
//! 2-D Line Integral Convolution (LIC), uniform seed point generation,
//! critical-point extraction, separatrix tracing, and vector field arrow
//! generation.

#![allow(dead_code)]

use rand::RngExt as _;
// ─────────────────────────────────────────────────────────────────────────────
// FlowField3D
// ─────────────────────────────────────────────────────────────────────────────

/// A 3-D Eulerian velocity field stored on a uniform Cartesian grid.
///
/// Velocity components are stored in flat row-major order:
/// index = `iz * ny * nx + iy * nx + ix`.
#[derive(Debug, Clone)]
pub struct FlowField3D {
    /// X-component of velocity at each grid point.
    pub u: Vec<f64>,
    /// Y-component of velocity at each grid point.
    pub v: Vec<f64>,
    /// W-component (Z) of velocity at each grid point.
    pub w: Vec<f64>,
    /// Number of grid points along X.
    pub nx: usize,
    /// Number of grid points along Y.
    pub ny: usize,
    /// Number of grid points along Z.
    pub nz: usize,
    /// Uniform grid spacing (isotropic).
    pub dx: f64,
}

impl FlowField3D {
    /// Create a zero-velocity field of the given dimensions.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        let n = nx * ny * nz;
        Self {
            u: vec![0.0; n],
            v: vec![0.0; n],
            w: vec![0.0; n],
            nx,
            ny,
            nz,
            dx,
        }
    }

    /// Total number of grid points.
    pub fn len(&self) -> usize {
        self.u.len()
    }

    /// Return `true` if the field has no grid points.
    pub fn is_empty(&self) -> bool {
        self.u.is_empty()
    }

    /// Flat index for grid coordinates `(ix, iy, iz)`, clamped to boundaries.
    pub fn idx_clamped(&self, ix: i64, iy: i64, iz: i64) -> usize {
        let ix = ix.clamp(0, self.nx as i64 - 1) as usize;
        let iy = iy.clamp(0, self.ny as i64 - 1) as usize;
        let iz = iz.clamp(0, self.nz as i64 - 1) as usize;
        iz * self.ny * self.nx + iy * self.nx + ix
    }

    /// Get velocity `(u, v, w)` at grid coordinates `(ix, iy, iz)`, clamped.
    pub fn vel_at(&self, ix: i64, iy: i64, iz: i64) -> [f64; 3] {
        let i = self.idx_clamped(ix, iy, iz);
        [self.u[i], self.v[i], self.w[i]]
    }

    /// Set velocity at grid coordinates `(ix, iy, iz)`, clamped.
    pub fn set_vel(&mut self, ix: i64, iy: i64, iz: i64, vel: [f64; 3]) {
        let i = self.idx_clamped(ix, iy, iz);
        self.u[i] = vel[0];
        self.v[i] = vel[1];
        self.w[i] = vel[2];
    }

    /// Compute the velocity magnitude at a grid point.
    pub fn speed_at(&self, ix: i64, iy: i64, iz: i64) -> f64 {
        let v = self.vel_at(ix, iy, iz);
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector arithmetic helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors: `a - b`.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Magnitude of a 3-vector.
#[inline]
fn mag3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalize a 3-vector (returns zero vector if magnitude is tiny).
#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let m = mag3(a);
    if m < 1e-30 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / m)
    }
}

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// interpolate_velocity
// ─────────────────────────────────────────────────────────────────────────────

/// Interpolate the velocity at a continuous position `pos` in grid units
/// using trilinear interpolation.
#[allow(clippy::too_many_arguments)]
pub fn interpolate_velocity(field: &FlowField3D, pos: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = pos;
    let x0 = x.floor() as i64;
    let y0 = y.floor() as i64;
    let z0 = z.floor() as i64;
    let tx = x - x0 as f64;
    let ty = y - y0 as f64;
    let tz = z - z0 as f64;

    let lerp_vel = |v000: [f64; 3],
                    v100: [f64; 3],
                    v010: [f64; 3],
                    v110: [f64; 3],
                    v001: [f64; 3],
                    v101: [f64; 3],
                    v011: [f64; 3],
                    v111: [f64; 3]|
     -> [f64; 3] {
        let mut out = [0.0; 3];
        for k in 0..3 {
            let c00 = v000[k] * (1.0 - tx) + v100[k] * tx;
            let c10 = v010[k] * (1.0 - tx) + v110[k] * tx;
            let c01 = v001[k] * (1.0 - tx) + v101[k] * tx;
            let c11 = v011[k] * (1.0 - tx) + v111[k] * tx;
            let c0 = c00 * (1.0 - ty) + c10 * ty;
            let c1 = c01 * (1.0 - ty) + c11 * ty;
            out[k] = c0 * (1.0 - tz) + c1 * tz;
        }
        out
    };

    lerp_vel(
        field.vel_at(x0, y0, z0),
        field.vel_at(x0 + 1, y0, z0),
        field.vel_at(x0, y0 + 1, z0),
        field.vel_at(x0 + 1, y0 + 1, z0),
        field.vel_at(x0, y0, z0 + 1),
        field.vel_at(x0 + 1, y0, z0 + 1),
        field.vel_at(x0, y0 + 1, z0 + 1),
        field.vel_at(x0 + 1, y0 + 1, z0 + 1),
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Streamline integration (RK4)
// ─────────────────────────────────────────────────────────────────────────────

/// Trace a streamline from `start` using 4th-order Runge-Kutta integration.
///
/// Positions are in grid units. Returns `n_steps + 1` points (seed + steps).
pub fn streamline_rk4(
    field: &FlowField3D,
    start: [f64; 3],
    dt: f64,
    n_steps: usize,
) -> Vec<[f64; 3]> {
    let mut path = Vec::with_capacity(n_steps + 1);
    let mut pos = start;
    path.push(pos);
    for _ in 0..n_steps {
        let k1 = interpolate_velocity(field, pos);
        let k2 = interpolate_velocity(field, add3(pos, scale3(k1, dt * 0.5)));
        let k3 = interpolate_velocity(field, add3(pos, scale3(k2, dt * 0.5)));
        let k4 = interpolate_velocity(field, add3(pos, scale3(k3, dt)));
        let vel = [
            (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
            (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
            (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0,
        ];
        pos = add3(pos, scale3(vel, dt));
        path.push(pos);
    }
    path
}

/// Trace a streamline backwards (against the velocity direction) using RK4.
///
/// Returns `n_steps + 1` points starting from `start`.
pub fn streamline_rk4_backward(
    field: &FlowField3D,
    start: [f64; 3],
    dt: f64,
    n_steps: usize,
) -> Vec<[f64; 3]> {
    streamline_rk4(field, start, -dt, n_steps)
}

// ─────────────────────────────────────────────────────────────────────────────
// Pathline integration (Euler)
// ─────────────────────────────────────────────────────────────────────────────

/// Trace a pathline from `start` using forward Euler integration.
///
/// Returns `n_steps + 1` points (seed + steps).
pub fn pathline_euler(
    field: &FlowField3D,
    start: [f64; 3],
    dt: f64,
    n_steps: usize,
) -> Vec<[f64; 3]> {
    let mut path = Vec::with_capacity(n_steps + 1);
    let mut pos = start;
    path.push(pos);
    for _ in 0..n_steps {
        let vel = interpolate_velocity(field, pos);
        pos = add3(pos, scale3(vel, dt));
        path.push(pos);
    }
    path
}

// ─────────────────────────────────────────────────────────────────────────────
// Streakline
// ─────────────────────────────────────────────────────────────────────────────

/// A streakline records the positions of particles released from a fixed
/// injection point at successive time steps.
#[derive(Debug, Clone)]
pub struct Streakline {
    /// Injection point (in grid units).
    pub injection_point: [f64; 3],
    /// Current positions of all released particles, oldest first.
    pub particles: Vec<[f64; 3]>,
}

impl Streakline {
    /// Create a new streakline from an injection point.
    pub fn new(injection_point: [f64; 3]) -> Self {
        Self {
            injection_point,
            particles: Vec::new(),
        }
    }

    /// Advance all existing particles by one Euler step and inject a new
    /// particle at the injection point.
    pub fn advance(&mut self, field: &FlowField3D, dt: f64) {
        // Advect existing particles
        for p in &mut self.particles {
            let vel = interpolate_velocity(field, *p);
            *p = add3(*p, scale3(vel, dt));
        }
        // Inject new particle
        self.particles.push(self.injection_point);
    }

    /// Number of particles currently in the streakline.
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// Return `true` if no particles have been released yet.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vorticity (curl of velocity)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the vorticity field omega = curl(v) using central finite differences.
///
/// Returns a new `FlowField3D` where `(u, v, w)` are the `(omega_x, omega_y, omega_z)`
/// components respectively. Boundary values use one-sided differences.
pub fn vorticity(field: &FlowField3D) -> FlowField3D {
    let nx = field.nx;
    let ny = field.ny;
    let nz = field.nz;
    let h = field.dx;
    let mut out = FlowField3D::new(nx, ny, nz, h);

    for iz in 0..nz as i64 {
        for iy in 0..ny as i64 {
            for ix in 0..nx as i64 {
                let dw_dy =
                    (field.vel_at(ix, iy + 1, iz)[2] - field.vel_at(ix, iy - 1, iz)[2]) / (2.0 * h);
                let dv_dz =
                    (field.vel_at(ix, iy, iz + 1)[1] - field.vel_at(ix, iy, iz - 1)[1]) / (2.0 * h);
                let du_dz =
                    (field.vel_at(ix, iy, iz + 1)[0] - field.vel_at(ix, iy, iz - 1)[0]) / (2.0 * h);
                let dw_dx =
                    (field.vel_at(ix + 1, iy, iz)[2] - field.vel_at(ix - 1, iy, iz)[2]) / (2.0 * h);
                let dv_dx =
                    (field.vel_at(ix + 1, iy, iz)[1] - field.vel_at(ix - 1, iy, iz)[1]) / (2.0 * h);
                let du_dy =
                    (field.vel_at(ix, iy + 1, iz)[0] - field.vel_at(ix, iy - 1, iz)[0]) / (2.0 * h);

                let idx = out.idx_clamped(ix, iy, iz);
                out.u[idx] = dw_dy - dv_dz;
                out.v[idx] = du_dz - dw_dx;
                out.w[idx] = dv_dx - du_dy;
            }
        }
    }
    out
}

/// Compute the vorticity magnitude field `|omega|` at each grid point.
///
/// Returns a flat scalar array of the same length as `field.u`.
pub fn vorticity_magnitude(field: &FlowField3D) -> Vec<f64> {
    let vor = vorticity(field);
    let n = vor.len();
    let mut mag = vec![0.0_f64; n];
    for i in 0..n {
        mag[i] = (vor.u[i] * vor.u[i] + vor.v[i] * vor.v[i] + vor.w[i] * vor.w[i]).sqrt();
    }
    mag
}

// ─────────────────────────────────────────────────────────────────────────────
// Divergence
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the divergence `div(v) = du/dx + dv/dy + dw/dz` using central
/// finite differences.
///
/// Returns a flat scalar array of the same length as `field.u`.
pub fn divergence(field: &FlowField3D) -> Vec<f64> {
    let nx = field.nx;
    let ny = field.ny;
    let nz = field.nz;
    let h = field.dx;
    let mut out = vec![0.0_f64; nx * ny * nz];

    for iz in 0..nz as i64 {
        for iy in 0..ny as i64 {
            for ix in 0..nx as i64 {
                let du_dx =
                    (field.vel_at(ix + 1, iy, iz)[0] - field.vel_at(ix - 1, iy, iz)[0]) / (2.0 * h);
                let dv_dy =
                    (field.vel_at(ix, iy + 1, iz)[1] - field.vel_at(ix, iy - 1, iz)[1]) / (2.0 * h);
                let dw_dz =
                    (field.vel_at(ix, iy, iz + 1)[2] - field.vel_at(ix, iy, iz - 1)[2]) / (2.0 * h);
                let idx = field.idx_clamped(ix, iy, iz);
                out[idx] = du_dx + dv_dy + dw_dz;
            }
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Velocity gradient tensor helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the 3x3 velocity gradient tensor `du_i/dx_j` at grid index `(ix, iy, iz)`.
/// Returns a flat 9-element array in row-major order.
fn velocity_gradient(field: &FlowField3D, ix: i64, iy: i64, iz: i64) -> [f64; 9] {
    let h = field.dx;
    let mut g = [0.0_f64; 9];
    g[0] = (field.vel_at(ix + 1, iy, iz)[0] - field.vel_at(ix - 1, iy, iz)[0]) / (2.0 * h);
    g[1] = (field.vel_at(ix, iy + 1, iz)[0] - field.vel_at(ix, iy - 1, iz)[0]) / (2.0 * h);
    g[2] = (field.vel_at(ix, iy, iz + 1)[0] - field.vel_at(ix, iy, iz - 1)[0]) / (2.0 * h);
    g[3] = (field.vel_at(ix + 1, iy, iz)[1] - field.vel_at(ix - 1, iy, iz)[1]) / (2.0 * h);
    g[4] = (field.vel_at(ix, iy + 1, iz)[1] - field.vel_at(ix, iy - 1, iz)[1]) / (2.0 * h);
    g[5] = (field.vel_at(ix, iy, iz + 1)[1] - field.vel_at(ix, iy, iz - 1)[1]) / (2.0 * h);
    g[6] = (field.vel_at(ix + 1, iy, iz)[2] - field.vel_at(ix - 1, iy, iz)[2]) / (2.0 * h);
    g[7] = (field.vel_at(ix, iy + 1, iz)[2] - field.vel_at(ix, iy - 1, iz)[2]) / (2.0 * h);
    g[8] = (field.vel_at(ix, iy, iz + 1)[2] - field.vel_at(ix, iy, iz - 1)[2]) / (2.0 * h);
    g
}

// ─────────────────────────────────────────────────────────────────────────────
// Q-criterion
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Q-criterion `Q = 0.5 * (|Omega|^2 - |S|^2)` at each grid point.
///
/// `Omega` is the anti-symmetric part (vorticity tensor) and `S` is the symmetric
/// part (strain-rate tensor) of the velocity gradient.
/// Vortices are identified where `Q > 0`.
pub fn q_criterion(field: &FlowField3D) -> Vec<f64> {
    let nx = field.nx;
    let ny = field.ny;
    let nz = field.nz;
    let mut out = vec![0.0_f64; nx * ny * nz];

    for iz in 0..nz as i64 {
        for iy in 0..ny as i64 {
            for ix in 0..nx as i64 {
                let g = velocity_gradient(field, ix, iy, iz);
                let mut s_frob_sq = 0.0_f64;
                let mut omega_frob_sq = 0.0_f64;
                for i in 0..3 {
                    for j in 0..3 {
                        let gij = g[i * 3 + j];
                        let gji = g[j * 3 + i];
                        let s = (gij + gji) / 2.0;
                        let om = (gij - gji) / 2.0;
                        s_frob_sq += s * s;
                        omega_frob_sq += om * om;
                    }
                }
                let idx = field.idx_clamped(ix, iy, iz);
                out[idx] = 0.5 * (omega_frob_sq - s_frob_sq);
            }
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Lambda-2 vortex criterion
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the lambda-2 vortex criterion at each grid point.
///
/// The lambda-2 criterion identifies vortex cores as regions where the second
/// eigenvalue of `S^2 + Omega^2` is negative, where `S` and `Omega` are the
/// symmetric and anti-symmetric parts of the velocity gradient tensor.
#[allow(clippy::too_many_arguments)]
pub fn lambda2_criterion(field: &FlowField3D) -> Vec<f64> {
    let nx = field.nx;
    let ny = field.ny;
    let nz = field.nz;
    let mut out = vec![0.0_f64; nx * ny * nz];

    for iz in 0..nz as i64 {
        for iy in 0..ny as i64 {
            for ix in 0..nx as i64 {
                let g = velocity_gradient(field, ix, iy, iz);
                let mut m = [0.0_f64; 9];
                for i in 0..3 {
                    for j in 0..3 {
                        let mut val = 0.0;
                        for kk in 0..3 {
                            let s_ik = (g[i * 3 + kk] + g[kk * 3 + i]) / 2.0;
                            let s_kj = (g[kk * 3 + j] + g[j * 3 + kk]) / 2.0;
                            let om_ik = (g[i * 3 + kk] - g[kk * 3 + i]) / 2.0;
                            let om_kj = (g[kk * 3 + j] - g[j * 3 + kk]) / 2.0;
                            val += s_ik * s_kj + om_ik * om_kj;
                        }
                        m[i * 3 + j] = val;
                    }
                }
                let lambda2 = middle_eigenvalue_3x3_sym(&m);
                let idx = field.idx_clamped(ix, iy, iz);
                out[idx] = lambda2;
            }
        }
    }
    out
}

/// Compute the middle (second) eigenvalue of a 3x3 symmetric matrix.
///
/// Uses the analytic method based on the characteristic polynomial via
/// trigonometric solution (Cardano/cosine method).
fn middle_eigenvalue_3x3_sym(m: &[f64; 9]) -> f64 {
    let tr = m[0] + m[4] + m[8];
    let p1 = m[1] * m[1] + m[2] * m[2] + m[5] * m[5];
    if p1.abs() < 1e-20 {
        let mut eigs = [m[0], m[4], m[8]];
        eigs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return eigs[1];
    }
    let q = (tr * tr
        - 3.0
            * (m[0] * m[4] + m[4] * m[8] + m[0] * m[8] - m[1] * m[1] - m[2] * m[2] - m[5] * m[5]))
        / 9.0;
    let r = (2.0 * tr.powi(3)
        - 9.0
            * tr
            * (m[0] * m[4] + m[4] * m[8] + m[0] * m[8] - m[1] * m[1] - m[2] * m[2] - m[5] * m[5])
        + 27.0
            * (m[0] * (m[4] * m[8] - m[5] * m[5]) - m[1] * (m[1] * m[8] - m[5] * m[2])
                + m[2] * (m[1] * m[5] - m[4] * m[2])))
        / 54.0;
    let q3 = q.powi(3);
    if q3 < 0.0 {
        return tr / 3.0;
    }
    let phi = if q3.sqrt() < 1e-30 {
        0.0
    } else {
        (r / q3.sqrt()).clamp(-1.0, 1.0).acos()
    };
    let sq = q.max(0.0).sqrt();
    let eig0 = tr / 3.0 - 2.0 * sq * (phi / 3.0).cos();
    let eig1 = tr / 3.0 - 2.0 * sq * ((phi - 2.0 * std::f64::consts::PI) / 3.0).cos();
    let eig2 = tr / 3.0 - 2.0 * sq * ((phi + 2.0 * std::f64::consts::PI) / 3.0).cos();
    let mut eigs = [eig0, eig1, eig2];
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    eigs[1]
}

// ─────────────────────────────────────────────────────────────────────────────
// Line Integral Convolution (LIC) - 2-D
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a 2-D Line Integral Convolution (LIC) texture.
///
/// `u` and `v` are the 2-D velocity components in row-major order
/// (index = `iy * nx + ix`). A white-noise texture is generated internally.
/// The LIC value at each pixel is the average of `n_steps` samples along the
/// streamline in both forward and backward directions.
///
/// Returns a `Vec`f64` of length `nx * ny` in the range `\[0, 1\]`.
pub fn lic_2d(u: &[f64], v: &[f64], nx: usize, ny: usize, n_steps: usize) -> Vec<f64> {
    assert_eq!(u.len(), nx * ny);
    assert_eq!(v.len(), nx * ny);

    let mut rng = rand::rng();
    let noise: Vec<f64> = (0..nx * ny)
        .map(|_| rng.random_range(0.0..1.0_f64))
        .collect();
    let mut out = vec![0.0_f64; nx * ny];

    let sample_noise = |px: f64, py: f64| -> f64 {
        let ix = (px as i64).clamp(0, nx as i64 - 1) as usize;
        let iy = (py as i64).clamp(0, ny as i64 - 1) as usize;
        noise[iy * nx + ix]
    };

    let sample_vel = |px: f64, py: f64| -> [f64; 2] {
        let x0 = (px.floor() as i64).clamp(0, nx as i64 - 1) as usize;
        let y0 = (py.floor() as i64).clamp(0, ny as i64 - 1) as usize;
        let tx = px - x0 as f64;
        let ty = py - y0 as f64;
        let x1 = (x0 + 1).min(nx - 1);
        let y1 = (y0 + 1).min(ny - 1);
        let u00 = u[y0 * nx + x0];
        let u10 = u[y0 * nx + x1];
        let u01 = u[y1 * nx + x0];
        let u11 = u[y1 * nx + x1];
        let v00 = v[y0 * nx + x0];
        let v10 = v[y0 * nx + x1];
        let v01 = v[y1 * nx + x0];
        let v11 = v[y1 * nx + x1];
        let ui = (u00 * (1.0 - tx) + u10 * tx) * (1.0 - ty) + (u01 * (1.0 - tx) + u11 * tx) * ty;
        let vi = (v00 * (1.0 - tx) + v10 * tx) * (1.0 - ty) + (v01 * (1.0 - tx) + v11 * tx) * ty;
        [ui, vi]
    };

    for iy in 0..ny {
        for ix in 0..nx {
            let mut acc = 0.0_f64;
            let mut count = 0usize;
            // Forward
            let mut px = ix as f64 + 0.5;
            let mut py = iy as f64 + 0.5;
            for _ in 0..n_steps {
                if px < 0.0 || px >= nx as f64 || py < 0.0 || py >= ny as f64 {
                    break;
                }
                acc += sample_noise(px, py);
                count += 1;
                let [ui, vi] = sample_vel(px, py);
                let len = (ui * ui + vi * vi).sqrt().max(1e-12);
                px += ui / len;
                py += vi / len;
            }
            // Backward
            let mut px = ix as f64 + 0.5;
            let mut py = iy as f64 + 0.5;
            for _ in 0..n_steps {
                if px < 0.0 || px >= nx as f64 || py < 0.0 || py >= ny as f64 {
                    break;
                }
                acc += sample_noise(px, py);
                count += 1;
                let [ui, vi] = sample_vel(px, py);
                let len = (ui * ui + vi * vi).sqrt().max(1e-12);
                px -= ui / len;
                py -= vi / len;
            }
            out[iy * nx + ix] = if count > 0 { acc / count as f64 } else { 0.0 };
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Seed point generation
// ─────────────────────────────────────────────────────────────────────────────

/// Generate `n` uniformly distributed seed points within an axis-aligned box.
///
/// `bounds` is `\[\[x_min, x_max\\], \[y_min, y_max\], \[z_min, z_max\]]`.
pub fn seeding_uniform(bounds: [[f64; 2]; 3], n: usize) -> Vec<[f64; 3]> {
    let mut rng = rand::rng();
    (0..n)
        .map(|_| {
            [
                rng.random_range(bounds[0][0]..bounds[0][1]),
                rng.random_range(bounds[1][0]..bounds[1][1]),
                rng.random_range(bounds[2][0]..bounds[2][1]),
            ]
        })
        .collect()
}

/// Generate seed points along a line from `start` to `end`.
///
/// Returns `n` evenly spaced points including both endpoints.
pub fn seeding_line(start: [f64; 3], end: [f64; 3], n: usize) -> Vec<[f64; 3]> {
    if n <= 1 {
        return vec![start];
    }
    let inv = 1.0 / (n - 1) as f64;
    (0..n)
        .map(|i| {
            let t = i as f64 * inv;
            [
                start[0] + t * (end[0] - start[0]),
                start[1] + t * (end[1] - start[1]),
                start[2] + t * (end[2] - start[2]),
            ]
        })
        .collect()
}

/// Generate seed points on a circle in the XY plane.
///
/// `center` is the center position, `radius` the circle radius, and `n`
/// the number of seed points.
pub fn seeding_circle(center: [f64; 3], radius: f64, n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            [
                center[0] + radius * theta.cos(),
                center[1] + radius * theta.sin(),
                center[2],
            ]
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector field arrows
// ─────────────────────────────────────────────────────────────────────────────

/// A single vector-field arrow defined by a tail position and a head position.
#[derive(Debug, Clone, Copy)]
pub struct VectorArrow {
    /// Arrow tail position (world coordinates).
    pub tail: [f64; 3],
    /// Arrow head position (world coordinates).
    pub head: [f64; 3],
    /// Velocity magnitude at the tail point.
    pub magnitude: f64,
}

/// Generate vector field arrows at regularly-spaced grid points.
///
/// `stride` controls subsampling: arrows are placed every `stride` grid points
/// along each axis. The `scale` factor controls arrow length.
pub fn vector_field_arrows(field: &FlowField3D, stride: usize, scale: f64) -> Vec<VectorArrow> {
    let stride = stride.max(1);
    let mut arrows = Vec::new();

    let mut iz = 0;
    while iz < field.nz {
        let mut iy = 0;
        while iy < field.ny {
            let mut ix = 0;
            while ix < field.nx {
                let vel = field.vel_at(ix as i64, iy as i64, iz as i64);
                let m = mag3(vel);
                if m > 1e-12 {
                    let tail = [
                        ix as f64 * field.dx,
                        iy as f64 * field.dx,
                        iz as f64 * field.dx,
                    ];
                    let head = [
                        tail[0] + vel[0] * scale,
                        tail[1] + vel[1] * scale,
                        tail[2] + vel[2] * scale,
                    ];
                    arrows.push(VectorArrow {
                        tail,
                        head,
                        magnitude: m,
                    });
                }
                ix += stride;
            }
            iy += stride;
        }
        iz += stride;
    }
    arrows
}

/// Generate a single arrow mesh as line segments (tail, head, and two barbs).
///
/// `barb_fraction` controls the relative size of the arrowhead barbs
/// (typically 0.2 to 0.3). Returns up to 3 line-segment pairs.
pub fn arrow_segments(
    tail: [f64; 3],
    head: [f64; 3],
    barb_fraction: f64,
) -> Vec<([f64; 3], [f64; 3])> {
    let dir = sub3(head, tail);
    let len = mag3(dir);
    if len < 1e-15 {
        return Vec::new();
    }
    let d = normalize3(dir);

    // Find a perpendicular direction
    let up = if d[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let perp = normalize3(cross3(d, up));

    let barb_len = len * barb_fraction;
    let barb_base = sub3(head, scale3(d, barb_len));
    let barb1 = add3(barb_base, scale3(perp, barb_len * 0.5));
    let barb2 = sub3(barb_base, scale3(perp, barb_len * 0.5));

    vec![(tail, head), (head, barb1), (head, barb2)]
}

// ─────────────────────────────────────────────────────────────────────────────
// Flow topology: critical points and separatrices
// ─────────────────────────────────────────────────────────────────────────────

/// Classification of a critical point in a 2-D or 3-D flow field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalPointType {
    /// An attracting node (sink).
    AttractingNode,
    /// A repelling node (source).
    RepellingNode,
    /// A saddle point.
    Saddle,
    /// An attracting focus (spiral sink).
    AttractingFocus,
    /// A repelling focus (spiral source).
    RepellingFocus,
    /// A center (closed orbits).
    Center,
    /// Degenerate or indeterminate type.
    Degenerate,
}

/// A critical point in the flow field with its classification.
#[derive(Debug, Clone)]
pub struct CriticalPoint {
    /// Position in world coordinates.
    pub position: [f64; 3],
    /// Grid coordinates (fractional).
    pub grid_coords: [f64; 3],
    /// Classification type.
    pub point_type: CriticalPointType,
    /// The velocity gradient tensor at this point (3x3, row-major).
    pub gradient: [f64; 9],
}

/// Find grid points where the velocity magnitude is approximately zero.
///
/// A point is classified as a critical point if `|v| < threshold` where
/// `threshold` defaults to `1e-4` (grid-unit velocity). Returns the
/// world-space positions of all such grid points.
pub fn flow_topology_critical_points(field: &FlowField3D) -> Vec<[f64; 3]> {
    let threshold = 1e-4_f64;
    let mut cps = Vec::new();
    for iz in 0..field.nz {
        for iy in 0..field.ny {
            for ix in 0..field.nx {
                let idx = field.idx_clamped(ix as i64, iy as i64, iz as i64);
                let u = field.u[idx];
                let v = field.v[idx];
                let w = field.w[idx];
                if (u * u + v * v + w * w).sqrt() < threshold {
                    cps.push([
                        ix as f64 * field.dx,
                        iy as f64 * field.dx,
                        iz as f64 * field.dx,
                    ]);
                }
            }
        }
    }
    cps
}

/// Find and classify critical points in the flow field.
///
/// Uses velocity magnitude threshold and the velocity gradient eigenvalue
/// structure to classify each critical point.
pub fn classify_critical_points(field: &FlowField3D) -> Vec<CriticalPoint> {
    let threshold = 1e-4_f64;
    let mut cps = Vec::new();
    for iz in 0..field.nz {
        for iy in 0..field.ny {
            for ix in 0..field.nx {
                let idx = field.idx_clamped(ix as i64, iy as i64, iz as i64);
                let u = field.u[idx];
                let v = field.v[idx];
                let w = field.w[idx];
                if (u * u + v * v + w * w).sqrt() < threshold {
                    let grad = velocity_gradient(field, ix as i64, iy as i64, iz as i64);
                    let pt = classify_from_gradient(&grad);
                    cps.push(CriticalPoint {
                        position: [
                            ix as f64 * field.dx,
                            iy as f64 * field.dx,
                            iz as f64 * field.dx,
                        ],
                        grid_coords: [ix as f64, iy as f64, iz as f64],
                        point_type: pt,
                        gradient: grad,
                    });
                }
            }
        }
    }
    cps
}

/// Classify a critical point from its 3x3 velocity gradient tensor.
///
/// Uses the trace (sum of real eigenvalue parts) and determinant as proxies.
fn classify_from_gradient(g: &[f64; 9]) -> CriticalPointType {
    let tr = g[0] + g[4] + g[8];
    let det = g[0] * (g[4] * g[8] - g[5] * g[7]) - g[1] * (g[3] * g[8] - g[5] * g[6])
        + g[2] * (g[3] * g[7] - g[4] * g[6]);

    // Discriminant of characteristic polynomial for 2D-like classification
    let p = -(g[0] * g[4] + g[4] * g[8] + g[0] * g[8] - g[1] * g[3] - g[2] * g[6] - g[5] * g[7]);
    let disc = tr * tr - 4.0 * p;

    if det.abs() < 1e-12 {
        return CriticalPointType::Degenerate;
    }

    if disc >= 0.0 {
        // Real eigenvalues
        if det > 0.0 && tr < -1e-10 {
            CriticalPointType::AttractingNode
        } else if det > 0.0 && tr > 1e-10 {
            CriticalPointType::RepellingNode
        } else if det < 0.0 {
            CriticalPointType::Saddle
        } else {
            CriticalPointType::Degenerate
        }
    } else {
        // Complex eigenvalues
        if tr < -1e-10 {
            CriticalPointType::AttractingFocus
        } else if tr > 1e-10 {
            CriticalPointType::RepellingFocus
        } else {
            CriticalPointType::Center
        }
    }
}

/// Trace separatrices from saddle points in the field.
///
/// For each saddle point, traces streamlines in both the stable and unstable
/// eigenvector directions using RK4. Returns pairs of (forward, backward) lines.
pub fn trace_separatrices(
    field: &FlowField3D,
    saddle_points: &[CriticalPoint],
    dt: f64,
    n_steps: usize,
) -> Vec<(Vec<[f64; 3]>, Vec<[f64; 3]>)> {
    let eps = field.dx * 0.01; // small offset from critical point
    let mut result = Vec::new();

    for cp in saddle_points {
        if cp.point_type != CriticalPointType::Saddle {
            continue;
        }
        // Use the gradient to determine principal directions
        let g = &cp.gradient;
        // Approximate eigenvector from first row
        let ev = normalize3([g[0], g[1], g[2]]);
        let seed_fwd = add3(cp.grid_coords, scale3(ev, eps / field.dx));
        let seed_bwd = sub3(cp.grid_coords, scale3(ev, eps / field.dx));

        let fwd = streamline_rk4(field, seed_fwd, dt, n_steps);
        let bwd = streamline_rk4(field, seed_bwd, -dt, n_steps);
        result.push((fwd, bwd));
    }
    result
}

/// Compute the strain-rate magnitude (Frobenius norm of symmetric part) field.
///
/// Returns a flat scalar array of the same length as `field.u`.
pub fn strain_rate_magnitude(field: &FlowField3D) -> Vec<f64> {
    let nx = field.nx;
    let ny = field.ny;
    let nz = field.nz;
    let mut out = vec![0.0_f64; nx * ny * nz];

    for iz in 0..nz as i64 {
        for iy in 0..ny as i64 {
            for ix in 0..nx as i64 {
                let g = velocity_gradient(field, ix, iy, iz);
                let mut s_frob_sq = 0.0_f64;
                for i in 0..3 {
                    for j in 0..3 {
                        let s = (g[i * 3 + j] + g[j * 3 + i]) / 2.0;
                        s_frob_sq += s * s;
                    }
                }
                let idx = field.idx_clamped(ix, iy, iz);
                out[idx] = s_frob_sq.sqrt();
            }
        }
    }
    out
}

/// Compute the enstrophy field (0.5 * |omega|^2) at each grid point.
///
/// Enstrophy is a measure of the intensity of vorticity in the flow.
pub fn enstrophy(field: &FlowField3D) -> Vec<f64> {
    let vor = vorticity(field);
    let n = vor.len();
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        out[i] = 0.5 * (vor.u[i] * vor.u[i] + vor.v[i] * vor.v[i] + vor.w[i] * vor.w[i]);
    }
    out
}

/// Compute the helicity field `H = v . omega` at each grid point.
///
/// Helicity measures the degree of linkage between vortex lines.
pub fn helicity(field: &FlowField3D) -> Vec<f64> {
    let vor = vorticity(field);
    let n = field.len();
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        out[i] = field.u[i] * vor.u[i] + field.v[i] * vor.v[i] + field.w[i] * vor.w[i];
    }
    out
}

/// Compute the speed (velocity magnitude) field.
pub fn speed_field(field: &FlowField3D) -> Vec<f64> {
    let n = field.len();
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        out[i] =
            (field.u[i] * field.u[i] + field.v[i] * field.v[i] + field.w[i] * field.w[i]).sqrt();
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a uniform flow field with constant velocity `(u0, v0, w0)`.
    #[allow(clippy::too_many_arguments)]
    fn uniform_field(
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        u0: f64,
        v0: f64,
        w0: f64,
    ) -> FlowField3D {
        let n = nx * ny * nz;
        FlowField3D {
            u: vec![u0; n],
            v: vec![v0; n],
            w: vec![w0; n],
            nx,
            ny,
            nz,
            dx,
        }
    }

    // ── FlowField3D basics ────────────────────────────────────────────────────

    #[test]
    fn test_flow_field_new_zero() {
        let f = FlowField3D::new(4, 4, 4, 1.0);
        assert_eq!(f.len(), 64);
        assert!(f.u.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_flow_field_empty() {
        let f = FlowField3D::new(0, 0, 0, 1.0);
        assert!(f.is_empty());
    }

    #[test]
    fn test_flow_field_vel_at_clamped() {
        let f = uniform_field(4, 4, 4, 1.0, 1.0, 2.0, 3.0);
        let v = f.vel_at(-1, 0, 0);
        assert_eq!(v, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_flow_field_set_vel() {
        let mut f = FlowField3D::new(4, 4, 4, 1.0);
        f.set_vel(1, 2, 3, [5.0, 6.0, 7.0]);
        let v = f.vel_at(1, 2, 3);
        assert_eq!(v, [5.0, 6.0, 7.0]);
    }

    #[test]
    fn test_flow_field_speed_at() {
        let f = uniform_field(4, 4, 4, 1.0, 3.0, 4.0, 0.0);
        let s = f.speed_at(0, 0, 0);
        assert!((s - 5.0).abs() < 1e-10);
    }

    // ── interpolate_velocity ──────────────────────────────────────────────────

    #[test]
    fn test_interpolate_uniform_field() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let v = interpolate_velocity(&f, [2.3, 2.3, 2.3]);
        assert!((v[0] - 1.0).abs() < 1e-10, "uniform field: {}", v[0]);
    }

    #[test]
    fn test_interpolate_zero_field() {
        let f = FlowField3D::new(5, 5, 5, 1.0);
        let v = interpolate_velocity(&f, [2.5, 2.5, 2.5]);
        assert_eq!(v, [0.0, 0.0, 0.0]);
    }

    // ── streamline_rk4 ────────────────────────────────────────────────────────

    #[test]
    fn test_streamline_rk4_uniform_field() {
        let f = uniform_field(10, 10, 10, 1.0, 1.0, 0.0, 0.0);
        let path = streamline_rk4(&f, [0.5, 5.0, 5.0], 0.5, 4);
        assert_eq!(path.len(), 5, "should have seed + 4 steps");
        assert!(path[4][0] > path[0][0], "should advance in x");
    }

    #[test]
    fn test_streamline_rk4_zero_field_stays_put() {
        let f = FlowField3D::new(10, 10, 10, 1.0);
        let start = [5.0, 5.0, 5.0];
        let path = streamline_rk4(&f, start, 0.1, 5);
        for p in &path {
            assert!((p[0] - start[0]).abs() < 1e-12);
            assert!((p[1] - start[1]).abs() < 1e-12);
            assert!((p[2] - start[2]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_streamline_rk4_step_count() {
        let f = uniform_field(10, 10, 10, 1.0, 0.0, 1.0, 0.0);
        let path = streamline_rk4(&f, [5.0, 0.5, 5.0], 0.1, 10);
        assert_eq!(path.len(), 11);
    }

    #[test]
    fn test_streamline_rk4_uniform_displacement() {
        let f = uniform_field(20, 20, 20, 1.0, 1.0, 0.0, 0.0);
        let path = streamline_rk4(&f, [5.0, 10.0, 10.0], 1.0, 3);
        assert!((path[1][0] - path[0][0] - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_streamline_rk4_backward() {
        let f = uniform_field(10, 10, 10, 1.0, 1.0, 0.0, 0.0);
        let path = streamline_rk4_backward(&f, [5.0, 5.0, 5.0], 0.5, 4);
        assert_eq!(path.len(), 5);
        // Should move backward in x
        assert!(path[4][0] < path[0][0], "should go backward in x");
    }

    // ── pathline_euler ────────────────────────────────────────────────────────

    #[test]
    fn test_pathline_euler_uniform_field() {
        let f = uniform_field(10, 10, 10, 1.0, 0.0, 1.0, 0.0);
        let path = pathline_euler(&f, [5.0, 0.5, 5.0], 0.5, 4);
        assert_eq!(path.len(), 5);
        assert!(path[4][1] > path[0][1], "should advance in y");
    }

    #[test]
    fn test_pathline_euler_zero_stays_put() {
        let f = FlowField3D::new(8, 8, 8, 1.0);
        let start = [4.0, 4.0, 4.0];
        let path = pathline_euler(&f, start, 0.1, 5);
        for p in &path {
            assert!((p[0] - start[0]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_pathline_euler_step_count() {
        let f = uniform_field(8, 8, 8, 1.0, 1.0, 0.0, 0.0);
        let path = pathline_euler(&f, [1.0, 4.0, 4.0], 0.1, 6);
        assert_eq!(path.len(), 7);
    }

    // ── streakline ───────────────────────────────────────────────────────────

    #[test]
    fn test_streakline_advance_adds_particles() {
        let f = uniform_field(10, 10, 10, 1.0, 1.0, 0.0, 0.0);
        let mut sl = Streakline::new([5.0, 5.0, 5.0]);
        assert!(sl.is_empty());
        sl.advance(&f, 0.1);
        assert_eq!(sl.len(), 1);
        sl.advance(&f, 0.1);
        assert_eq!(sl.len(), 2);
    }

    #[test]
    fn test_streakline_particles_move() {
        let f = uniform_field(10, 10, 10, 1.0, 1.0, 0.0, 0.0);
        let mut sl = Streakline::new([5.0, 5.0, 5.0]);
        sl.advance(&f, 0.1);
        sl.advance(&f, 0.1);
        // First particle was injected at step 0, advected at step 1
        // It should have moved in x
        let first = sl.particles[0];
        assert!(first[0] > 5.0, "first particle should have moved in x");
    }

    // ── vorticity ─────────────────────────────────────────────────────────────

    #[test]
    fn test_vorticity_uniform_field_is_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 2.0, 3.0);
        let vor = vorticity(&f);
        for i in 0..vor.len() {
            assert!(
                vor.u[i].abs() < 1e-10,
                "omega_x should be 0, got {}",
                vor.u[i]
            );
            assert!(
                vor.v[i].abs() < 1e-10,
                "omega_y should be 0, got {}",
                vor.v[i]
            );
            assert!(
                vor.w[i].abs() < 1e-10,
                "omega_z should be 0, got {}",
                vor.w[i]
            );
        }
    }

    #[test]
    fn test_vorticity_same_dimensions() {
        let f = uniform_field(4, 5, 6, 1.0, 1.0, 0.0, 0.0);
        let vor = vorticity(&f);
        assert_eq!(vor.nx, 4);
        assert_eq!(vor.ny, 5);
        assert_eq!(vor.nz, 6);
    }

    #[test]
    fn test_vorticity_magnitude_uniform_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let mag = vorticity_magnitude(&f);
        for &m in &mag {
            assert!(
                m.abs() < 1e-10,
                "vorticity magnitude should be 0 in uniform, got {m}"
            );
        }
    }

    // ── divergence ────────────────────────────────────────────────────────────

    #[test]
    fn test_divergence_uniform_field_is_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 2.0, 3.0);
        let div = divergence(&f);
        for &d in &div {
            assert!(
                d.abs() < 1e-10,
                "uniform field should have zero divergence, got {d}"
            );
        }
    }

    #[test]
    fn test_divergence_length() {
        let f = FlowField3D::new(4, 4, 4, 1.0);
        let div = divergence(&f);
        assert_eq!(div.len(), 64);
    }

    #[test]
    fn test_divergence_linear_x_field() {
        let nx = 6;
        let ny = 6;
        let nz = 6;
        let mut f = FlowField3D::new(nx, ny, nz, 1.0);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = iz * ny * nx + iy * nx + ix;
                    f.u[idx] = ix as f64;
                }
            }
        }
        let div = divergence(&f);
        let idx = f.idx_clamped(3, 3, 3);
        assert!(
            (div[idx] - 1.0).abs() < 1e-10,
            "linear x-field div = {}",
            div[idx]
        );
    }

    // ── q_criterion ───────────────────────────────────────────────────────────

    #[test]
    fn test_q_criterion_uniform_is_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let q = q_criterion(&f);
        for &qi in &q {
            assert!(qi.abs() < 1e-10, "uniform flow has Q=0, got {qi}");
        }
    }

    #[test]
    fn test_q_criterion_length() {
        let f = FlowField3D::new(3, 3, 3, 1.0);
        let q = q_criterion(&f);
        assert_eq!(q.len(), 27);
    }

    // ── lambda2_criterion ─────────────────────────────────────────────────────

    #[test]
    fn test_lambda2_criterion_length() {
        let f = FlowField3D::new(4, 4, 4, 1.0);
        let l2 = lambda2_criterion(&f);
        assert_eq!(l2.len(), 64);
    }

    #[test]
    fn test_lambda2_uniform_field_no_vortex() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let l2 = lambda2_criterion(&f);
        for &v in &l2 {
            assert!(v.abs() < 1e-8, "uniform field: lambda2 = {v}");
        }
    }

    // ── lic_2d ────────────────────────────────────────────────────────────────

    #[test]
    fn test_lic_2d_length() {
        let nx = 8;
        let ny = 8;
        let u = vec![1.0_f64; nx * ny];
        let v_comp = vec![0.0_f64; nx * ny];
        let lic = lic_2d(&u, &v_comp, nx, ny, 3);
        assert_eq!(lic.len(), nx * ny);
    }

    #[test]
    fn test_lic_2d_values_in_range() {
        let nx = 4;
        let ny = 4;
        let u = vec![1.0_f64; nx * ny];
        let v_comp = vec![0.0_f64; nx * ny];
        let lic = lic_2d(&u, &v_comp, nx, ny, 5);
        for &val in &lic {
            assert!((0.0..=1.0).contains(&val), "LIC value out of [0,1]: {val}");
        }
    }

    // ── seeding ─────────────────────────────────────────────────────────────

    #[test]
    fn test_seeding_uniform_count() {
        let seeds = seeding_uniform([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], 50);
        assert_eq!(seeds.len(), 50);
    }

    #[test]
    fn test_seeding_uniform_in_bounds() {
        let bounds = [[1.0, 3.0], [2.0, 5.0], [0.0, 1.0]];
        let seeds = seeding_uniform(bounds, 100);
        for p in &seeds {
            assert!(p[0] >= 1.0 && p[0] < 3.0);
            assert!(p[1] >= 2.0 && p[1] < 5.0);
            assert!(p[2] >= 0.0 && p[2] < 1.0);
        }
    }

    #[test]
    fn test_seeding_uniform_zero_seeds() {
        let seeds = seeding_uniform([[0.0, 1.0]; 3], 0);
        assert!(seeds.is_empty());
    }

    #[test]
    fn test_seeding_line_endpoints() {
        let seeds = seeding_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 5);
        assert_eq!(seeds.len(), 5);
        assert!((seeds[0][0] - 0.0).abs() < 1e-10);
        assert!((seeds[4][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_seeding_circle_count() {
        let seeds = seeding_circle([0.0, 0.0, 0.0], 1.0, 8);
        assert_eq!(seeds.len(), 8);
        // All should be at radius 1 in XY plane
        for p in &seeds {
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!((r - 1.0).abs() < 1e-10, "point not on circle: r={r}");
        }
    }

    // ── vector field arrows ──────────────────────────────────────────────────

    #[test]
    fn test_vector_field_arrows_uniform() {
        let f = uniform_field(4, 4, 4, 1.0, 1.0, 0.0, 0.0);
        let arrows = vector_field_arrows(&f, 1, 1.0);
        assert!(!arrows.is_empty());
        for a in &arrows {
            assert!((a.magnitude - 1.0).abs() < 1e-10);
            assert!(a.head[0] > a.tail[0], "arrow should point in +x");
        }
    }

    #[test]
    fn test_vector_field_arrows_zero_field() {
        let f = FlowField3D::new(4, 4, 4, 1.0);
        let arrows = vector_field_arrows(&f, 1, 1.0);
        assert!(arrows.is_empty(), "zero field should produce no arrows");
    }

    #[test]
    fn test_arrow_segments_basic() {
        let segs = arrow_segments([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.3);
        assert_eq!(segs.len(), 3, "shaft + 2 barbs");
    }

    #[test]
    fn test_arrow_segments_zero_length() {
        let segs = arrow_segments([1.0, 1.0, 1.0], [1.0, 1.0, 1.0], 0.3);
        assert!(segs.is_empty());
    }

    // ── flow topology ────────────────────────────────────────────────────────

    #[test]
    fn test_critical_points_zero_field() {
        let f = FlowField3D::new(3, 3, 3, 1.0);
        let cps = flow_topology_critical_points(&f);
        assert_eq!(cps.len(), 27, "all zero-field grid points are critical");
    }

    #[test]
    fn test_critical_points_uniform_field_none() {
        let f = uniform_field(4, 4, 4, 1.0, 5.0, 0.0, 0.0);
        let cps = flow_topology_critical_points(&f);
        assert!(
            cps.is_empty(),
            "uniform non-zero field has no critical points"
        );
    }

    #[test]
    fn test_critical_points_partial() {
        let mut f = FlowField3D::new(3, 3, 3, 1.0);
        f.u[13] = 1.0;
        let cps = flow_topology_critical_points(&f);
        assert_eq!(cps.len(), 26);
    }

    #[test]
    fn test_critical_points_positions_in_bounds() {
        let f = FlowField3D::new(4, 4, 4, 0.5);
        let cps = flow_topology_critical_points(&f);
        let extent = 4.0 * 0.5;
        for p in &cps {
            assert!(p[0] < extent + 1e-6);
            assert!(p[1] < extent + 1e-6);
            assert!(p[2] < extent + 1e-6);
        }
    }

    #[test]
    fn test_classify_critical_points_zero_field() {
        let f = FlowField3D::new(3, 3, 3, 1.0);
        let cps = classify_critical_points(&f);
        assert_eq!(cps.len(), 27);
        // All zero gradient -> degenerate
        for cp in &cps {
            assert_eq!(cp.point_type, CriticalPointType::Degenerate);
        }
    }

    #[test]
    fn test_separatrices_no_saddles() {
        let f = FlowField3D::new(5, 5, 5, 1.0);
        let cps = classify_critical_points(&f);
        let seps = trace_separatrices(&f, &cps, 0.1, 10);
        assert!(seps.is_empty(), "no saddles -> no separatrices");
    }

    // ── derived fields ───────────────────────────────────────────────────────

    #[test]
    fn test_strain_rate_magnitude_uniform_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let sr = strain_rate_magnitude(&f);
        for &val in &sr {
            assert!(
                val.abs() < 1e-10,
                "uniform field strain rate should be 0, got {val}"
            );
        }
    }

    #[test]
    fn test_enstrophy_uniform_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let ens = enstrophy(&f);
        for &val in &ens {
            assert!(
                val.abs() < 1e-10,
                "uniform field enstrophy should be 0, got {val}"
            );
        }
    }

    #[test]
    fn test_helicity_uniform_zero() {
        let f = uniform_field(5, 5, 5, 1.0, 1.0, 0.0, 0.0);
        let h = helicity(&f);
        for &val in &h {
            assert!(
                val.abs() < 1e-10,
                "uniform field helicity should be 0, got {val}"
            );
        }
    }

    #[test]
    fn test_speed_field_uniform() {
        let f = uniform_field(4, 4, 4, 1.0, 3.0, 4.0, 0.0);
        let spd = speed_field(&f);
        for &val in &spd {
            assert!((val - 5.0).abs() < 1e-10, "speed should be 5.0, got {val}");
        }
    }

    #[test]
    fn test_speed_field_length() {
        let f = FlowField3D::new(3, 3, 3, 1.0);
        let spd = speed_field(&f);
        assert_eq!(spd.len(), 27);
    }

    // ── vector helpers ───────────────────────────────────────────────────────

    #[test]
    fn test_vec_helpers() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert_eq!(add3(a, b), [5.0, 7.0, 9.0]);
        assert_eq!(sub3(a, b), [-3.0, -3.0, -3.0]);
        assert_eq!(scale3(a, 2.0), [2.0, 4.0, 6.0]);
        assert!((dot3(a, b) - 32.0).abs() < 1e-10);
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3_unit() {
        let n = normalize3([3.0, 4.0, 0.0]);
        assert!((mag3(n) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3_zero() {
        let n = normalize3([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }
}
