// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-side parallel flux and advection compute kernels.
//!
//! Provides scalar and vector flux computation on structured 3-D Cartesian
//! grids: upwind advection, central-difference divergence, Lax-Friedrichs
//! flux splitting, Godunov interface reconstruction, limiter functions
//! (minmod, superbee, van Leer), and a simple finite-volume time-stepping
//! helper.  All kernels use Rayon for data-parallel execution.

use rayon::prelude::*;

// ---------------------------------------------------------------------------
// 1. FluxGrid3D — a 3-D scalar field on a Cartesian grid
// ---------------------------------------------------------------------------

/// A 3-D scalar field stored in row-major order on a uniform Cartesian grid.
///
/// Index layout: `values[i*(ny*nz) + j*nz + k]` for cell `(i, j, k)`.
#[derive(Debug, Clone)]
pub struct FluxGrid3D {
    /// Cells in the x-direction.
    pub nx: usize,
    /// Cells in the y-direction.
    pub ny: usize,
    /// Cells in the z-direction.
    pub nz: usize,
    /// Uniform cell spacing.
    pub dx: f64,
    /// Scalar field values.
    pub values: Vec<f64>,
}

impl FluxGrid3D {
    /// Create a grid filled with `fill`.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, fill: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            values: vec![fill; nx * ny * nz],
        }
    }

    /// Flat index for cell `(i, j, k)`.
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i * self.ny * self.nz + j * self.nz + k
    }

    /// Read cell `(i, j, k)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize, k: usize) -> f64 {
        self.values[self.idx(i, j, k)]
    }

    /// Write cell `(i, j, k)`.
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, k: usize, v: f64) {
        let idx = self.idx(i, j, k);
        self.values[idx] = v;
    }

    /// Total number of cells.
    pub fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// True when the grid has no cells.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Clamp-to-boundary accessor: returns the value at `(i, j, k)` with
    /// zero-gradient boundary conditions.
    pub fn get_clamped(&self, i: isize, j: isize, k: isize) -> f64 {
        let ci = i.clamp(0, self.nx as isize - 1) as usize;
        let cj = j.clamp(0, self.ny as isize - 1) as usize;
        let ck = k.clamp(0, self.nz as isize - 1) as usize;
        self.get(ci, cj, ck)
    }

    /// L2 norm of all values.
    pub fn l2_norm(&self) -> f64 {
        let sum_sq: f64 = self.values.iter().map(|&v| v * v).sum();
        sum_sq.sqrt()
    }

    /// Element-wise maximum value.
    pub fn max_val(&self) -> f64 {
        self.values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Element-wise minimum value.
    pub fn min_val(&self) -> f64 {
        self.values.iter().copied().fold(f64::INFINITY, f64::min)
    }
}

// ---------------------------------------------------------------------------
// 2. VectorFluxGrid3D — three-component (u, v, w) velocity field
// ---------------------------------------------------------------------------

/// Three-component 3-D velocity / flux field.
///
/// Components `u` (x), `v` (y), `w` (z) are stored as separate [`FluxGrid3D`]
/// instances sharing the same dimensions and spacing.
#[derive(Debug, Clone)]
pub struct VectorFluxGrid3D {
    /// x-component.
    pub u: FluxGrid3D,
    /// y-component.
    pub v: FluxGrid3D,
    /// z-component.
    pub w: FluxGrid3D,
}

impl VectorFluxGrid3D {
    /// Create a zero-filled vector field.
    pub fn zeros(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self {
            u: FluxGrid3D::new(nx, ny, nz, dx, 0.0),
            v: FluxGrid3D::new(nx, ny, nz, dx, 0.0),
            w: FluxGrid3D::new(nx, ny, nz, dx, 0.0),
        }
    }

    /// Set velocity at cell `(i, j, k)`.
    pub fn set_vel(&mut self, i: usize, j: usize, k: usize, vel: [f64; 3]) {
        self.u.set(i, j, k, vel[0]);
        self.v.set(i, j, k, vel[1]);
        self.w.set(i, j, k, vel[2]);
    }

    /// Get velocity at cell `(i, j, k)`.
    pub fn get_vel(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.u.get(i, j, k),
            self.v.get(i, j, k),
            self.w.get(i, j, k),
        ]
    }
}

// ---------------------------------------------------------------------------
// 3. Slope limiter functions
// ---------------------------------------------------------------------------

/// Minmod slope limiter: smallest slope among `a` and `b` with same sign.
///
/// Returns zero when `a` and `b` have opposite signs.
#[inline]
pub fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

/// Superbee limiter: most compressive TVD limiter.
pub fn superbee(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        return 0.0;
    }
    let s = if a >= 0.0 { 1.0 } else { -1.0 };
    s * a.abs().max(b.abs()).min(2.0 * a.abs().min(b.abs()))
}

/// Van Leer limiter: smooth TVD limiter.
pub fn van_leer(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        return 0.0;
    }
    2.0 * a * b / (a + b)
}

/// MC (Monotonised Central) limiter.
pub fn mc_limiter(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        return 0.0;
    }
    let centered = 0.5 * (a + b);
    let s = if a >= 0.0 { 1.0 } else { -1.0 };
    s * centered.abs().min(2.0 * a.abs()).min(2.0 * b.abs())
}

// ---------------------------------------------------------------------------
// 4. Upwind advection (1-D per axis, assembled into 3-D)
// ---------------------------------------------------------------------------

/// First-order upwind advection flux for a scalar `phi` transported by
/// uniform velocity `vel` in one dimension.
///
/// Returns `F_right - F_left` (net flux into cell `i`) for each cell.
///
/// Uses zero-gradient boundary conditions.
pub fn upwind_flux_1d(phi: &[f64], vel: f64, dx: f64, dt: f64) -> Vec<f64> {
    let n = phi.len();
    let mut flux = vec![0.0; n];
    let cfl = vel * dt / dx;
    for i in 0..n {
        let phi_l = if i == 0 { phi[0] } else { phi[i - 1] };
        let phi_r = if i + 1 >= n { phi[n - 1] } else { phi[i + 1] };
        if vel >= 0.0 {
            // Upwind from left
            flux[i] = -cfl * (phi[i] - phi_l);
        } else {
            // Upwind from right
            flux[i] = -cfl * (phi_r - phi[i]);
        }
    }
    flux
}

/// Apply first-order upwind advection to a [`FluxGrid3D`] for one time step.
///
/// Returns a new grid with updated values.
pub fn advect_upwind_3d(phi: &FluxGrid3D, vel: &VectorFluxGrid3D, dt: f64) -> FluxGrid3D {
    let _nx = phi.nx;
    let ny = phi.ny;
    let nz = phi.nz;
    let dx = phi.dx;

    let mut out = phi.clone();

    out.values.par_iter_mut().enumerate().for_each(|(idx, v)| {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;

        let ux = vel.u.get(i, j, k);
        let uy = vel.v.get(i, j, k);
        let uz = vel.w.get(i, j, k);

        let phi_xm = phi.get_clamped(i as isize - 1, j as isize, k as isize);
        let phi_xp = phi.get_clamped(i as isize + 1, j as isize, k as isize);
        let phi_ym = phi.get_clamped(i as isize, j as isize - 1, k as isize);
        let phi_yp = phi.get_clamped(i as isize, j as isize + 1, k as isize);
        let phi_zm = phi.get_clamped(i as isize, j as isize, k as isize - 1);
        let phi_zp = phi.get_clamped(i as isize, j as isize, k as isize + 1);
        let phi_c = phi.get(i, j, k);

        let flux_x = if ux >= 0.0 {
            ux * (phi_c - phi_xm) / dx
        } else {
            ux * (phi_xp - phi_c) / dx
        };
        let flux_y = if uy >= 0.0 {
            uy * (phi_c - phi_ym) / dx
        } else {
            uy * (phi_yp - phi_c) / dx
        };
        let flux_z = if uz >= 0.0 {
            uz * (phi_c - phi_zm) / dx
        } else {
            uz * (phi_zp - phi_c) / dx
        };

        *v = phi_c - dt * (flux_x + flux_y + flux_z);
    });

    out
}

// ---------------------------------------------------------------------------
// 5. Lax-Friedrichs flux splitting
// ---------------------------------------------------------------------------

/// Lax-Friedrichs numerical flux at the interface between cells `i` and `i+1`.
///
/// `f(u)` is the physical flux function, `alpha` is the maximum wave speed.
#[inline]
pub fn lax_friedrichs_flux(u_l: f64, u_r: f64, f_l: f64, f_r: f64, alpha: f64) -> f64 {
    0.5 * (f_l + f_r) - 0.5 * alpha * (u_r - u_l)
}

/// Lax-Friedrichs advection update for a 1-D scalar array.
///
/// Physical flux: `f(u) = vel * u`.
/// Returns the updated scalar array after one time step.
pub fn lax_friedrichs_advect_1d(u: &[f64], vel: f64, dx: f64, dt: f64) -> Vec<f64> {
    let n = u.len();
    if n == 0 {
        return Vec::new();
    }
    let alpha = vel.abs();
    let mut u_new = vec![0.0; n];
    for i in 0..n {
        let u_l = if i == 0 { u[0] } else { u[i - 1] };
        let u_r = if i + 1 >= n { u[n - 1] } else { u[i + 1] };
        let f_l = vel * u_l;
        let f_r = vel * u_r;
        let f_left = lax_friedrichs_flux(u_l, u[i], f_l, vel * u[i], alpha);
        let f_right = lax_friedrichs_flux(u[i], u_r, vel * u[i], f_r, alpha);
        u_new[i] = u[i] - dt / dx * (f_right - f_left);
    }
    u_new
}

// ---------------------------------------------------------------------------
// 6. Divergence and gradient on the 3-D grid
// ---------------------------------------------------------------------------

/// Compute the divergence of a vector field at every interior cell.
///
/// Uses central differences: `div(F) ≈ (Fx_{i+1} - Fx_{i-1}) / (2dx) + ...`
pub fn divergence_3d(vel: &VectorFluxGrid3D) -> FluxGrid3D {
    let nx = vel.u.nx;
    let ny = vel.u.ny;
    let nz = vel.u.nz;
    let dx = vel.u.dx;

    let mut div = FluxGrid3D::new(nx, ny, nz, dx, 0.0);
    div.values.par_iter_mut().enumerate().for_each(|(idx, d)| {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let ii = i as isize;
        let jj = j as isize;
        let kk = k as isize;
        let du_dx =
            (vel.u.get_clamped(ii + 1, jj, kk) - vel.u.get_clamped(ii - 1, jj, kk)) / (2.0 * dx);
        let dv_dy =
            (vel.v.get_clamped(ii, jj + 1, kk) - vel.v.get_clamped(ii, jj - 1, kk)) / (2.0 * dx);
        let dw_dz =
            (vel.w.get_clamped(ii, jj, kk + 1) - vel.w.get_clamped(ii, jj, kk - 1)) / (2.0 * dx);
        *d = du_dx + dv_dy + dw_dz;
    });
    div
}

/// Compute the gradient of a scalar field at every cell.
///
/// Returns a [`VectorFluxGrid3D`] where each component is a partial derivative.
pub fn gradient_3d(phi: &FluxGrid3D) -> VectorFluxGrid3D {
    let nx = phi.nx;
    let ny = phi.ny;
    let nz = phi.nz;
    let dx = phi.dx;

    let mut grad = VectorFluxGrid3D::zeros(nx, ny, nz, dx);

    let u_vals: Vec<f64> = (0..nx * ny * nz)
        .into_par_iter()
        .map(|idx| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let ii = i as isize;
            let jj = j as isize;
            let kk = k as isize;
            (phi.get_clamped(ii + 1, jj, kk) - phi.get_clamped(ii - 1, jj, kk)) / (2.0 * dx)
        })
        .collect();

    let v_vals: Vec<f64> = (0..nx * ny * nz)
        .into_par_iter()
        .map(|idx| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let ii = i as isize;
            let jj = j as isize;
            let kk = k as isize;
            (phi.get_clamped(ii, jj + 1, kk) - phi.get_clamped(ii, jj - 1, kk)) / (2.0 * dx)
        })
        .collect();

    let w_vals: Vec<f64> = (0..nx * ny * nz)
        .into_par_iter()
        .map(|idx| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let ii = i as isize;
            let jj = j as isize;
            let kk = k as isize;
            (phi.get_clamped(ii, jj, kk + 1) - phi.get_clamped(ii, jj, kk - 1)) / (2.0 * dx)
        })
        .collect();

    grad.u.values = u_vals;
    grad.v.values = v_vals;
    grad.w.values = w_vals;
    grad
}

// ---------------------------------------------------------------------------
// 7. Finite-volume Euler step helper
// ---------------------------------------------------------------------------

/// Advance `phi` by one Euler time step using explicit upwind advection.
///
/// Convenience wrapper: `phi_{n+1} = advect_upwind_3d(phi_n, vel, dt)`.
pub fn euler_step_advect(phi: &FluxGrid3D, vel: &VectorFluxGrid3D, dt: f64) -> FluxGrid3D {
    advect_upwind_3d(phi, vel, dt)
}

/// Compute the CFL-limited maximum stable time step.
///
/// `dt ≤ CFL_factor * dx / max(|u|, |v|, |w|)`.
pub fn cfl_dt(vel: &VectorFluxGrid3D, cfl_factor: f64) -> f64 {
    let max_u = vel
        .u
        .values
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let max_v = vel
        .v
        .values
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let max_w = vel
        .w
        .values
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let max_vel = max_u.max(max_v).max(max_w);
    let dx = vel.u.dx;
    if max_vel < 1e-300 {
        f64::INFINITY
    } else {
        cfl_factor * dx / max_vel
    }
}

// ---------------------------------------------------------------------------
// 8. Conservative state vector for Euler equations
// ---------------------------------------------------------------------------

/// Conservative variable vector for the Euler equations: \[rho, rho*u, rho*v, rho*w, E\].
#[derive(Debug, Clone, Copy)]
pub struct EulerState {
    /// Density.
    pub rho: f64,
    /// x-momentum.
    pub rho_u: f64,
    /// y-momentum.
    pub rho_v: f64,
    /// z-momentum.
    pub rho_w: f64,
    /// Total energy density.
    pub e: f64,
}

impl EulerState {
    /// Create from primitives `(rho, u, v, w, p)` with ratio of specific heats `gamma`.
    pub fn from_primitives(rho: f64, u: f64, v: f64, w: f64, p: f64, gamma: f64) -> Self {
        let ke = 0.5 * rho * (u * u + v * v + w * w);
        let e = p / (gamma - 1.0) + ke;
        Self {
            rho,
            rho_u: rho * u,
            rho_v: rho * v,
            rho_w: rho * w,
            e,
        }
    }

    /// Reconstruct velocity components.
    pub fn velocity(&self) -> [f64; 3] {
        let rho = if self.rho.abs() > 1e-30 {
            self.rho
        } else {
            1e-30
        };
        [self.rho_u / rho, self.rho_v / rho, self.rho_w / rho]
    }

    /// Reconstruct pressure from conserved variables.
    pub fn pressure(&self, gamma: f64) -> f64 {
        let u = self.velocity();
        let ke = 0.5 * self.rho * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        (gamma - 1.0) * (self.e - ke)
    }

    /// Sound speed.
    pub fn sound_speed(&self, gamma: f64) -> f64 {
        let p = self.pressure(gamma);
        let rho = if self.rho.abs() > 1e-30 {
            self.rho
        } else {
            1e-30
        };
        (gamma * p / rho).sqrt().max(0.0)
    }

    /// x-direction Euler flux vector.
    pub fn flux_x(&self, gamma: f64) -> [f64; 5] {
        let u = self.velocity();
        let p = self.pressure(gamma);
        [
            self.rho_u,
            self.rho_u * u[0] + p,
            self.rho_v * u[0],
            self.rho_w * u[0],
            (self.e + p) * u[0],
        ]
    }

    /// Linear interpolation between two states.
    pub fn lerp(&self, other: &EulerState, t: f64) -> EulerState {
        EulerState {
            rho: self.rho + t * (other.rho - self.rho),
            rho_u: self.rho_u + t * (other.rho_u - self.rho_u),
            rho_v: self.rho_v + t * (other.rho_v - self.rho_v),
            rho_w: self.rho_w + t * (other.rho_w - self.rho_w),
            e: self.e + t * (other.e - self.e),
        }
    }
}

// ---------------------------------------------------------------------------
// 9. Van Albada slope limiter
// ---------------------------------------------------------------------------

/// Van Albada slope limiter: smooth TVD limiter with quadratic shape.
///
/// `φ(r) = (r² + r) / (r² + 1)` for `r = a/b`.
pub fn van_albada(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        return 0.0;
    }
    let r = a / b;
    b * (r * r + r) / (r * r + 1.0)
}

// ---------------------------------------------------------------------------
// 10. MUSCL reconstruction (piecewise-linear, limited)
// ---------------------------------------------------------------------------

/// Reconstruct left and right interface values from cell-centred values
/// using the MUSCL scheme with a given limiter.
///
/// Returns `(u_L, u_R)` at the interface between cells `i` and `i+1`.
///
/// Uses: `u_L = u[i]   + 0.5 * phi(delta_m, delta_p) * delta_m`
///       `u_R = u[i+1] - 0.5 * phi(delta_p, delta_m) * delta_p`
/// where `delta_m = u[i] - u[i-1]`, `delta_p = u[i+1] - u[i]`.
pub fn muscl_reconstruct(u: &[f64], i: usize, limiter: fn(f64, f64) -> f64) -> (f64, f64) {
    let n = u.len();
    let u_im = if i == 0 { u[0] } else { u[i - 1] };
    let u_i = u[i];
    let u_ip = if i + 1 < n { u[i + 1] } else { u[n - 1] };
    let u_ipp = if i + 2 < n { u[i + 2] } else { u[n - 1] };

    let delta_m = u_i - u_im;
    let delta_p = u_ip - u_i;
    let delta_pp = u_ipp - u_ip;

    let sigma_l = limiter(delta_m, delta_p);
    let sigma_r = limiter(delta_p, delta_pp);

    let u_l = u_i + 0.5 * sigma_l;
    let u_r = u_ip - 0.5 * sigma_r;
    (u_l, u_r)
}

/// Apply MUSCL reconstruction across all interfaces in a 1-D array.
///
/// Returns a `Vec<(f64, f64)>` of `(u_L, u_R)` at interfaces `0..n-1`.
pub fn muscl_reconstruct_all(u: &[f64], limiter: fn(f64, f64) -> f64) -> Vec<(f64, f64)> {
    let n = u.len();
    if n < 2 {
        return Vec::new();
    }
    (0..n - 1)
        .map(|i| muscl_reconstruct(u, i, limiter))
        .collect()
}

// ---------------------------------------------------------------------------
// 11. Godunov flux (exact Riemann solver for linear advection)
// ---------------------------------------------------------------------------

/// Godunov flux for a scalar advection equation `u_t + a * u_x = 0`.
///
/// Uses upwinding: the Godunov flux is `F = a * u_L` if `a >= 0`, else `a * u_R`.
pub fn godunov_flux_advection(u_l: f64, u_r: f64, wave_speed: f64) -> f64 {
    if wave_speed >= 0.0 {
        wave_speed * u_l
    } else {
        wave_speed * u_r
    }
}

/// Godunov-type flux for Burgers' equation `u_t + (u²/2)_x = 0`.
///
/// Handles the sonic point (where the characteristic speed changes sign).
pub fn godunov_flux_burgers(u_l: f64, u_r: f64) -> f64 {
    // Rankine-Hugoniot shock speed: s = (u_l + u_r) / 2
    // Entropy fix: rarefaction fans may straddle u=0
    if u_l >= u_r {
        // Shock
        let s = 0.5 * (u_l + u_r);
        if s >= 0.0 {
            0.5 * u_l * u_l
        } else {
            0.5 * u_r * u_r
        }
    } else {
        // Rarefaction
        if u_l >= 0.0 {
            0.5 * u_l * u_l
        } else if u_r <= 0.0 {
            0.5 * u_r * u_r
        } else {
            0.0 // entropy-satisfying: sonic point at u=0
        }
    }
}

// ---------------------------------------------------------------------------
// 12. Roe flux for the linear-advection system
// ---------------------------------------------------------------------------

/// Roe-averaged scalar flux for a 1-D scalar conservation law.
///
/// For `u_t + f(u)_x = 0` with `f(u) = a * u`:
/// `F_Roe = 0.5*(f_L + f_R) - 0.5*|a_Roe|*(u_R - u_L)`
pub fn roe_flux_scalar(u_l: f64, u_r: f64, a_roe: f64) -> f64 {
    let f_l = a_roe * u_l;
    let f_r = a_roe * u_r;
    0.5 * (f_l + f_r) - 0.5 * a_roe.abs() * (u_r - u_l)
}

/// Roe flux for 1-D Euler equations (x-direction).
///
/// Returns the interface flux `[F_rho, F_rhou, F_e]` using the Roe-average state.
/// `gamma` is the ratio of specific heats (1.4 for air).
pub fn roe_flux_euler_1d(
    rho_l: f64,
    u_l: f64,
    p_l: f64,
    rho_r: f64,
    u_r: f64,
    p_r: f64,
    gamma: f64,
) -> [f64; 3] {
    // Roe averages
    let sqrt_rl = rho_l.sqrt();
    let sqrt_rr = rho_r.sqrt();
    let denom = sqrt_rl + sqrt_rr;

    let u_roe = (sqrt_rl * u_l + sqrt_rr * u_r) / denom;
    let h_l = (p_l / (rho_l * (gamma - 1.0)) + p_l / rho_l) + 0.5 * u_l * u_l;
    let h_r = (p_r / (rho_r * (gamma - 1.0)) + p_r / rho_r) + 0.5 * u_r * u_r;
    let h_roe = (sqrt_rl * h_l + sqrt_rr * h_r) / denom;

    let a2 = (gamma - 1.0) * (h_roe - 0.5 * u_roe * u_roe);
    let a_roe = if a2 > 0.0 { a2.sqrt() } else { 0.0 };

    // Conservative variables
    let e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * u_l * u_l;
    let e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * u_r * u_r;

    // Physical fluxes
    let fl = [rho_l * u_l, rho_l * u_l * u_l + p_l, (e_l + p_l) * u_l];
    let fr = [rho_r * u_r, rho_r * u_r * u_r + p_r, (e_r + p_r) * u_r];

    // Delta conserved variables
    let du = [rho_r - rho_l, rho_r * u_r - rho_l * u_l, e_r - e_l];

    // Eigenvalues
    let lam = [u_roe - a_roe, u_roe, u_roe + a_roe];

    // Roe dissipation matrix applied to du (simplified diagonal form)
    let r_inv_du = [
        (du[0] * u_roe - du[1]) / (2.0 * a_roe.max(1e-30)),
        du[0] - (du[0] * u_roe - du[1]) / (a_roe.max(1e-30)),
        (du[0] * u_roe + du[1]) / (2.0 * a_roe.max(1e-30)),
    ];

    // Right eigenvectors (simplified form)
    let r_mat = [
        [1.0, 1.0, 1.0],
        [u_roe - a_roe, u_roe, u_roe + a_roe],
        [
            h_roe - u_roe * a_roe,
            0.5 * u_roe * u_roe,
            h_roe + u_roe * a_roe,
        ],
    ];

    let mut dissipation = [0.0f64; 3];
    for i in 0..3 {
        for k in 0..3 {
            dissipation[i] += lam[k].abs() * r_inv_du[k] * r_mat[i][k];
        }
    }

    [
        0.5 * (fl[0] + fr[0]) - 0.5 * dissipation[0],
        0.5 * (fl[1] + fr[1]) - 0.5 * dissipation[1],
        0.5 * (fl[2] + fr[2]) - 0.5 * dissipation[2],
    ]
}

// ---------------------------------------------------------------------------
// 13. HLL flux
// ---------------------------------------------------------------------------

/// HLL (Harten-Lax-van Leer) numerical flux for a scalar conservation law.
///
/// Wave speed estimates `s_l` and `s_r` are the left and right running speeds.
/// `f_l = f(u_l)`, `f_r = f(u_r)`.
pub fn hll_flux(u_l: f64, u_r: f64, f_l: f64, f_r: f64, s_l: f64, s_r: f64) -> f64 {
    if s_l >= 0.0 {
        f_l
    } else if s_r <= 0.0 {
        f_r
    } else {
        (s_r * f_l - s_l * f_r + s_l * s_r * (u_r - u_l)) / (s_r - s_l)
    }
}

/// Estimate HLL wave speeds for the 1-D Euler equations.
///
/// Uses Davis' estimates: `s_l = min(u_l - a_l, u_r - a_r)`, `s_r = max(u_l + a_l, u_r + a_r)`.
pub fn hll_wave_speeds(u_l: f64, a_l: f64, u_r: f64, a_r: f64) -> (f64, f64) {
    let s_l = (u_l - a_l).min(u_r - a_r);
    let s_r = (u_l + a_l).max(u_r + a_r);
    (s_l, s_r)
}

/// HLL flux for 1-D Euler equations.
///
/// Returns `[F_rho, F_rhou, F_e]`.
pub fn hll_flux_euler_1d(
    rho_l: f64,
    u_l: f64,
    p_l: f64,
    rho_r: f64,
    u_r: f64,
    p_r: f64,
    gamma: f64,
) -> [f64; 3] {
    let a_l = (gamma * p_l / rho_l).sqrt().max(0.0);
    let a_r = (gamma * p_r / rho_r).sqrt().max(0.0);
    let (s_l, s_r) = hll_wave_speeds(u_l, a_l, u_r, a_r);

    let e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * u_l * u_l;
    let e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * u_r * u_r;

    let fl = [rho_l * u_l, rho_l * u_l * u_l + p_l, (e_l + p_l) * u_l];
    let fr = [rho_r * u_r, rho_r * u_r * u_r + p_r, (e_r + p_r) * u_r];
    let ul = [rho_l, rho_l * u_l, e_l];
    let ur = [rho_r, rho_r * u_r, e_r];

    let mut f = [0.0f64; 3];
    for i in 0..3 {
        f[i] = hll_flux(ul[i], ur[i], fl[i], fr[i], s_l, s_r);
    }
    f
}

// ---------------------------------------------------------------------------
// 14. HLLC flux for 1-D Euler equations
// ---------------------------------------------------------------------------

/// HLLC (Harten-Lax-van Leer-Contact) flux for 1-D Euler equations.
///
/// The HLLC flux restores the contact wave missing in HLL.
/// Returns `[F_rho, F_rhou, F_e]`.
pub fn hllc_flux_euler_1d(
    rho_l: f64,
    u_l: f64,
    p_l: f64,
    rho_r: f64,
    u_r: f64,
    p_r: f64,
    gamma: f64,
) -> [f64; 3] {
    let a_l = (gamma * p_l / rho_l.max(1e-30)).sqrt();
    let a_r = (gamma * p_r / rho_r.max(1e-30)).sqrt();

    // Davis wave speed estimates
    let (s_l, s_r) = hll_wave_speeds(u_l, a_l, u_r, a_r);

    // Contact wave speed (HLLC star state)
    let s_star = (p_r - p_l + rho_l * u_l * (s_l - u_l) - rho_r * u_r * (s_r - u_r))
        / (rho_l * (s_l - u_l) - rho_r * (s_r - u_r) + 1e-200);

    let e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * u_l * u_l;
    let e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * u_r * u_r;

    let fl = [rho_l * u_l, rho_l * u_l * u_l + p_l, (e_l + p_l) * u_l];
    let fr = [rho_r * u_r, rho_r * u_r * u_r + p_r, (e_r + p_r) * u_r];

    let hllc_flux_side = |rho: f64, u: f64, e: f64, p: f64, f: &[f64; 3], s: f64| -> [f64; 3] {
        let coeff = rho * (s - u) / (s - s_star + 1e-200);
        let u_star = [
            coeff,
            coeff * s_star,
            coeff * (e / rho + (s_star - u) * (s_star + p / (rho * (s - u) + 1e-200))),
        ];
        [
            f[0] + s * (u_star[0] - rho),
            f[1] + s * (u_star[1] - rho * u),
            f[2] + s * (u_star[2] - e),
        ]
    };

    if s_l >= 0.0 {
        fl
    } else if s_star >= 0.0 {
        hllc_flux_side(rho_l, u_l, e_l, p_l, &fl, s_l)
    } else if s_r >= 0.0 {
        hllc_flux_side(rho_r, u_r, e_r, p_r, &fr, s_r)
    } else {
        fr
    }
}

// ---------------------------------------------------------------------------
// 15. Characteristic decomposition for 1-D advection systems
// ---------------------------------------------------------------------------

/// Project a vector of conservative variables onto the characteristic variables
/// for a simple 2-component advection system.
///
/// Eigenvalues `lambda = [a, -a]` (rightward and leftward waves).
/// Returns `(w_plus, w_minus)` — the characteristic amplitudes.
pub fn characteristic_decompose_2wave(u0: f64, u1: f64, a: f64) -> (f64, f64) {
    let w_plus = 0.5 * (u0 + u1 / a);
    let w_minus = 0.5 * (u0 - u1 / a);
    (w_plus, w_minus)
}

/// Reconstruct conservative variables from characteristic variables.
pub fn characteristic_recompose_2wave(w_plus: f64, w_minus: f64, a: f64) -> (f64, f64) {
    let u0 = w_plus + w_minus;
    let u1 = a * (w_plus - w_minus);
    (u0, u1)
}

/// Apply a characteristic-based limiter to 1-D reconstruction.
///
/// Projects to characteristic variables, limits each wave separately, then
/// reconstructs in physical space.  This is important for multi-component
/// systems to prevent spurious oscillations across contact waves.
pub fn characteristic_limited_reconstruct(
    u0: &[f64],
    u1: &[f64],
    i: usize,
    a: f64,
    limiter: fn(f64, f64) -> f64,
) -> ([f64; 2], [f64; 2]) {
    let n = u0.len();
    let get = |arr: &[f64], j: usize| -> f64 { arr[j.min(n - 1)] };

    // Project cell values around i to characteristic space
    let (wm_im, _) = characteristic_decompose_2wave(
        get(u0, i.saturating_sub(1)),
        get(u1, i.saturating_sub(1)),
        a,
    );
    let (wm_i, wp_i) = characteristic_decompose_2wave(get(u0, i), get(u1, i), a);
    let (wm_ip, wp_ip) = characteristic_decompose_2wave(get(u0, i + 1), get(u1, i + 1), a);
    let (wm_ipp, wp_ipp) = if i + 2 < n {
        characteristic_decompose_2wave(get(u0, i + 2), get(u1, i + 2), a)
    } else {
        (wm_ip, wp_ip)
    };
    let (_, wp_im) = characteristic_decompose_2wave(
        get(u0, i.saturating_sub(1)),
        get(u1, i.saturating_sub(1)),
        a,
    );

    // Limit each characteristic wave
    let sig_p_l = limiter(wp_i - wp_im, wp_ip - wp_i);
    let sig_p_r = limiter(wp_ip - wp_i, wp_ipp - wp_ip);
    let sig_m_l = limiter(wm_i - wm_im, wm_ip - wm_i);
    let sig_m_r = limiter(wm_ip - wm_i, wm_ipp - wm_ip);

    let wp_l = wp_i + 0.5 * sig_p_l;
    let wp_r = wp_ip - 0.5 * sig_p_r;
    let wm_l = wm_i + 0.5 * sig_m_l;
    let wm_r = wm_ip - 0.5 * sig_m_r;

    let (u0_l, u1_l) = characteristic_recompose_2wave(wp_l, wm_l, a);
    let (u0_r, u1_r) = characteristic_recompose_2wave(wp_r, wm_r, a);

    ([u0_l, u1_l], [u0_r, u1_r])
}

// ---------------------------------------------------------------------------
// 16. Second-order Runge-Kutta (TVD RK2) time stepping
// ---------------------------------------------------------------------------

/// Advance `phi` using the TVD Runge-Kutta 2 (Heun) scheme.
///
/// `phi_{n+1} = 0.5 * (phi_n + phi_n^(2))`
/// `phi_n^(1) = phi_n + dt * L(phi_n)` (Euler predictor)
/// `phi_n^(2) = phi_n^(1) + dt * L(phi_n^(1))` (corrector)
///
/// where `L` is the upwind advection operator.
pub fn tvd_rk2_advect(phi: &FluxGrid3D, vel: &VectorFluxGrid3D, dt: f64) -> FluxGrid3D {
    // Stage 1: phi^(1) = phi + dt * L(phi)
    let phi_1 = advect_upwind_3d(phi, vel, dt);
    // Stage 2: phi^(2) = phi^(1) + dt * L(phi^(1))
    let phi_2 = advect_upwind_3d(&phi_1, vel, dt);
    // Combine: phi_{n+1} = 0.5*(phi + phi^(2))
    let mut result = phi.clone();
    result
        .values
        .iter_mut()
        .zip(phi.values.iter().zip(phi_2.values.iter()))
        .for_each(|(r, (&a, &b))| *r = 0.5 * (a + b));
    result
}

// ---------------------------------------------------------------------------
// 17. Flux limiter in ratio form
// ---------------------------------------------------------------------------

/// Compute the minmod limiter in ratio form: `phi(r) = max(0, min(1, r))`.
///
/// Input `r = delta_{i-1/2} / delta_{i+1/2}` (ratio of consecutive differences).
pub fn minmod_ratio(r: f64) -> f64 {
    if r <= 0.0 {
        0.0
    } else if r <= 1.0 {
        r
    } else {
        1.0
    }
}

/// Superbee limiter in ratio form: `phi(r) = max(0, max(min(2r,1), min(r,2)))`.
pub fn superbee_ratio(r: f64) -> f64 {
    if r <= 0.0 {
        0.0
    } else {
        (2.0 * r).min(1.0_f64).max(r.min(2.0_f64)).max(0.0_f64)
    }
}

/// Van Leer limiter in ratio form: `phi(r) = (r + |r|) / (1 + |r|)`.
pub fn van_leer_ratio(r: f64) -> f64 {
    (r + r.abs()) / (1.0 + r.abs())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod flux_compute_tests {
    use super::*;

    #[test]
    fn test_flux_grid3d_get_set() {
        let mut g = FluxGrid3D::new(4, 4, 4, 0.1, 0.0);
        g.set(1, 2, 3, 7.5);
        assert!((g.get(1, 2, 3) - 7.5).abs() < 1e-12);
    }

    #[test]
    fn test_flux_grid3d_clamped_boundary() {
        let g = FluxGrid3D::new(5, 5, 5, 0.1, 3.0);
        // Negative indices and out-of-bounds indices should clamp.
        assert!((g.get_clamped(-1, 0, 0) - 3.0).abs() < 1e-12);
        assert!((g.get_clamped(100, 0, 0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_l2_norm_uniform() {
        // A 2×2×2 grid filled with 1.0: l2 = sqrt(8).
        let g = FluxGrid3D::new(2, 2, 2, 0.5, 1.0);
        let expected = (8.0_f64).sqrt();
        assert!((g.l2_norm() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_minmod_same_sign() {
        assert!(
            (minmod(2.0, 3.0) - 2.0).abs() < 1e-12,
            "minmod picks smaller"
        );
        assert!((minmod(-1.0, -4.0) - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_minmod_opposite_sign() {
        assert!(
            (minmod(2.0, -1.0)).abs() < 1e-12,
            "minmod returns 0 for opposite signs"
        );
    }

    #[test]
    fn test_van_leer_smooth() {
        // van_leer(a, b) = 2ab/(a+b) for same-sign inputs.
        let vl = van_leer(2.0, 2.0);
        assert!(
            (vl - 2.0).abs() < 1e-12,
            "symmetric inputs: van_leer = value"
        );
        let vl2 = van_leer(1.0, 3.0);
        let expected = 2.0 * 1.0 * 3.0 / (1.0 + 3.0);
        assert!((vl2 - expected).abs() < 1e-12);
    }

    #[test]
    fn test_upwind_flux_1d_positive_velocity() {
        // Uniform phi = 1.0, positive velocity: no flux.
        let phi = vec![1.0; 10];
        let flux = upwind_flux_1d(&phi, 1.0, 0.1, 0.05);
        for f in &flux {
            assert!(
                f.abs() < 1e-12,
                "uniform field with positive vel: zero flux"
            );
        }
    }

    #[test]
    fn test_lax_friedrichs_1d_preserves_mass() {
        // Total mass ≈ integral of u dx. Advection should conserve mass
        // (up to boundary corrections). Check sum is approximately preserved.
        let u: Vec<f64> = (0..20).map(|i| (i as f64 * 0.1 - 1.0).powi(2)).collect();
        let sum_before: f64 = u.iter().sum();
        let u_new = lax_friedrichs_advect_1d(&u, 0.5, 0.1, 0.01);
        let sum_after: f64 = u_new.iter().sum();
        // Small boundary effect allowed.
        assert!(
            (sum_after - sum_before).abs() < sum_before * 0.05 + 1.0,
            "mass should be approximately conserved"
        );
    }

    #[test]
    fn test_divergence_constant_field_is_zero() {
        // Uniform velocity field → divergence = 0 everywhere.
        let n = 6;
        let mut vel = VectorFluxGrid3D::zeros(n, n, n, 0.1);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    vel.set_vel(i, j, k, [1.0, 2.0, 3.0]);
                }
            }
        }
        let div = divergence_3d(&vel);
        for &d in &div.values {
            assert!(
                d.abs() < 1e-10,
                "divergence of uniform field must be 0, got {d}"
            );
        }
    }

    #[test]
    fn test_gradient_linear_field() {
        // phi(i,j,k) = x = i*dx. Gradient_x should be ~1.0 everywhere.
        let n = 8;
        let dx = 0.25;
        let mut phi = FluxGrid3D::new(n, n, n, dx, 0.0);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    phi.set(i, j, k, i as f64 * dx);
                }
            }
        }
        let grad = gradient_3d(&phi);
        // Check interior cells.
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                for k in 1..n - 1 {
                    let gx = grad.u.get(i, j, k);
                    assert!(
                        (gx - 1.0).abs() < 1e-10,
                        "gradient_x of linear field should be 1, got {gx} at ({i},{j},{k})"
                    );
                }
            }
        }
    }

    #[test]
    fn test_cfl_dt_uniform_velocity() {
        let mut vel = VectorFluxGrid3D::zeros(4, 4, 4, 0.1);
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    vel.set_vel(i, j, k, [2.0, 0.0, 0.0]);
                }
            }
        }
        // CFL dt = 0.5 * dx / 2.0 = 0.025
        let dt = cfl_dt(&vel, 0.5);
        assert!(
            (dt - 0.025).abs() < 1e-10,
            "cfl_dt should be 0.025, got {dt}"
        );
    }

    #[test]
    fn test_cfl_dt_zero_velocity() {
        let vel = VectorFluxGrid3D::zeros(4, 4, 4, 0.1);
        let dt = cfl_dt(&vel, 0.5);
        assert!(dt.is_infinite(), "zero velocity → infinite dt");
    }

    // ── Additional FluxGrid3D tests ──────────────────────────────────────────

    #[test]
    fn test_flux_grid3d_len_empty() {
        let g = FluxGrid3D::new(0, 5, 5, 0.1, 0.0);
        assert_eq!(g.len(), 0);
        assert!(g.is_empty());
    }

    #[test]
    fn test_flux_grid3d_len_nonempty() {
        let g = FluxGrid3D::new(3, 4, 5, 0.1, 0.0);
        assert_eq!(g.len(), 60);
        assert!(!g.is_empty());
    }

    #[test]
    fn test_flux_grid3d_max_val() {
        let mut g = FluxGrid3D::new(2, 2, 2, 0.1, 0.0);
        g.set(0, 0, 0, 5.0);
        g.set(1, 1, 1, -3.0);
        assert!((g.max_val() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_flux_grid3d_min_val() {
        let mut g = FluxGrid3D::new(2, 2, 2, 0.1, 0.0);
        g.set(0, 0, 0, 5.0);
        g.set(1, 1, 1, -3.0);
        assert!((g.min_val() - (-3.0)).abs() < 1e-12);
    }

    #[test]
    fn test_flux_grid3d_l2_norm_zeros() {
        let g = FluxGrid3D::new(3, 3, 3, 0.1, 0.0);
        assert!(g.l2_norm() < 1e-15);
    }

    #[test]
    fn test_flux_grid3d_index_layout() {
        let g = FluxGrid3D::new(3, 4, 5, 0.1, 0.0);
        // Check index formula: i*ny*nz + j*nz + k
        assert_eq!(g.idx(0, 0, 0), 0);
        assert_eq!(g.idx(1, 0, 0), 20); // 1 * 4 * 5
        assert_eq!(g.idx(0, 1, 0), 5); // 0 * 20 + 1 * 5
        assert_eq!(g.idx(0, 0, 1), 1);
    }

    #[test]
    fn test_flux_grid3d_clamped_low_indices() {
        let mut g = FluxGrid3D::new(5, 5, 5, 0.1, 0.0);
        g.set(0, 0, 0, 99.0);
        // Clamping -1 → 0
        assert!((g.get_clamped(-1, 0, 0) - 99.0).abs() < 1e-12);
        assert!((g.get_clamped(0, -5, 0) - 99.0).abs() < 1e-12);
    }

    // ── VectorFluxGrid3D tests ───────────────────────────────────────────────

    #[test]
    fn test_vector_flux_grid_zeros() {
        let v = VectorFluxGrid3D::zeros(3, 3, 3, 0.1);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let vel = v.get_vel(i, j, k);
                    assert_eq!(vel, [0.0, 0.0, 0.0]);
                }
            }
        }
    }

    #[test]
    fn test_vector_flux_grid_set_get_vel() {
        let mut v = VectorFluxGrid3D::zeros(4, 4, 4, 0.1);
        v.set_vel(1, 2, 3, [1.5, 2.5, 3.5]);
        let vel = v.get_vel(1, 2, 3);
        assert!((vel[0] - 1.5).abs() < 1e-12);
        assert!((vel[1] - 2.5).abs() < 1e-12);
        assert!((vel[2] - 3.5).abs() < 1e-12);
    }

    // ── Limiter tests ────────────────────────────────────────────────────────

    #[test]
    fn test_superbee_same_sign_positive() {
        // superbee(a, b) = max(min(2a, b), min(a, 2b)) for same-sign
        let s = superbee(1.0, 2.0);
        assert!(s > 0.0, "superbee of positives should be positive, got {s}");
    }

    #[test]
    fn test_superbee_opposite_signs() {
        assert_eq!(superbee(1.0, -1.0), 0.0);
        assert_eq!(superbee(-2.0, 3.0), 0.0);
    }

    #[test]
    fn test_superbee_both_zero() {
        assert_eq!(superbee(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_mc_limiter_same_sign() {
        let mc = mc_limiter(2.0, 4.0);
        // centered = 3, limited by 2*min(2,4) = 4, min(2*2, 2*4) = 4, answer ≤ 3
        assert!(mc > 0.0 && mc <= 3.0, "mc_limiter = {mc}");
    }

    #[test]
    fn test_mc_limiter_opposite_sign() {
        assert_eq!(mc_limiter(1.0, -1.0), 0.0);
    }

    #[test]
    fn test_van_leer_zero_sum() {
        // When a + b ≈ 0 this should handle gracefully
        // van_leer(a, -a) = 0 because a * -a < 0
        let vl = van_leer(2.0, -2.0);
        assert_eq!(vl, 0.0);
    }

    #[test]
    fn test_minmod_equal_values() {
        // minmod(a, a) = a
        assert!((minmod(3.0, 3.0) - 3.0).abs() < 1e-12);
        assert!((minmod(-2.0, -2.0) - (-2.0)).abs() < 1e-12);
    }

    // ── Upwind flux tests ────────────────────────────────────────────────────

    #[test]
    fn test_upwind_flux_1d_negative_velocity() {
        // Uniform phi = 1.0, negative velocity: no flux.
        let phi = vec![1.0; 8];
        let flux = upwind_flux_1d(&phi, -1.0, 0.1, 0.05);
        for f in &flux {
            assert!(
                f.abs() < 1e-12,
                "uniform field with negative vel: zero flux"
            );
        }
    }

    #[test]
    fn test_upwind_flux_1d_step_function() {
        // Step function: 0 left, 1 right. Positive velocity should transport left.
        let phi = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let flux = upwind_flux_1d(&phi, 1.0, 0.1, 0.01);
        // At boundary cell (index 3), phi_l = 0, phi_c = 1 → flux negative
        assert!(
            flux[3] < 0.0,
            "step up with positive vel: negative flux at boundary"
        );
    }

    #[test]
    fn test_upwind_flux_1d_length() {
        let phi = vec![1.0; 7];
        let flux = upwind_flux_1d(&phi, 1.0, 0.1, 0.01);
        assert_eq!(flux.len(), 7);
    }

    // ── Lax-Friedrichs tests ─────────────────────────────────────────────────

    #[test]
    fn test_lax_friedrichs_flux_zero_difference() {
        // u_l == u_r, f_l == f_r: LF flux = f_l = f_r
        let f = lax_friedrichs_flux(1.0, 1.0, 2.0, 2.0, 1.0);
        assert!((f - 2.0).abs() < 1e-12, "LF flux = {f}");
    }

    #[test]
    fn test_lax_friedrichs_advect_1d_empty() {
        let result = lax_friedrichs_advect_1d(&[], 1.0, 0.1, 0.01);
        assert!(result.is_empty());
    }

    #[test]
    fn test_lax_friedrichs_advect_1d_length() {
        let u = vec![1.0; 6];
        let u_new = lax_friedrichs_advect_1d(&u, 0.5, 0.1, 0.01);
        assert_eq!(u_new.len(), 6);
    }

    #[test]
    fn test_lax_friedrichs_uniform_field_stable() {
        // Uniform field should not change much.
        let u = vec![1.0; 10];
        let u_new = lax_friedrichs_advect_1d(&u, 1.0, 0.1, 0.05);
        let diff: f64 = u_new.iter().zip(u.iter()).map(|(a, b)| (a - b).abs()).sum();
        assert!(
            diff < 1e-10,
            "uniform field should not change under LF advection"
        );
    }

    // ── Divergence / gradient tests ──────────────────────────────────────────

    #[test]
    fn test_divergence_size_matches_input() {
        let vel = VectorFluxGrid3D::zeros(3, 4, 5, 0.1);
        let div = divergence_3d(&vel);
        assert_eq!(div.nx, 3);
        assert_eq!(div.ny, 4);
        assert_eq!(div.nz, 5);
    }

    #[test]
    fn test_gradient_size_matches_input() {
        let phi = FluxGrid3D::new(3, 4, 5, 0.1, 0.0);
        let grad = gradient_3d(&phi);
        assert_eq!(grad.u.nx, 3);
        assert_eq!(grad.v.ny, 4);
        assert_eq!(grad.w.nz, 5);
    }

    #[test]
    fn test_gradient_constant_field_is_zero() {
        let phi = FluxGrid3D::new(5, 5, 5, 0.1, 7.0);
        let grad = gradient_3d(&phi);
        for &g in &grad.u.values {
            assert!(g.abs() < 1e-10, "gradient of constant = 0");
        }
    }

    #[test]
    fn test_divergence_linear_vx_field() {
        // u(i,j,k) = i*dx, v=w=0 → div = ∂u/∂x = 1 at interior cells
        let n = 8;
        let dx = 0.5;
        let mut vel = VectorFluxGrid3D::zeros(n, n, n, dx);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    vel.u.set(i, j, k, i as f64 * dx);
                }
            }
        }
        let div = divergence_3d(&vel);
        // Interior cells: central diff of i*dx → 1
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                for k in 1..n - 1 {
                    let d = div.get(i, j, k);
                    assert!(
                        (d - 1.0).abs() < 1e-10,
                        "div at ({i},{j},{k}) = {d}, expected 1"
                    );
                }
            }
        }
    }

    // ── Euler step test ──────────────────────────────────────────────────────

    #[test]
    fn test_euler_step_advect_preserves_size() {
        let phi = FluxGrid3D::new(4, 4, 4, 0.1, 1.0);
        let mut vel = VectorFluxGrid3D::zeros(4, 4, 4, 0.1);
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    vel.set_vel(i, j, k, [1.0, 0.0, 0.0]);
                }
            }
        }
        let phi_new = euler_step_advect(&phi, &vel, 0.01);
        assert_eq!(phi_new.nx, phi.nx);
        assert_eq!(phi_new.ny, phi.ny);
        assert_eq!(phi_new.nz, phi.nz);
    }

    #[test]
    fn test_euler_step_constant_phi_no_change() {
        // Uniform phi, uniform velocity → upwind flux is zero → phi unchanged
        let phi = FluxGrid3D::new(5, 5, 5, 0.1, 3.0);
        let mut vel = VectorFluxGrid3D::zeros(5, 5, 5, 0.1);
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    vel.set_vel(i, j, k, [1.0, 1.0, 1.0]);
                }
            }
        }
        let phi_new = euler_step_advect(&phi, &vel, 0.01);
        for (&a, &b) in phi.values.iter().zip(phi_new.values.iter()) {
            assert!(
                (a - b).abs() < 1e-10,
                "constant phi should not change: {a} vs {b}"
            );
        }
    }

    // ── CFL dt tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_cfl_dt_negative_velocity() {
        // Negative velocity magnitude should still work.
        let mut vel = VectorFluxGrid3D::zeros(4, 4, 4, 0.2);
        vel.set_vel(0, 0, 0, [0.0, -4.0, 0.0]);
        let dt = cfl_dt(&vel, 1.0);
        assert!((dt - 0.05).abs() < 1e-10, "cfl_dt = {dt}, expected 0.05");
    }

    #[test]
    fn test_cfl_dt_multiple_components() {
        let mut vel = VectorFluxGrid3D::zeros(3, 3, 3, 0.1);
        // Set u=1, v=2, w=3 somewhere — max is 3
        vel.set_vel(1, 1, 1, [1.0, 2.0, 3.0]);
        let dt = cfl_dt(&vel, 1.0);
        assert!((dt - 0.1 / 3.0).abs() < 1e-10, "cfl_dt = {dt}");
    }

    // ── EulerState tests ─────────────────────────────────────────────────────

    #[test]
    fn test_euler_state_from_primitives() {
        let gamma = 1.4;
        let s = EulerState::from_primitives(1.0, 0.1, 0.0, 0.0, 1.0, gamma);
        assert!((s.rho - 1.0).abs() < 1e-12);
        let p_back = s.pressure(gamma);
        assert!((p_back - 1.0).abs() < 1e-10, "pressure roundtrip: {p_back}");
    }

    #[test]
    fn test_euler_state_velocity() {
        let gamma = 1.4;
        let s = EulerState::from_primitives(2.0, 3.0, 4.0, 5.0, 1.0, gamma);
        let u = s.velocity();
        assert!((u[0] - 3.0).abs() < 1e-12);
        assert!((u[1] - 4.0).abs() < 1e-12);
        assert!((u[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_euler_state_sound_speed() {
        let gamma = 1.4;
        // c = sqrt(gamma*p/rho) = sqrt(1.4 * 1.0 / 1.0) ≈ 1.1832
        let s = EulerState::from_primitives(1.0, 0.0, 0.0, 0.0, 1.0, gamma);
        let c = s.sound_speed(gamma);
        let expected = (gamma).sqrt();
        assert!(
            (c - expected).abs() < 1e-10,
            "sound speed = {c}, expected {expected}"
        );
    }

    #[test]
    fn test_euler_state_lerp() {
        let gamma = 1.4;
        let s0 = EulerState::from_primitives(1.0, 0.0, 0.0, 0.0, 1.0, gamma);
        let s1 = EulerState::from_primitives(2.0, 0.0, 0.0, 0.0, 2.0, gamma);
        let s05 = s0.lerp(&s1, 0.5);
        assert!((s05.rho - 1.5).abs() < 1e-12, "lerp rho = {}", s05.rho);
    }

    #[test]
    fn test_euler_state_flux_x() {
        let gamma = 1.4;
        // At rest: F = [0, p, 0, 0, 0]
        let s = EulerState::from_primitives(1.0, 0.0, 0.0, 0.0, 1.0, gamma);
        let f = s.flux_x(gamma);
        assert!(f[0].abs() < 1e-12, "mass flux at rest = {}", f[0]);
        assert!((f[1] - 1.0).abs() < 1e-10, "pressure flux = {}", f[1]);
    }

    // ── Van Albada tests ──────────────────────────────────────────────────────

    #[test]
    fn test_van_albada_same_sign() {
        // van_albada should return positive for same-sign inputs
        let v = van_albada(2.0, 4.0);
        assert!(v > 0.0, "van_albada(2,4) = {v}");
        // Should be TVD: |result| <= min(2|a|, 2|b|, |a+b|/2...)
        assert!(v <= 4.0);
    }

    #[test]
    fn test_van_albada_opposite_signs() {
        assert_eq!(van_albada(1.0, -1.0), 0.0);
        assert_eq!(van_albada(-2.0, 3.0), 0.0);
    }

    #[test]
    fn test_van_albada_symmetry() {
        // van_albada(a, b) vs van_albada(b, a): not strictly symmetric but both > 0
        let v1 = van_albada(1.0, 3.0);
        let v2 = van_albada(3.0, 1.0);
        assert!(v1 > 0.0 && v2 > 0.0);
    }

    // ── MUSCL tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_muscl_reconstruct_uniform() {
        // Uniform field: reconstruction should give exact values
        let u = vec![1.0_f64; 10];
        let (u_l, u_r) = muscl_reconstruct(&u, 3, minmod);
        assert!((u_l - 1.0).abs() < 1e-12, "uniform MUSCL u_l = {u_l}");
        assert!((u_r - 1.0).abs() < 1e-12, "uniform MUSCL u_r = {u_r}");
    }

    #[test]
    fn test_muscl_reconstruct_all_length() {
        let u = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let pairs = muscl_reconstruct_all(&u, minmod);
        assert_eq!(pairs.len(), 4, "MUSCL all: expected 4 interfaces");
    }

    #[test]
    fn test_muscl_reconstruct_all_monotone() {
        // Linear field: reconstruction should preserve ordering
        let u: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let pairs = muscl_reconstruct_all(&u, minmod);
        for (i, &(u_l, u_r)) in pairs.iter().enumerate() {
            assert!(u_l <= u_r + 1e-10, "interface {i}: u_l={u_l} > u_r={u_r}");
        }
    }

    // ── Godunov flux tests ────────────────────────────────────────────────────

    #[test]
    fn test_godunov_flux_advection_positive_speed() {
        // With a > 0, flux = a * u_l
        let f = godunov_flux_advection(2.0, 3.0, 1.0);
        assert!((f - 2.0).abs() < 1e-12, "Godunov advection +speed: {f}");
    }

    #[test]
    fn test_godunov_flux_advection_negative_speed() {
        let f = godunov_flux_advection(2.0, 3.0, -1.0);
        assert!((f - (-3.0)).abs() < 1e-12, "Godunov advection -speed: {f}");
    }

    #[test]
    fn test_godunov_flux_burgers_shock() {
        // Shock: u_l > u_r, s > 0 → f = 0.5 * u_l^2
        let f = godunov_flux_burgers(2.0, 1.0);
        assert!((f - 0.5 * 4.0).abs() < 1e-12, "Burgers shock: {f}");
    }

    #[test]
    fn test_godunov_flux_burgers_rarefaction_sonic() {
        // Rarefaction spanning u=0: u_l < 0 < u_r → entropy-satisfying flux = 0
        let f = godunov_flux_burgers(-1.0, 1.0);
        assert!(f.abs() < 1e-12, "Burgers sonic: {f}");
    }

    #[test]
    fn test_godunov_flux_burgers_rarefaction_positive() {
        // Rarefaction, all positive: u_l > 0, u_r > u_l → f = 0.5*u_l^2
        let f = godunov_flux_burgers(1.0, 2.0);
        assert!((f - 0.5).abs() < 1e-12, "Burgers pos rarefaction: {f}");
    }

    // ── Roe flux tests ────────────────────────────────────────────────────────

    #[test]
    fn test_roe_flux_scalar_symmetry() {
        // For uniform state, Roe flux = physical flux
        let f = roe_flux_scalar(2.0, 2.0, 1.5);
        assert!((f - 1.5 * 2.0).abs() < 1e-12, "Roe scalar uniform: {f}");
    }

    #[test]
    fn test_roe_flux_scalar_upwind() {
        // For a > 0 and jump, Roe flux should propagate from left
        let f = roe_flux_scalar(1.0, 2.0, 1.0);
        // f = 0.5*(1+2) - 0.5*1*(2-1) = 1.5 - 0.5 = 1.0
        assert!((f - 1.0).abs() < 1e-12, "Roe scalar upwind: {f}");
    }

    #[test]
    fn test_roe_flux_euler_1d_finite() {
        let f = roe_flux_euler_1d(1.0, 0.1, 1.0, 1.0, 0.1, 1.0, 1.4);
        for (i, &fi) in f.iter().enumerate() {
            assert!(fi.is_finite(), "Roe Euler flux[{i}] not finite: {fi}");
        }
    }

    // ── HLL flux tests ────────────────────────────────────────────────────────

    #[test]
    fn test_hll_flux_subsonic() {
        // s_l < 0 < s_r: active HLL formula
        let f = hll_flux(1.0, 2.0, 1.0, 2.0, -1.0, 1.0);
        assert!(f.is_finite(), "HLL subsonic: {f}");
    }

    #[test]
    fn test_hll_flux_supersonic_right() {
        // s_l > 0: use f_l
        let f = hll_flux(1.0, 2.0, 3.0, 4.0, 1.0, 2.0);
        assert!((f - 3.0).abs() < 1e-12, "HLL supersonic right: {f}");
    }

    #[test]
    fn test_hll_flux_supersonic_left() {
        // s_r < 0: use f_r
        let f = hll_flux(1.0, 2.0, 3.0, 4.0, -2.0, -1.0);
        assert!((f - 4.0).abs() < 1e-12, "HLL supersonic left: {f}");
    }

    #[test]
    fn test_hll_flux_euler_1d_finite() {
        let f = hll_flux_euler_1d(1.0, 0.1, 1.0, 1.2, 0.05, 1.1, 1.4);
        for &fi in &f {
            assert!(fi.is_finite(), "HLL Euler flux not finite");
        }
    }

    #[test]
    fn test_hll_wave_speeds() {
        let (s_l, s_r) = hll_wave_speeds(0.5, 1.0, -0.5, 1.0);
        assert!(s_l < 0.0 && s_r > 0.0, "wave speeds s_l={s_l} s_r={s_r}");
        assert!(s_l <= s_r);
    }

    // ── HLLC flux tests ───────────────────────────────────────────────────────

    #[test]
    fn test_hllc_flux_euler_1d_finite() {
        let f = hllc_flux_euler_1d(1.0, 0.1, 1.0, 1.2, 0.05, 1.1, 1.4);
        for &fi in &f {
            assert!(fi.is_finite(), "HLLC Euler flux not finite: {fi}");
        }
    }

    #[test]
    fn test_hllc_agrees_with_hll_at_uniform_state() {
        // For uniform state (no jump), HLL and HLLC should give similar answers
        let f_hll = hll_flux_euler_1d(1.0, 0.1, 1.0, 1.0, 0.1, 1.0, 1.4);
        let f_hllc = hllc_flux_euler_1d(1.0, 0.1, 1.0, 1.0, 0.1, 1.0, 1.4);
        for i in 0..3 {
            assert!(
                (f_hll[i] - f_hllc[i]).abs() < 1e-6,
                "HLL vs HLLC at [{}]: {} vs {}",
                i,
                f_hll[i],
                f_hllc[i]
            );
        }
    }

    // ── Characteristic decomposition tests ───────────────────────────────────

    #[test]
    fn test_characteristic_roundtrip() {
        let u0 = 1.5;
        let u1 = 0.7;
        let a = 2.0;
        let (wp, wm) = characteristic_decompose_2wave(u0, u1, a);
        let (u0_back, u1_back) = characteristic_recompose_2wave(wp, wm, a);
        assert!((u0_back - u0).abs() < 1e-12, "roundtrip u0: {u0_back}");
        assert!((u1_back - u1).abs() < 1e-12, "roundtrip u1: {u1_back}");
    }

    #[test]
    fn test_characteristic_limited_reconstruct_finite() {
        let u0: Vec<f64> = (0..8).map(|i| i as f64 * 0.1).collect();
        let u1: Vec<f64> = vec![0.5; 8];
        let (left, right) = characteristic_limited_reconstruct(&u0, &u1, 3, 1.0, minmod);
        for &v in left.iter().chain(right.iter()) {
            assert!(v.is_finite(), "characteristic reconstruct not finite: {v}");
        }
    }

    // ── TVD RK2 tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_tvd_rk2_advect_uniform_unchanged() {
        let phi = FluxGrid3D::new(5, 5, 5, 0.1, 2.0);
        let mut vel = VectorFluxGrid3D::zeros(5, 5, 5, 0.1);
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    vel.set_vel(i, j, k, [1.0, 0.0, 0.0]);
                }
            }
        }
        let phi_new = tvd_rk2_advect(&phi, &vel, 0.01);
        for (&a, &b) in phi.values.iter().zip(phi_new.values.iter()) {
            assert!((a - b).abs() < 1e-10, "TVD RK2 uniform: {a} vs {b}");
        }
    }

    #[test]
    fn test_tvd_rk2_advect_size_preserved() {
        let phi = FluxGrid3D::new(4, 3, 2, 0.1, 1.0);
        let vel = VectorFluxGrid3D::zeros(4, 3, 2, 0.1);
        let phi_new = tvd_rk2_advect(&phi, &vel, 0.01);
        assert_eq!(phi_new.nx, phi.nx);
        assert_eq!(phi_new.ny, phi.ny);
        assert_eq!(phi_new.nz, phi.nz);
    }

    // ── Limiter ratio-form tests ──────────────────────────────────────────────

    #[test]
    fn test_minmod_ratio_clamp() {
        assert_eq!(minmod_ratio(-0.5), 0.0);
        assert!((minmod_ratio(0.5) - 0.5).abs() < 1e-12);
        assert!((minmod_ratio(2.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_superbee_ratio_basic() {
        assert_eq!(superbee_ratio(-1.0), 0.0);
        let v = superbee_ratio(1.5);
        assert!(v > 0.0 && v <= 2.0, "superbee_ratio(1.5) = {v}");
    }

    #[test]
    fn test_van_leer_ratio_basic() {
        // van_leer_ratio(1.0) = 2*1/(1+1) = 1.0
        let v = van_leer_ratio(1.0);
        assert!((v - 1.0).abs() < 1e-12, "van_leer_ratio(1) = {v}");
        // Negative r → 0
        assert_eq!(van_leer_ratio(-0.5), 0.0);
    }

    #[test]
    fn test_van_leer_ratio_monotone() {
        // phi(r) should be monotone increasing for r > 0
        let v1 = van_leer_ratio(0.5);
        let v2 = van_leer_ratio(1.0);
        let v3 = van_leer_ratio(2.0);
        assert!(v1 <= v2, "van_leer not monotone: {v1} > {v2}");
        assert!(v2 <= v3 + 1e-12, "van_leer not monotone: {v2} > {v3}");
    }
}
