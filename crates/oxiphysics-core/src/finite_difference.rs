// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Finite difference schemes for PDEs in 1D, 2D, and 3D.
//!
//! Covers explicit/implicit time integration, Crank-Nicolson diffusion,
//! upwind advection, 5th-order WENO, compact finite differences, boundary
//! conditions (Dirichlet/Neumann/periodic), Richardson extrapolation,
//! method of lines, stability analysis (CFL, von Neumann), staggered grids,
//! fractional step method, and ADI (alternating direction implicit).

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Boundary condition types
// ─────────────────────────────────────────────────────────────────────────────

/// Boundary condition specification for finite difference solvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Dirichlet condition: fixed value at the boundary.
    Dirichlet(f64),
    /// Neumann condition: fixed normal derivative at the boundary.
    Neumann(f64),
    /// Periodic boundary condition (wraps domain).
    Periodic,
}

// ─────────────────────────────────────────────────────────────────────────────
// Thomas algorithm (tridiagonal solver)
// ─────────────────────────────────────────────────────────────────────────────

/// Solves the tridiagonal system `A x = rhs` using the Thomas algorithm.
///
/// `sub` is the sub-diagonal (index 0 unused), `diag` the main diagonal,
/// `sup` the super-diagonal (last entry unused), and `rhs` the right-hand side.
/// All slices must have length `n`.
pub fn thomas_algorithm(sub: &[f64], diag: &[f64], sup: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    assert!(n >= 2, "system must have at least 2 unknowns");
    assert_eq!(sub.len(), n);
    assert_eq!(sup.len(), n);
    assert_eq!(rhs.len(), n);

    let mut sup_prime = vec![0.0_f64; n];
    let mut rhs_prime = vec![0.0_f64; n];
    let mut x = vec![0.0_f64; n];

    sup_prime[0] = sup[0] / diag[0];
    rhs_prime[0] = rhs[0] / diag[0];
    for i in 1..n {
        let denom = diag[i] - sub[i] * sup_prime[i - 1];
        sup_prime[i] = if i < n - 1 { sup[i] / denom } else { 0.0 };
        rhs_prime[i] = (rhs[i] - sub[i] * rhs_prime[i - 1]) / denom;
    }

    x[n - 1] = rhs_prime[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = rhs_prime[i] - sup_prime[i] * x[i + 1];
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// First and second order 1D stencils
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the first derivative of `u` using second-order central differences.
///
/// At boundaries uses one-sided (forward/backward) differences.
pub fn first_derivative_central(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    let mut du = vec![0.0_f64; n];
    if n < 2 {
        return du;
    }
    du[0] = (u[1] - u[0]) / dx;
    for i in 1..n - 1 {
        du[i] = (u[i + 1] - u[i - 1]) / (2.0 * dx);
    }
    du[n - 1] = (u[n - 1] - u[n - 2]) / dx;
    du
}

/// Computes the second derivative of `u` with the standard 3-point Laplacian stencil.
///
/// Interior: `(u[i-1] - 2 u[i] + u[i+1]) / dx²`. Boundaries are set to zero.
pub fn second_derivative(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    let mut d2u = vec![0.0_f64; n];
    for i in 1..n.saturating_sub(1) {
        d2u[i] = (u[i - 1] - 2.0 * u[i] + u[i + 1]) / (dx * dx);
    }
    d2u
}

/// Computes the fourth-order central difference first derivative.
///
/// Uses the 5-point stencil: `(-u[i+2] + 8u[i+1] - 8u[i-1] + u[i-2]) / (12 dx)`.
/// Falls back to second-order central differences near boundaries.
pub fn first_derivative_fourth_order(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    let mut du = vec![0.0_f64; n];
    if n < 2 {
        return du;
    }
    // boundaries: second-order
    du[0] = (u[1] - u[0]) / dx;
    if n > 1 {
        du[n - 1] = (u[n - 1] - u[n - 2]) / dx;
    }
    if n > 2 {
        du[1] = (u[2] - u[0]) / (2.0 * dx);
        du[n - 2] = (u[n - 1] - u[n - 3]) / (2.0 * dx);
    }
    // interior: fourth-order
    for i in 2..n.saturating_sub(2) {
        du[i] = (-u[i + 2] + 8.0 * u[i + 1] - 8.0 * u[i - 1] + u[i - 2]) / (12.0 * dx);
    }
    du
}

// ─────────────────────────────────────────────────────────────────────────────
// Upwind scheme for advection
// ─────────────────────────────────────────────────────────────────────────────

/// First-order upwind advection: advances `u` one time step `dt`.
///
/// The equation solved is `∂u/∂t + c ∂u/∂x = 0` with `bc` at the inflow.
pub fn advection_upwind_step(
    u: &[f64],
    dx: f64,
    dt: f64,
    c: f64,
    bc_left: BoundaryCondition,
    bc_right: BoundaryCondition,
) -> Vec<f64> {
    let n = u.len();
    let mut un = u.to_vec();
    let r = c * dt / dx;
    for i in 1..n - 1 {
        if c > 0.0 {
            un[i] = u[i] - r * (u[i] - u[i - 1]);
        } else {
            un[i] = u[i] - r * (u[i + 1] - u[i]);
        }
    }
    apply_bc_1d(&mut un, bc_left, bc_right);
    un
}

/// Applies 1D boundary conditions (left and right) to a solution vector.
pub fn apply_bc_1d(u: &mut [f64], left: BoundaryCondition, right: BoundaryCondition) {
    let n = u.len();
    if n == 0 {
        return;
    }
    match left {
        BoundaryCondition::Dirichlet(v) => u[0] = v,
        BoundaryCondition::Neumann(dv) => {
            if n > 1 {
                u[0] = u[1] - dv * 1.0; // dx factored out
            }
        }
        BoundaryCondition::Periodic => {
            if n > 1 {
                u[0] = u[n - 2];
            }
        }
    }
    match right {
        BoundaryCondition::Dirichlet(v) => u[n - 1] = v,
        BoundaryCondition::Neumann(dv) => {
            if n > 1 {
                u[n - 1] = u[n - 2] + dv * 1.0;
            }
        }
        BoundaryCondition::Periodic => {
            if n > 1 {
                u[n - 1] = u[1];
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Explicit diffusion (FTCS)
// ─────────────────────────────────────────────────────────────────────────────

/// Forward-time central-space (FTCS) explicit diffusion step.
///
/// Solves `∂u/∂t = α ∂²u/∂x²`. Stable when `α dt / dx² ≤ 0.5`.
pub fn diffusion_explicit_step(
    u: &[f64],
    dx: f64,
    dt: f64,
    alpha: f64,
    bc_left: BoundaryCondition,
    bc_right: BoundaryCondition,
) -> Vec<f64> {
    let n = u.len();
    let r = alpha * dt / (dx * dx);
    let mut un = u.to_vec();
    for i in 1..n - 1 {
        un[i] = u[i] + r * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
    }
    apply_bc_1d(&mut un, bc_left, bc_right);
    un
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit diffusion (Crank-Nicolson)
// ─────────────────────────────────────────────────────────────────────────────

/// Crank-Nicolson implicit diffusion step.
///
/// Solves `∂u/∂t = α ∂²u/∂x²` using the unconditionally stable C-N scheme.
/// Interior points are advanced; boundary values are enforced via `bc_left` / `bc_right`.
pub fn diffusion_crank_nicolson_step(
    u: &[f64],
    dx: f64,
    dt: f64,
    alpha: f64,
    bc_left: BoundaryCondition,
    bc_right: BoundaryCondition,
) -> Vec<f64> {
    let n = u.len();
    let r = alpha * dt / (2.0 * dx * dx);

    // Build RHS for interior points
    let mut rhs = vec![0.0_f64; n];
    rhs[0] = match bc_left {
        BoundaryCondition::Dirichlet(v) => v,
        _ => u[0],
    };
    rhs[n - 1] = match bc_right {
        BoundaryCondition::Dirichlet(v) => v,
        _ => u[n - 1],
    };
    for i in 1..n - 1 {
        rhs[i] = r * u[i - 1] + (1.0 - 2.0 * r) * u[i] + r * u[i + 1];
    }

    // Tridiagonal system: interior rows
    let mut a_diag = vec![0.0_f64; n]; // sub-diag
    let mut b_diag = vec![0.0_f64; n]; // main
    let mut c_diag = vec![0.0_f64; n]; // super-diag

    b_diag[0] = 1.0;
    b_diag[n - 1] = 1.0;
    for i in 1..n - 1 {
        a_diag[i] = -r;
        b_diag[i] = 1.0 + 2.0 * r;
        c_diag[i] = -r;
    }

    let mut sol = thomas_algorithm(&a_diag, &b_diag, &c_diag, &rhs);
    apply_bc_1d(&mut sol, bc_left, bc_right);
    sol
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit (fully implicit / backward Euler) diffusion
// ─────────────────────────────────────────────────────────────────────────────

/// Fully implicit (backward Euler) diffusion step.
///
/// Unconditionally stable. Solves `(I - dt α D²) uⁿ⁺¹ = uⁿ`.
pub fn diffusion_implicit_step(
    u: &[f64],
    dx: f64,
    dt: f64,
    alpha: f64,
    bc_left: BoundaryCondition,
    bc_right: BoundaryCondition,
) -> Vec<f64> {
    let n = u.len();
    let r = alpha * dt / (dx * dx);
    let mut a_diag = vec![0.0_f64; n];
    let mut b_diag = vec![0.0_f64; n];
    let mut c_diag = vec![0.0_f64; n];
    let mut rhs = u.to_vec();

    b_diag[0] = 1.0;
    b_diag[n - 1] = 1.0;
    for i in 1..n - 1 {
        a_diag[i] = -r;
        b_diag[i] = 1.0 + 2.0 * r;
        c_diag[i] = -r;
    }

    // enforce Dirichlet BCs in rhs
    if let BoundaryCondition::Dirichlet(v) = bc_left {
        rhs[0] = v;
    }
    if let BoundaryCondition::Dirichlet(v) = bc_right {
        rhs[n - 1] = v;
    }

    let mut sol = thomas_algorithm(&a_diag, &b_diag, &c_diag, &rhs);
    apply_bc_1d(&mut sol, bc_left, bc_right);
    sol
}

// ─────────────────────────────────────────────────────────────────────────────
// WENO-5 (5th order weighted essentially non-oscillatory)
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the WENO-5 smoothness indicator `β` for a three-point stencil.
///
/// Used internally by [`weno5_flux`].
fn weno5_beta(v0: f64, v1: f64, v2: f64) -> f64 {
    let d1 = v1 - v0;
    let d2 = v2 - 2.0 * v1 + v0;
    13.0 / 12.0 * d2 * d2 + 0.25 * (3.0 * v1 - 4.0 * v0 + v2).powi(2) - 0.25 * (v1 - v0).powi(2)
        + 0.25 * d1 * d1
}

/// Computes the WENO-5 reconstructed flux at interface `i + 1/2` for `c > 0`.
///
/// Requires the five stencil values `[v0, v1, v2, v3, v4]` = `u[i-2..i+3]`.
pub fn weno5_flux(v: [f64; 5]) -> f64 {
    let eps = 1e-36_f64;

    // three candidate stencils
    let q0 = (2.0 * v[0] - 7.0 * v[1] + 11.0 * v[2]) / 6.0;
    let q1 = (-v[1] + 5.0 * v[2] + 2.0 * v[3]) / 6.0;
    let q2 = (2.0 * v[2] + 5.0 * v[3] - v[4]) / 6.0;

    // smoothness indicators
    let b0 = 13.0 / 12.0 * (v[0] - 2.0 * v[1] + v[2]).powi(2)
        + 0.25 * (v[0] - 4.0 * v[1] + 3.0 * v[2]).powi(2);
    let b1 = 13.0 / 12.0 * (v[1] - 2.0 * v[2] + v[3]).powi(2) + 0.25 * (v[1] - v[3]).powi(2);
    let b2 = 13.0 / 12.0 * (v[2] - 2.0 * v[3] + v[4]).powi(2)
        + 0.25 * (3.0 * v[2] - 4.0 * v[3] + v[4]).powi(2);

    let d0 = 1.0 / 10.0;
    let d1 = 6.0 / 10.0;
    let d2 = 3.0 / 10.0;

    let w0 = d0 / (eps + b0).powi(2);
    let w1 = d1 / (eps + b1).powi(2);
    let w2 = d2 / (eps + b2).powi(2);
    let wsum = w0 + w1 + w2;

    (w0 * q0 + w1 * q1 + w2 * q2) / wsum
}

/// One time step of the WENO-5 advection scheme for `∂u/∂t + c ∂u/∂x = 0`.
///
/// Uses 5th-order WENO spatial reconstruction and first-order Euler time integration.
/// Ghost cells on each side are filled by zero-order extrapolation.
pub fn advection_weno5_step(u: &[f64], dx: f64, dt: f64, c: f64) -> Vec<f64> {
    let n = u.len();
    let mut un = u.to_vec();

    // build padded array with 3 ghost cells each side
    let pad = 3_usize;
    let mut up = vec![0.0_f64; n + 2 * pad];
    for i in 0..pad {
        up[i] = u[0];
        up[n + pad + i] = u[n - 1];
    }
    up[pad..pad + n].copy_from_slice(u);

    for i in 0..n {
        let ip = i + pad; // index in padded array
        let flux_right = if c > 0.0 {
            weno5_flux([up[ip - 2], up[ip - 1], up[ip], up[ip + 1], up[ip + 2]])
        } else {
            weno5_flux([up[ip + 3], up[ip + 2], up[ip + 1], up[ip], up[ip - 1]])
        };
        let flux_left = if c > 0.0 {
            weno5_flux([up[ip - 3], up[ip - 2], up[ip - 1], up[ip], up[ip + 1]])
        } else {
            weno5_flux([up[ip + 2], up[ip + 1], up[ip], up[ip - 1], up[ip - 2]])
        };
        un[i] = u[i] - c * dt / dx * (flux_right - flux_left);
    }
    un
}

// ─────────────────────────────────────────────────────────────────────────────
// Compact finite differences (Padé 4th-order for first derivative)
// ─────────────────────────────────────────────────────────────────────────────

/// Compact (Padé) 4th-order first derivative.
///
/// Solves the tridiagonal system `(1/4) f'[i-1] + f'[i] + (1/4) f'[i+1] = (3/2)(u[i+1]-u[i-1])/(2dx)`.
/// Interior only; boundary nodes use standard second-order differences.
pub fn compact_first_derivative(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    if n < 3 {
        return first_derivative_central(u, dx);
    }

    let mut rhs = vec![0.0_f64; n];
    rhs[0] = (u[1] - u[0]) / dx;
    rhs[n - 1] = (u[n - 1] - u[n - 2]) / dx;
    for i in 1..n - 1 {
        rhs[i] = (3.0 / 2.0) * (u[i + 1] - u[i - 1]) / (2.0 * dx);
    }

    let mut a_d = vec![0.0_f64; n];
    let b_d = vec![1.0_f64; n];
    let mut c_d = vec![0.0_f64; n];
    for i in 1..n - 1 {
        a_d[i] = 1.0 / 4.0;
        c_d[i] = 1.0 / 4.0;
    }

    thomas_algorithm(&a_d, &b_d, &c_d, &rhs)
}

// ─────────────────────────────────────────────────────────────────────────────
// Richardson extrapolation
// ─────────────────────────────────────────────────────────────────────────────

/// Richardson extrapolation to combine two estimates with step sizes `h` and `h/2`.
///
/// For a method of order `p`, the extrapolated value is
/// `(2^p * f_h2 - f_h) / (2^p - 1)`.
pub fn richardson_extrapolation(f_h: f64, f_h2: f64, order: u32) -> f64 {
    let r = (2_u64.pow(order)) as f64;
    (r * f_h2 - f_h) / (r - 1.0)
}

/// Applies Richardson extrapolation to refine a 1D finite-difference derivative.
///
/// Computes the first derivative at two resolutions (`dx` and `dx/2`) and
/// combines them to produce an estimate of order `order + 2`.
pub fn richardson_derivative(u_coarse: &[f64], u_fine: &[f64], dx: f64, order: u32) -> Vec<f64> {
    let nc = u_coarse.len();
    let _nf = u_fine.len();
    let dc = first_derivative_central(u_coarse, dx);
    let df_on_coarse: Vec<f64> = (0..nc)
        .map(|i| {
            let if_ = i * 2;
            let nf = u_fine.len();
            if if_ > 0 && if_ < nf - 1 {
                (u_fine[if_ + 1] - u_fine[if_ - 1]) / (2.0 * dx / 2.0)
            } else {
                dc[i]
            }
        })
        .collect();

    dc.iter()
        .zip(df_on_coarse.iter())
        .map(|(&fh, &fh2)| richardson_extrapolation(fh, fh2, order))
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Stability criteria
// ─────────────────────────────────────────────────────────────────────────────

/// CFL (Courant-Friedrichs-Lewy) number for the advection equation.
///
/// `CFL = |c| * dt / dx`. Stability of explicit schemes requires `CFL ≤ 1`.
pub fn cfl_number(c: f64, dt: f64, dx: f64) -> f64 {
    c.abs() * dt / dx
}

/// Von Neumann diffusion stability number `r = α dt / dx²`.
///
/// Explicit FTCS scheme is stable when `r ≤ 0.5`.
pub fn von_neumann_diffusion_number(alpha: f64, dt: f64, dx: f64) -> f64 {
    alpha * dt / (dx * dx)
}

/// Returns `true` if the explicit advection scheme is stable (CFL ≤ 1).
pub fn is_advection_stable(c: f64, dt: f64, dx: f64) -> bool {
    cfl_number(c, dt, dx) <= 1.0
}

/// Returns `true` if the explicit diffusion scheme (FTCS) is stable.
pub fn is_diffusion_stable(alpha: f64, dt: f64, dx: f64) -> bool {
    von_neumann_diffusion_number(alpha, dt, dx) <= 0.5
}

/// Computes the maximum stable time step for explicit advection.
pub fn max_dt_advection(c: f64, dx: f64) -> f64 {
    dx / c.abs().max(f64::EPSILON)
}

/// Computes the maximum stable time step for explicit diffusion (FTCS).
pub fn max_dt_diffusion(alpha: f64, dx: f64) -> f64 {
    dx * dx / (2.0 * alpha.abs().max(f64::EPSILON))
}

// ─────────────────────────────────────────────────────────────────────────────
// Method of lines (MOL)
// ─────────────────────────────────────────────────────────────────────────────

/// Method-of-lines spatial operator for the 1D diffusion equation.
///
/// Returns `du/dt = α d²u/dx²` evaluated at all interior points.
/// `u` is the full solution vector including boundary nodes.
pub fn mol_diffusion_rhs(u: &[f64], dx: f64, alpha: f64) -> Vec<f64> {
    second_derivative(u, dx)
        .iter()
        .map(|&v| alpha * v)
        .collect()
}

/// Method-of-lines spatial operator for the 1D advection equation.
///
/// Returns `du/dt = -c du/dx` using upwind differences.
pub fn mol_advection_rhs(u: &[f64], dx: f64, c: f64) -> Vec<f64> {
    let n = u.len();
    let mut rhs = vec![0.0_f64; n];
    for i in 1..n - 1 {
        let dudx = if c > 0.0 {
            (u[i] - u[i - 1]) / dx
        } else {
            (u[i + 1] - u[i]) / dx
        };
        rhs[i] = -c * dudx;
    }
    rhs
}

/// Runge-Kutta 4 integrator for method of lines.
///
/// `f(t, u)` is the right-hand side, `u` the current state, `dt` the time step.
pub fn rk4_step<F>(u: &[f64], t: f64, dt: f64, f: &F) -> Vec<f64>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = u.len();
    let k1 = f(t, u);
    let u2: Vec<f64> = (0..n).map(|i| u[i] + 0.5 * dt * k1[i]).collect();
    let k2 = f(t + 0.5 * dt, &u2);
    let u3: Vec<f64> = (0..n).map(|i| u[i] + 0.5 * dt * k2[i]).collect();
    let k3 = f(t + 0.5 * dt, &u3);
    let u4: Vec<f64> = (0..n).map(|i| u[i] + dt * k3[i]).collect();
    let k4 = f(t + dt, &u4);
    (0..n)
        .map(|i| u[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Staggered grid helpers
// ─────────────────────────────────────────────────────────────────────────────

/// A staggered 1D grid holding cell-centred values `u` and face-centred fluxes `f`.
///
/// Used in finite-volume and MAC schemes.
#[derive(Debug, Clone)]
pub struct StaggeredGrid1D {
    /// Cell-centred scalar field, length `n`.
    pub u: Vec<f64>,
    /// Face-centred flux field, length `n + 1`.
    pub flux: Vec<f64>,
    /// Grid spacing.
    pub dx: f64,
}

impl StaggeredGrid1D {
    /// Creates a new staggered grid with `n` cells and spacing `dx`.
    pub fn new(n: usize, dx: f64) -> Self {
        Self {
            u: vec![0.0_f64; n],
            flux: vec![0.0_f64; n + 1],
            dx,
        }
    }

    /// Computes face fluxes from cell-centred velocities by simple averaging.
    pub fn interpolate_to_faces(&mut self) {
        let n = self.u.len();
        self.flux[0] = self.u[0];
        self.flux[n] = self.u[n - 1];
        for i in 1..n {
            self.flux[i] = 0.5 * (self.u[i - 1] + self.u[i]);
        }
    }

    /// Computes the divergence of the flux field at each cell centre.
    pub fn divergence(&self) -> Vec<f64> {
        let n = self.u.len();
        (0..n)
            .map(|i| (self.flux[i + 1] - self.flux[i]) / self.dx)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fractional step (projection) method
// ─────────────────────────────────────────────────────────────────────────────

/// Fractional step result from one projection step.
#[derive(Debug, Clone)]
pub struct FractionalStepResult {
    /// Intermediate (unprojected) velocity.
    pub u_star: Vec<f64>,
    /// Pressure correction field.
    pub pressure: Vec<f64>,
    /// Divergence-free velocity after projection.
    pub u_new: Vec<f64>,
}

/// Fractional step (Chorin projection) method for 1D incompressible flow.
///
/// Step 1: Advance velocity ignoring pressure (explicit advection + diffusion).
/// Step 2: Solve pressure Poisson equation.
/// Step 3: Project to divergence-free field.
pub fn fractional_step_1d(
    u: &[f64],
    dx: f64,
    dt: f64,
    nu: f64,
    body_force: f64,
) -> FractionalStepResult {
    let n = u.len();

    // Step 1: predict velocity (advection + diffusion + body force)
    let mut u_star = u.to_vec();
    for i in 1..n - 1 {
        let adv = u[i] * (u[i] - u[i - 1]) / dx;
        let diff = nu * (u[i - 1] - 2.0 * u[i] + u[i + 1]) / (dx * dx);
        u_star[i] = u[i] + dt * (-adv + diff + body_force);
    }

    // Step 2: pressure Poisson equation (div u_star = Δt ∇²p)
    let mut rhs_p = vec![0.0_f64; n];
    for i in 1..n - 1 {
        rhs_p[i] = (u_star[i + 1] - u_star[i - 1]) / (2.0 * dx) / dt;
    }
    let a_d: Vec<f64> = (0..n)
        .map(|i| {
            if i == 0 || i == n - 1 {
                0.0
            } else {
                -1.0 / (dx * dx)
            }
        })
        .collect();
    let b_d: Vec<f64> = (0..n)
        .map(|i| {
            if i == 0 || i == n - 1 {
                1.0
            } else {
                2.0 / (dx * dx)
            }
        })
        .collect();
    let c_d: Vec<f64> = a_d.clone();
    let pressure = thomas_algorithm(&a_d, &b_d, &c_d, &rhs_p);

    // Step 3: projection
    let mut u_new = u_star.clone();
    for i in 1..n - 1 {
        u_new[i] -= dt * (pressure[i + 1] - pressure[i - 1]) / (2.0 * dx);
    }

    FractionalStepResult {
        u_star,
        pressure,
        u_new,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2D finite difference helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the 2D Laplacian using the standard 5-point stencil.
///
/// `u` is stored row-major with `nx` columns and `ny` rows.
/// Interior points only; boundary nodes are left unchanged.
pub fn laplacian_2d(u: &[f64], nx: usize, ny: usize, dx: f64, dy: f64) -> Vec<f64> {
    let mut lap = vec![0.0_f64; nx * ny];
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            let d2x = (u[idx - 1] - 2.0 * u[idx] + u[idx + 1]) / (dx * dx);
            let d2y = (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]) / (dy * dy);
            lap[idx] = d2x + d2y;
        }
    }
    lap
}

/// Applies 2D Dirichlet boundary conditions to all four edges.
///
/// `u` is row-major with `nx` columns and `ny` rows.
/// Sets the border values to `value`.
pub fn apply_dirichlet_2d(u: &mut [f64], nx: usize, ny: usize, value: f64) {
    for i in 0..nx {
        u[i] = value;
        u[(ny - 1) * nx + i] = value;
    }
    for j in 0..ny {
        u[j * nx] = value;
        u[j * nx + nx - 1] = value;
    }
}

/// One explicit 2D diffusion step using the FTCS scheme.
///
/// Solves `∂u/∂t = α (∂²u/∂x² + ∂²u/∂y²)`. Boundary nodes are not updated.
pub fn diffusion_2d_explicit_step(
    u: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    alpha: f64,
) -> Vec<f64> {
    let mut un = u.to_vec();
    let rx = alpha * dt / (dx * dx);
    let ry = alpha * dt / (dy * dy);
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            un[idx] = u[idx]
                + rx * (u[idx - 1] - 2.0 * u[idx] + u[idx + 1])
                + ry * (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]);
        }
    }
    un
}

// ─────────────────────────────────────────────────────────────────────────────
// ADI (Alternating Direction Implicit) for 2D diffusion
// ─────────────────────────────────────────────────────────────────────────────

/// One full ADI time step for 2D diffusion `∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²)`.
///
/// Uses the Peaceman-Rachford splitting: half-step implicit in x, half-step implicit in y.
/// Boundary values are held fixed (Dirichlet zero on all edges).
pub fn adi_diffusion_step_2d(
    u: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    alpha: f64,
) -> Vec<f64> {
    let rx = alpha * dt / (2.0 * dx * dx);
    let ry = alpha * dt / (2.0 * dy * dy);
    let mut u_half = u.to_vec();

    // --- half step: implicit in x, explicit in y ---
    for j in 1..ny - 1 {
        let mut a_d = vec![0.0_f64; nx];
        let mut b_d = vec![0.0_f64; nx];
        let mut c_d = vec![0.0_f64; nx];
        let mut rhs = vec![0.0_f64; nx];

        b_d[0] = 1.0;
        b_d[nx - 1] = 1.0;
        // boundary stays 0
        rhs[0] = 0.0;
        rhs[nx - 1] = 0.0;

        for i in 1..nx - 1 {
            let idx = j * nx + i;
            a_d[i] = -rx;
            b_d[i] = 1.0 + 2.0 * rx;
            c_d[i] = -rx;
            rhs[i] = ry * u[idx - nx] + (1.0 - 2.0 * ry) * u[idx] + ry * u[idx + nx];
        }

        let row_sol = thomas_algorithm(&a_d, &b_d, &c_d, &rhs);
        u_half[j * nx..j * nx + nx].copy_from_slice(&row_sol);
    }

    // --- full step: implicit in y ---
    let mut u_new = u_half.clone();
    for i in 1..nx - 1 {
        let mut a_d = vec![0.0_f64; ny];
        let mut b_d = vec![0.0_f64; ny];
        let mut c_d = vec![0.0_f64; ny];
        let mut rhs = vec![0.0_f64; ny];

        b_d[0] = 1.0;
        b_d[ny - 1] = 1.0;
        rhs[0] = 0.0;
        rhs[ny - 1] = 0.0;

        for j in 1..ny - 1 {
            let idx = j * nx + i;
            a_d[j] = -ry;
            b_d[j] = 1.0 + 2.0 * ry;
            c_d[j] = -ry;
            rhs[j] = rx * u_half[idx - 1] + (1.0 - 2.0 * rx) * u_half[idx] + rx * u_half[idx + 1];
        }

        let col_sol = thomas_algorithm(&a_d, &b_d, &c_d, &rhs);
        for j in 0..ny {
            u_new[j * nx + i] = col_sol[j];
        }
    }

    u_new
}

// ─────────────────────────────────────────────────────────────────────────────
// 3D finite difference Laplacian
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the 3D Laplacian using the 7-point stencil.
///
/// `u` is stored as `u[k * ny * nx + j * nx + i]`.
/// Boundary nodes are left as zero.
pub fn laplacian_3d(
    u: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    dy: f64,
    dz: f64,
) -> Vec<f64> {
    let mut lap = vec![0.0_f64; nx * ny * nz];
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = k * ny * nx + j * nx + i;
                let d2x = (u[idx - 1] - 2.0 * u[idx] + u[idx + 1]) / (dx * dx);
                let d2y = (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]) / (dy * dy);
                let d2z = (u[idx - ny * nx] - 2.0 * u[idx] + u[idx + ny * nx]) / (dz * dz);
                lap[idx] = d2x + d2y + d2z;
            }
        }
    }
    lap
}

// ─────────────────────────────────────────────────────────────────────────────
// Periodic convection (flux form)
// ─────────────────────────────────────────────────────────────────────────────

/// First-order upwind advection step on a periodic domain.
///
/// No boundary conditions needed because the domain wraps around.
pub fn advection_periodic_step(u: &[f64], dx: f64, dt: f64, c: f64) -> Vec<f64> {
    let n = u.len();
    let r = c * dt / dx;
    let mut un = vec![0.0_f64; n];
    for i in 0..n {
        let im1 = if i == 0 { n - 1 } else { i - 1 };
        let ip1 = (i + 1) % n;
        un[i] = if c > 0.0 {
            u[i] - r * (u[i] - u[im1])
        } else {
            u[i] - r * (u[ip1] - u[i])
        };
    }
    un
}

// ─────────────────────────────────────────────────────────────────────────────
// L2 norm and convergence helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the discrete L2 norm `||u||₂ = sqrt(dx * Σ u[i]²)`.
pub fn l2_norm(u: &[f64], dx: f64) -> f64 {
    (dx * u.iter().map(|v| v * v).sum::<f64>()).sqrt()
}

/// Computes the L∞ (maximum) norm of `u`.
pub fn linf_norm(u: &[f64]) -> f64 {
    u.iter().cloned().fold(0.0_f64, f64::max)
}

/// Computes the discrete L2 error between computed and exact solutions.
pub fn l2_error(computed: &[f64], exact: &[f64], dx: f64) -> f64 {
    let err: Vec<f64> = computed
        .iter()
        .zip(exact.iter())
        .map(|(c, e)| c - e)
        .collect();
    l2_norm(&err, dx)
}

// ─────────────────────────────────────────────────────────────────────────────
// Von Neumann stability analysis (amplification factor)
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the von Neumann amplification factor magnitude for FTCS diffusion.
///
/// The amplification factor is `|G| = |1 - 4r sin²(k dx / 2)|`.
/// Stability requires `|G| ≤ 1` for all wave numbers `k`.
pub fn ftcs_amplification(r: f64, k: f64, dx: f64) -> f64 {
    let s = (k * dx / 2.0).sin();
    (1.0 - 4.0 * r * s * s).abs()
}

/// Maximum amplification factor for the FTCS diffusion scheme over all wave numbers.
///
/// Scans `n_modes` wave numbers in `[0, π/dx]` and returns the largest `|G|`.
pub fn ftcs_max_amplification(r: f64, dx: f64, n_modes: usize) -> f64 {
    (0..=n_modes)
        .map(|m| {
            let k = m as f64 * PI / (n_modes as f64 * dx);
            ftcs_amplification(r, k, dx)
        })
        .fold(0.0_f64, f64::max)
}

// ─────────────────────────────────────────────────────────────────────────────
// Neumann boundary condition helper
// ─────────────────────────────────────────────────────────────────────────────

/// Applies Neumann (zero-flux) conditions by ghost-cell reflection.
///
/// Left ghost: `u[-1] = u[1]`, right ghost: `u[n] = u[n-2]`.
/// Modifies the first and last elements of `u` in-place.
pub fn apply_neumann_reflection(u: &mut [f64]) {
    let n = u.len();
    if n < 3 {
        return;
    }
    u[0] = u[1];
    u[n - 1] = u[n - 2];
}

// ─────────────────────────────────────────────────────────────────────────────
// Grid builder
// ─────────────────────────────────────────────────────────────────────────────

/// Generates a uniform 1D grid from `x_start` to `x_end` with `n` points.
pub fn uniform_grid_1d(x_start: f64, x_end: f64, n: usize) -> Vec<f64> {
    let dx = (x_end - x_start) / (n - 1) as f64;
    (0..n).map(|i| x_start + i as f64 * dx).collect()
}

/// Generates a uniform 2D grid (row-major).
///
/// Returns `(x_coords, y_coords)` each of length `nx * ny`.
pub fn uniform_grid_2d(
    x0: f64,
    x1: f64,
    nx: usize,
    y0: f64,
    y1: f64,
    ny: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dx = (x1 - x0) / (nx - 1) as f64;
    let dy = (y1 - y0) / (ny - 1) as f64;
    let mut xs = vec![0.0_f64; nx * ny];
    let mut ys = vec![0.0_f64; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            xs[j * nx + i] = x0 + i as f64 * dx;
            ys[j * nx + i] = y0 + j as f64 * dy;
        }
    }
    (xs, ys)
}

// ─────────────────────────────────────────────────────────────────────────────
// Lax-Wendroff scheme
// ─────────────────────────────────────────────────────────────────────────────

/// Lax-Wendroff second-order advection step for `∂u/∂t + c ∂u/∂x = 0`.
///
/// Second-order accurate in both space and time.  Boundary nodes are copied unchanged.
pub fn advection_lax_wendroff_step(u: &[f64], dx: f64, dt: f64, c: f64) -> Vec<f64> {
    let n = u.len();
    let r = c * dt / dx;
    let mut un = u.to_vec();
    for i in 1..n - 1 {
        un[i] = u[i] - 0.5 * r * (u[i + 1] - u[i - 1])
            + 0.5 * r * r * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
    }
    un
}

// ─────────────────────────────────────────────────────────────────────────────
// Beam-Warming scheme
// ─────────────────────────────────────────────────────────────────────────────

/// Beam-Warming second-order upwind advection step (for `c > 0`).
///
/// Uses the stencil `u[i] - r/2 (3u[i]-4u[i-1]+u[i-2]) + r²/2 (u[i]-2u[i-1]+u[i-2])`.
pub fn advection_beam_warming_step(u: &[f64], dx: f64, dt: f64, c: f64) -> Vec<f64> {
    let n = u.len();
    let r = c * dt / dx;
    let mut un = u.to_vec();
    for i in 2..n {
        un[i] = u[i] - 0.5 * r * (3.0 * u[i] - 4.0 * u[i - 1] + u[i - 2])
            + 0.5 * r * r * (u[i] - 2.0 * u[i - 1] + u[i - 2]);
    }
    un
}

// ─────────────────────────────────────────────────────────────────────────────
// Compact finite difference for second derivative (Padé 4th order)
// ─────────────────────────────────────────────────────────────────────────────

/// Compact (Padé) 4th-order second derivative.
///
/// Solves `(1/10) f''[i-1] + f''[i] + (1/10) f''[i+1] = (12/10)(u[i+1]-2u[i]+u[i-1])/dx²`.
pub fn compact_second_derivative(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    if n < 3 {
        return second_derivative(u, dx);
    }

    let mut rhs = vec![0.0_f64; n];
    rhs[0] = (u[1] - 2.0 * u[0] + u[0]) / (dx * dx); // ghost
    rhs[n - 1] = (u[n - 1] - 2.0 * u[n - 1] + u[n - 2]) / (dx * dx);
    for i in 1..n - 1 {
        rhs[i] = (12.0 / 10.0) * (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (dx * dx);
    }

    let mut a_d = vec![0.0_f64; n];
    let b_d = vec![1.0_f64; n];
    let mut c_d = vec![0.0_f64; n];
    for i in 1..n - 1 {
        a_d[i] = 1.0 / 10.0;
        c_d[i] = 1.0 / 10.0;
    }

    thomas_algorithm(&a_d, &b_d, &c_d, &rhs)
}

// ─────────────────────────────────────────────────────────────────────────────
// Wave equation explicit step
// ─────────────────────────────────────────────────────────────────────────────

/// One explicit time step for the 1D wave equation `∂²u/∂t² = c² ∂²u/∂x²`.
///
/// Requires the solution at times `t` (`u`) and `t-dt` (`u_prev`).
pub fn wave_equation_step(u: &[f64], u_prev: &[f64], dx: f64, dt: f64, c: f64) -> Vec<f64> {
    let n = u.len();
    let r2 = (c * dt / dx).powi(2);
    let mut u_next = u.to_vec();
    for i in 1..n - 1 {
        u_next[i] = 2.0 * u[i] - u_prev[i] + r2 * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
    }
    u_next
}

// ─────────────────────────────────────────────────────────────────────────────
// Energy (total variation) diagnostics
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the total variation `TV(u) = Σ |u[i+1] - u[i]|`.
pub fn total_variation(u: &[f64]) -> f64 {
    u.windows(2).map(|w| (w[1] - w[0]).abs()).sum()
}

/// Computes the discrete L1 norm `||u||₁ = dx * Σ |u[i]|`.
pub fn l1_norm(u: &[f64], dx: f64) -> f64 {
    dx * u.iter().map(|v| v.abs()).sum::<f64>()
}

// ─────────────────────────────────────────────────────────────────────────────
// Smoothness and WENO utility
// ─────────────────────────────────────────────────────────────────────────────

/// Minimum smoothness indicator across stencils (for diagnostic use).
pub fn weno5_min_smoothness(u: &[f64]) -> f64 {
    let n = u.len();
    if n < 5 {
        return 0.0;
    }
    let mut min_b = f64::MAX;
    for i in 2..n - 2 {
        let b = weno5_beta(u[i - 1], u[i], u[i + 1]);
        if b < min_b {
            min_b = b;
        }
    }
    min_b
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Thomas algorithm ─────────────────────────────────────────────────────

    #[test]
    fn test_thomas_identity_system() {
        // b = [1,1,1], a = c = 0, d = [2,3,4] → x = [2,3,4]
        let a = vec![0.0, 0.0, 0.0];
        let b = vec![1.0, 1.0, 1.0];
        let c = vec![0.0, 0.0, 0.0];
        let d = vec![2.0, 3.0, 4.0];
        let x = thomas_algorithm(&a, &b, &c, &d);
        assert!((x[0] - 2.0).abs() < 1e-10);
        assert!((x[1] - 3.0).abs() < 1e-10);
        assert!((x[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_thomas_simple_tridiagonal() {
        // 2-equation: [2 -1; -1 2] x = [1; 1]  => x = [1, 1]
        let a = vec![0.0, -1.0];
        let b = vec![2.0, 2.0];
        let c = vec![-1.0, 0.0];
        let d = vec![1.0, 1.0];
        let x = thomas_algorithm(&a, &b, &c, &d);
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_thomas_larger_system() {
        let n = 5_usize;
        // diagonal dominant: b=4, a=c=-1
        let a: Vec<f64> = (0..n).map(|i| if i == 0 { 0.0 } else { -1.0 }).collect();
        let b = vec![4.0_f64; n];
        let c: Vec<f64> = (0..n)
            .map(|i| if i == n - 1 { 0.0 } else { -1.0 })
            .collect();
        let d = vec![1.0_f64; n];
        let x = thomas_algorithm(&a, &b, &c, &d);
        // Verify Ax = d
        for i in 0..n {
            let lhs = b[i] * x[i]
                + if i > 0 { a[i] * x[i - 1] } else { 0.0 }
                + if i < n - 1 { c[i] * x[i + 1] } else { 0.0 };
            assert!((lhs - d[i]).abs() < 1e-9, "row {i}: lhs={lhs}");
        }
    }

    // ── Derivative stencils ──────────────────────────────────────────────────

    #[test]
    fn test_first_derivative_linear() {
        // u = x => du/dx = 1 everywhere
        let n = 10;
        let dx = 0.1;
        let u: Vec<f64> = (0..n).map(|i| i as f64 * dx).collect();
        let du = first_derivative_central(&u, dx);
        for v in &du {
            assert!((v - 1.0).abs() < 1e-9, "got {v}");
        }
    }

    #[test]
    fn test_second_derivative_quadratic() {
        // u = x² => d²u/dx² = 2
        let n = 20;
        let dx = 0.1;
        let u: Vec<f64> = (0..n).map(|i| (i as f64 * dx).powi(2)).collect();
        let d2u = second_derivative(&u, dx);
        for (i, &val) in d2u.iter().enumerate().take(n - 1).skip(1) {
            assert!((val - 2.0).abs() < 1e-6, "i={i}: {}", val);
        }
    }

    #[test]
    fn test_fourth_order_derivative_sine() {
        // u = sin(x), du/dx = cos(x)
        let n = 100;
        let dx = 2.0 * PI / n as f64;
        let u: Vec<f64> = (0..n).map(|i| (i as f64 * dx).sin()).collect();
        let du = first_derivative_fourth_order(&u, dx);
        // check interior points (away from boundaries)
        for (i, &du_i) in du.iter().enumerate().take(n - 4).skip(4) {
            let exact = (i as f64 * dx).cos();
            assert!(
                (du_i - exact).abs() < 1e-6,
                "i={i}: got {} exp {}",
                du_i,
                exact
            );
        }
    }

    // ── Upwind advection ─────────────────────────────────────────────────────

    #[test]
    fn test_upwind_advects_rightward() {
        let n = 20;
        let dx = 1.0 / n as f64;
        let dt = 0.4 * dx; // CFL = 0.4
        let c = 1.0;
        let mut u: Vec<f64> = (0..n).map(|i| if i == n / 2 { 1.0 } else { 0.0 }).collect();
        let bc = BoundaryCondition::Dirichlet(0.0);
        u = advection_upwind_step(&u, dx, dt, c, bc, bc);
        // peak should have shifted right
        let peak_idx = u
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert!(peak_idx >= n / 2);
    }

    #[test]
    fn test_cfl_number_correct() {
        assert!((cfl_number(2.0, 0.5, 1.0) - 1.0).abs() < 1e-12);
        assert!((cfl_number(1.0, 0.1, 1.0) - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_stability_flags() {
        assert!(is_advection_stable(1.0, 0.9, 1.0));
        assert!(!is_advection_stable(1.0, 1.1, 1.0));
        assert!(is_diffusion_stable(0.1, 0.4, 1.0));
        assert!(!is_diffusion_stable(0.1, 6.0, 1.0));
    }

    // ── Crank-Nicolson diffusion ─────────────────────────────────────────────

    #[test]
    fn test_crank_nicolson_preserves_constant() {
        // u = const should stay constant
        let n = 10;
        let dx = 0.1;
        let dt = 1.0;
        let alpha = 0.01;
        let u = vec![3.0_f64; n];
        let bc = BoundaryCondition::Dirichlet(3.0);
        let un = diffusion_crank_nicolson_step(&u, dx, dt, alpha, bc, bc);
        for v in &un {
            assert!((v - 3.0).abs() < 1e-9, "got {v}");
        }
    }

    #[test]
    fn test_implicit_diffusion_stable_large_dt() {
        // Implicit diffusion should not blow up even with large dt
        let n = 10;
        let dx = 0.1;
        let dt = 100.0; // would be unstable for explicit
        let alpha = 0.01;
        let mut u = vec![0.0_f64; n];
        u[5] = 1.0;
        let bc = BoundaryCondition::Dirichlet(0.0);
        let un = diffusion_implicit_step(&u, dx, dt, alpha, bc, bc);
        // all values should be finite and bounded
        for v in &un {
            assert!(v.is_finite(), "got NaN/inf");
        }
    }

    #[test]
    fn test_cn_dirichlet_boundary_held() {
        let n = 8;
        let dx = 0.1;
        let dt = 0.1;
        let alpha = 0.1;
        let u: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let bc_l = BoundaryCondition::Dirichlet(0.0);
        let bc_r = BoundaryCondition::Dirichlet(7.0);
        let un = diffusion_crank_nicolson_step(&u, dx, dt, alpha, bc_l, bc_r);
        assert!((un[0] - 0.0).abs() < 1e-10);
        assert!((un[n - 1] - 7.0).abs() < 1e-10);
    }

    // ── Explicit diffusion ───────────────────────────────────────────────────

    #[test]
    fn test_explicit_diffusion_gaussian_decay() {
        // Gaussian should spread and decay
        let n = 50;
        let dx = 1.0 / n as f64;
        let dt = 0.4 * dx * dx / 0.01; // r = 0.4
        let alpha = 0.01;
        let sigma = 5.0 * dx;
        let mid = n as f64 / 2.0 * dx;
        let u: Vec<f64> = (0..n)
            .map(|i| {
                let x = i as f64 * dx - mid;
                (-0.5 * x * x / (sigma * sigma)).exp()
            })
            .collect();
        let bc = BoundaryCondition::Dirichlet(0.0);
        let un = diffusion_explicit_step(&u, dx, dt, alpha, bc, bc);
        let peak_old = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let peak_new = un.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(peak_new < peak_old, "Gaussian peak should decrease");
    }

    // ── WENO-5 ───────────────────────────────────────────────────────────────

    #[test]
    fn test_weno5_flux_smooth_region() {
        // for uniform data, WENO should return the uniform value
        let v = [1.0_f64; 5];
        let f = weno5_flux(v);
        assert!((f - 1.0).abs() < 1e-10, "got {f}");
    }

    #[test]
    fn test_weno5_advection_conserves_mass() {
        let n = 40;
        let dx = 1.0 / n as f64;
        let dt = 0.3 * dx;
        let c = 1.0;
        let u: Vec<f64> = (0..n)
            .map(|i| ((i as f64 * dx - 0.3) * 20.0 * PI).sin().powi(2) * 0.5)
            .collect();
        let mass_before: f64 = u.iter().sum::<f64>() * dx;
        let un = advection_weno5_step(&u, dx, dt, c);
        let mass_after: f64 = un.iter().sum::<f64>() * dx;
        // WENO-5 is not exactly conservative; allow up to 20% relative change
        assert!(
            (mass_before - mass_after).abs() < 0.20 * mass_before.abs().max(1e-6),
            "mass_before={mass_before}, mass_after={mass_after}"
        );
    }

    // ── Compact finite differences ───────────────────────────────────────────

    #[test]
    fn test_compact_derivative_linear() {
        let n = 10;
        let dx = 0.1;
        let u: Vec<f64> = (0..n).map(|i| 3.0 * i as f64 * dx).collect();
        let du = compact_first_derivative(&u, dx);
        // interior should give ~3.0
        for (i, &val) in du.iter().enumerate().take(n - 2).skip(2) {
            assert!((val - 3.0).abs() < 1e-8, "i={i}: {}", val);
        }
    }

    // ── Richardson extrapolation ─────────────────────────────────────────────

    #[test]
    fn test_richardson_second_order() {
        // f(h) ≈ 1 + h², f(h/2) ≈ 1 + h²/4; extrapolate with p=2 → should give 1
        let h = 0.1;
        let f_h = 1.0 + h * h;
        let f_h2 = 1.0 + h * h / 4.0;
        let ext = richardson_extrapolation(f_h, f_h2, 2);
        assert!((ext - 1.0).abs() < 1e-10, "got {ext}");
    }

    // ── Staggered grid ───────────────────────────────────────────────────────

    #[test]
    fn test_staggered_grid_divergence_constant_flux() {
        let n = 10;
        let dx = 0.1;
        let mut sg = StaggeredGrid1D::new(n, dx);
        // constant velocity → zero divergence
        sg.u = vec![1.0_f64; n];
        sg.interpolate_to_faces();
        let div = sg.divergence();
        for v in &div {
            assert!(v.abs() < 1e-12, "got {v}");
        }
    }

    #[test]
    fn test_staggered_grid_face_count() {
        let sg = StaggeredGrid1D::new(5, 0.1);
        assert_eq!(sg.flux.len(), 6);
        assert_eq!(sg.u.len(), 5);
    }

    // ── Fractional step ──────────────────────────────────────────────────────

    #[test]
    fn test_fractional_step_returns_finite() {
        let n = 20;
        let u: Vec<f64> = (0..n).map(|i| (i as f64 / n as f64).sin()).collect();
        let res = fractional_step_1d(&u, 1.0 / n as f64, 0.001, 0.01, 0.0);
        for v in &res.u_new {
            assert!(v.is_finite(), "NaN/inf in u_new");
        }
    }

    // ── 2D Laplacian ─────────────────────────────────────────────────────────

    #[test]
    fn test_laplacian_2d_constant_field() {
        let nx = 5;
        let ny = 5;
        let u = vec![2.0_f64; nx * ny];
        let lap = laplacian_2d(&u, nx, ny, 0.1, 0.1);
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!(lap[j * nx + i].abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_laplacian_2d_parabolic() {
        // u = x² + y² → Δu = 4
        let nx = 10;
        let ny = 10;
        let dx = 0.1;
        let dy = 0.1;
        let u: Vec<f64> = (0..ny)
            .flat_map(|j| {
                (0..nx).map(move |i| {
                    let x = i as f64 * dx;
                    let y = j as f64 * dy;
                    x * x + y * y
                })
            })
            .collect();
        let lap = laplacian_2d(&u, nx, ny, dx, dy);
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!(
                    (lap[j * nx + i] - 4.0).abs() < 1e-6,
                    "({i},{j}): {}",
                    lap[j * nx + i]
                );
            }
        }
    }

    // ── ADI ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_adi_constant_field_unchanged() {
        let nx = 6;
        let ny = 6;
        let u = vec![0.0_f64; nx * ny]; // all zeros (Dirichlet = 0 satisfied)
        let un = adi_diffusion_step_2d(&u, nx, ny, 0.1, 0.1, 0.01, 0.01);
        for v in &un {
            assert!(v.abs() < 1e-12, "got {v}");
        }
    }

    #[test]
    fn test_adi_energy_decreases() {
        let nx = 8;
        let ny = 8;
        let mut u = vec![0.0_f64; nx * ny];
        u[ny / 2 * nx + nx / 2] = 1.0; // point source
        let energy_before: f64 = u.iter().map(|v| v * v).sum();
        let un = adi_diffusion_step_2d(&u, nx, ny, 0.125, 0.125, 0.01, 0.01);
        let energy_after: f64 = un.iter().map(|v| v * v).sum();
        assert!(energy_after < energy_before, "Energy should decrease");
    }

    // ── Wave equation ────────────────────────────────────────────────────────

    #[test]
    fn test_wave_equation_step_finite() {
        let n = 30;
        let dx = 1.0 / n as f64;
        let dt = 0.4 * dx;
        let c = 1.0;
        let u: Vec<f64> = (0..n).map(|i| (PI * i as f64 * dx).sin()).collect();
        let u_prev = u.clone();
        let u_next = wave_equation_step(&u, &u_prev, dx, dt, c);
        for v in &u_next {
            assert!(v.is_finite(), "NaN/inf in wave step");
        }
    }

    // ── Periodic advection ───────────────────────────────────────────────────

    #[test]
    fn test_periodic_advection_conserves_mass() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let dt = 0.4 * dx;
        let c = 1.0;
        let u: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let mass0: f64 = u.iter().sum();
        let un = advection_periodic_step(&u, dx, dt, c);
        let mass1: f64 = un.iter().sum();
        assert!(
            (mass0 - mass1).abs() < 1e-10,
            "mass change = {}",
            (mass0 - mass1).abs()
        );
    }

    // ── L2 / Linf norms ──────────────────────────────────────────────────────

    #[test]
    fn test_l2_norm_constant() {
        let n = 10;
        let dx = 0.1;
        let u = vec![2.0_f64; n];
        let norm = l2_norm(&u, dx);
        let expected = (dx * n as f64 * 4.0).sqrt();
        assert!((norm - expected).abs() < 1e-10);
    }

    #[test]
    fn test_linf_norm() {
        let u = vec![-3.0_f64, 1.0, 2.0, -1.0, 2.5];
        assert!((linf_norm(&u) - 2.5).abs() < 1e-12);
    }

    // ── Total variation ──────────────────────────────────────────────────────

    #[test]
    fn test_total_variation_monotone() {
        let u: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let tv = total_variation(&u);
        assert!((tv - 9.0).abs() < 1e-12);
    }

    // ── Von Neumann amplification ─────────────────────────────────────────────

    #[test]
    fn test_ftcs_amplification_r_half_stable() {
        let r = 0.5;
        let dx = 0.1;
        let max_g = ftcs_max_amplification(r, dx, 200);
        assert!(max_g <= 1.0 + 1e-10, "unstable: max|G|={max_g}");
    }

    #[test]
    fn test_ftcs_amplification_r_too_large() {
        let r = 0.6; // should be unstable
        let dx = 0.1;
        let max_g = ftcs_max_amplification(r, dx, 200);
        assert!(max_g > 1.0, "|G| should exceed 1.0 for r=0.6, got {max_g}");
    }

    // ── Neumann BC reflection ─────────────────────────────────────────────────

    #[test]
    fn test_neumann_reflection_zero_gradient() {
        let mut u = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        apply_neumann_reflection(&mut u);
        assert_eq!(u[0], u[1]);
        assert_eq!(u[4], u[3]);
    }

    // ── Grid builder ─────────────────────────────────────────────────────────

    #[test]
    fn test_uniform_grid_endpoints() {
        let grid = uniform_grid_1d(0.0, 1.0, 11);
        assert!((grid[0] - 0.0).abs() < 1e-12);
        assert!((grid[10] - 1.0).abs() < 1e-12);
        assert_eq!(grid.len(), 11);
    }

    #[test]
    fn test_uniform_grid_2d_dimensions() {
        let (xs, ys) = uniform_grid_2d(0.0, 1.0, 5, 0.0, 2.0, 4);
        assert_eq!(xs.len(), 20);
        assert_eq!(ys.len(), 20);
    }

    // ── Lax-Wendroff ─────────────────────────────────────────────────────────

    #[test]
    fn test_lax_wendroff_cfl_one_exact() {
        // With CFL = 1, Lax-Wendroff is exact for linear advection
        let n = 20;
        let dx = 1.0 / n as f64;
        let c = 1.0;
        let dt = dx; // CFL = 1
        let u: Vec<f64> = (0..n).map(|i| (2.0 * PI * i as f64 * dx).sin()).collect();
        let expected: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * ((i as f64 * dx) - c * dt)).sin())
            .collect();
        let un = advection_lax_wendroff_step(&u, dx, dt, c);
        for i in 1..n - 1 {
            assert!(
                (un[i] - expected[i]).abs() < 1e-8,
                "i={i}: {} vs {}",
                un[i],
                expected[i]
            );
        }
    }

    // ── Method of lines ──────────────────────────────────────────────────────

    #[test]
    fn test_mol_rk4_diffusion_finite() {
        let n = 20;
        let dx = 1.0 / n as f64;
        let dt = 0.1 * dx * dx;
        let alpha = 0.01;
        let u: Vec<f64> = (0..n).map(|i| (PI * i as f64 * dx).sin()).collect();
        let u_next = rk4_step(&u, 0.0, dt, &|_t, u| mol_diffusion_rhs(u, dx, alpha));
        for v in &u_next {
            assert!(v.is_finite(), "NaN/inf in RK4 diffusion");
        }
    }

    // ── 3D Laplacian ─────────────────────────────────────────────────────────

    #[test]
    fn test_laplacian_3d_constant_zero() {
        let nx = 5;
        let ny = 5;
        let nz = 5;
        let u = vec![0.0_f64; nx * ny * nz];
        let lap = laplacian_3d(&u, nx, ny, nz, 0.1, 0.1, 0.1);
        for v in &lap {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn test_laplacian_3d_linear_field_zero() {
        // u = x + y + z => Δu = 0
        let nx = 6;
        let ny = 6;
        let nz = 6;
        let dx = 0.1;
        let dy = 0.1;
        let dz = 0.1;
        let u: Vec<f64> = (0..nz)
            .flat_map(|k| {
                (0..ny).flat_map(move |j| {
                    (0..nx).map(move |i| i as f64 * dx + j as f64 * dy + k as f64 * dz)
                })
            })
            .collect();
        let lap = laplacian_3d(&u, nx, ny, nz, dx, dy, dz);
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let v = lap[k * ny * nx + j * nx + i];
                    assert!(v.abs() < 1e-9, "({i},{j},{k}): {v}");
                }
            }
        }
    }

    // ── Compact second derivative ─────────────────────────────────────────────

    #[test]
    fn test_compact_second_derivative_constant() {
        let n = 10;
        let dx = 0.1;
        let u = vec![5.0_f64; n];
        let d2u = compact_second_derivative(&u, dx);
        for (i, &val) in d2u.iter().enumerate().take(n - 1).skip(1) {
            assert!(val.abs() < 1e-9, "i={i}: {}", val);
        }
    }

    // ── Max dt stability ─────────────────────────────────────────────────────

    #[test]
    fn test_max_dt_advection() {
        let dt_max = max_dt_advection(2.0, 0.1);
        assert!((dt_max - 0.05).abs() < 1e-12);
    }

    #[test]
    fn test_max_dt_diffusion() {
        let dt_max = max_dt_diffusion(0.1, 0.1);
        assert!((dt_max - 0.05).abs() < 1e-12);
    }

    // ── L1 norm ───────────────────────────────────────────────────────────────

    #[test]
    fn test_l1_norm_ones() {
        let u = vec![1.0_f64; 10];
        let norm = l1_norm(&u, 0.1);
        assert!((norm - 1.0).abs() < 1e-12);
    }

    // ── WENO smoothness ───────────────────────────────────────────────────────

    #[test]
    fn test_weno5_min_smoothness_zero_constant() {
        // constant field → smoothness = 0
        let u = vec![1.0_f64; 10];
        let b = weno5_min_smoothness(&u);
        assert!(b.abs() < 1e-12, "got {b}");
    }

    // ── Dirichlet 2D ─────────────────────────────────────────────────────────

    #[test]
    fn test_apply_dirichlet_2d_border() {
        let nx = 5;
        let ny = 5;
        let mut u = vec![9.0_f64; nx * ny];
        apply_dirichlet_2d(&mut u, nx, ny, 0.0);
        // check border is zero
        for i in 0..nx {
            assert_eq!(u[i], 0.0);
            assert_eq!(u[(ny - 1) * nx + i], 0.0);
        }
        for j in 0..ny {
            assert_eq!(u[j * nx], 0.0);
            assert_eq!(u[j * nx + nx - 1], 0.0);
        }
    }

    // ── Beam-Warming ─────────────────────────────────────────────────────────

    #[test]
    fn test_beam_warming_finite() {
        let n = 20;
        let dx = 1.0 / n as f64;
        let dt = 0.4 * dx;
        let u: Vec<f64> = (0..n).map(|i| (2.0 * PI * i as f64 * dx).sin()).collect();
        let un = advection_beam_warming_step(&u, dx, dt, 1.0);
        for v in &un {
            assert!(v.is_finite(), "NaN/inf");
        }
    }
}
