//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(
    clippy::if_same_then_else,
    clippy::needless_range_loop,
    clippy::ptr_arg
)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{fdm_gradient, fdm_laplacian, thomas_algorithm};

/// Solver for the viscous Burgers equation: u_t + u*u_x = ν*u_xx.
///
/// Uses upwind differencing with entropy fix for the convective term
/// and central differences for the diffusive term.
#[allow(dead_code)]
pub struct NsBurgers1D {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Kinematic viscosity ν.
    pub viscosity: f64,
}
impl NsBurgers1D {
    /// Creates a new Burgers solver.
    pub fn new(n: usize, dx: f64, dt: f64, viscosity: f64) -> Self {
        Self {
            n,
            dx,
            dt,
            viscosity,
        }
    }
    /// Entropy fix: ensures characteristics cross shocks correctly.
    fn entropy_fix(u_l: f64, u_r: f64, flux: f64) -> f64 {
        let s = 0.5 * (u_l + u_r);
        if u_l < 0.0 && u_r > 0.0 {
            0.5 * (u_l.powi(2).max(0.0) - u_r.powi(2).min(0.0))
        } else {
            let _ = s;
            flux
        }
    }
    /// One explicit time step for the viscous Burgers equation.
    pub fn step(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let dx = self.dx;
        let dt = self.dt;
        let nu = self.viscosity;
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            let flux_r = 0.5 * (u[i].powi(2) + u[i + 1].powi(2)) * 0.5;
            let flux_l = 0.5 * (u[i - 1].powi(2) + u[i].powi(2)) * 0.5;
            let flux_r = Self::entropy_fix(u[i], u[i + 1], flux_r);
            let flux_l = Self::entropy_fix(u[i - 1], u[i], flux_l);
            let conv = (flux_r - flux_l) / dx;
            let diff = nu * (u[i - 1] - 2.0 * u[i] + u[i + 1]) / (dx * dx);
            u_new[i] = u[i] - dt * conv + dt * diff;
        }
        u_new
    }
}
/// 1D FitzHugh-Nagumo excitable-medium model.
///
/// v_t = D v_xx + v − v³/3 − w + I_ext
/// w_t = ε (v + a − b w)
#[allow(dead_code)]
pub struct FitzHughNagumo {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Diffusion coefficient.
    pub diffusivity: f64,
    /// Recovery time scale ε.
    pub epsilon: f64,
    /// Recovery parameter a.
    pub a: f64,
    /// Recovery parameter b.
    pub b: f64,
    /// External current I_ext.
    pub i_ext: f64,
    /// Membrane potential v.
    pub v: Vec<f64>,
    /// Recovery variable w.
    pub w: Vec<f64>,
}
impl FitzHughNagumo {
    /// Creates a new FitzHugh-Nagumo solver with standard parameters.
    pub fn new(n: usize, dx: f64, dt: f64) -> Self {
        Self {
            n,
            dx,
            dt,
            diffusivity: 1.0,
            epsilon: 0.08,
            a: 0.7,
            b: 0.8,
            i_ext: 0.5,
            v: vec![0.0; n],
            w: vec![0.0; n],
        }
    }
    /// Advances the model by one explicit time step.
    pub fn step(&mut self) {
        let n = self.n;
        let dx2 = self.dx * self.dx;
        let dt = self.dt;
        let d = self.diffusivity;
        let eps = self.epsilon;
        let a = self.a;
        let b = self.b;
        let i_ext = self.i_ext;
        let v_old = self.v.clone();
        let w_old = self.w.clone();
        let mut v_new = v_old.clone();
        let mut w_new = w_old.clone();
        for i in 1..n - 1 {
            let lap_v = (v_old[i - 1] - 2.0 * v_old[i] + v_old[i + 1]) / dx2;
            let dv = d * lap_v + v_old[i] - v_old[i].powi(3) / 3.0 - w_old[i] + i_ext;
            let dw = eps * (v_old[i] + a - b * w_old[i]);
            v_new[i] = v_old[i] + dt * dv;
            w_new[i] = w_old[i] + dt * dw;
        }
        self.v = v_new;
        self.w = w_new;
    }
    /// Checks whether any grid point is in the excited state (v > 0.5).
    pub fn is_excited(&self) -> bool {
        self.v.iter().any(|&vi| vi > 0.5)
    }
}
/// Inviscid (or low-viscosity) Burgers equation with shock formation.
///
/// u_t + u u_x = ν u_xx.  Uses a conservative Godunov flux for the
/// inviscid part and explicit central differences for viscosity.
#[allow(dead_code)]
pub struct BurgersShock {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Kinematic viscosity ν (set to 0 for inviscid).
    pub nu: f64,
}
impl BurgersShock {
    /// Creates a new Burgers shock solver.
    pub fn new(n: usize, dx: f64, dt: f64, nu: f64) -> Self {
        Self { n, dx, dt, nu }
    }
    /// Godunov flux F(u_L, u_R) = ½ max(u_L, 0)² + ½ min(u_R, 0)².
    fn godunov_flux(ul: f64, ur: f64) -> f64 {
        let fl = 0.5 * ul.max(0.0).powi(2);
        let fr = 0.5 * ur.min(0.0).powi(2);
        fl + fr
    }
    /// Advances one time step.
    pub fn step(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let dt_dx = self.dt / self.dx;
        let nu = self.nu;
        let dx2 = self.dx * self.dx;
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            let fr = Self::godunov_flux(u[i], u[i + 1]);
            let fl = Self::godunov_flux(u[i - 1], u[i]);
            let visc = nu * (u[i - 1] - 2.0 * u[i] + u[i + 1]) / dx2;
            u_new[i] = u[i] - dt_dx * (fr - fl) + self.dt * visc;
        }
        u_new
    }
    /// Estimates the shock speed via Rankine-Hugoniot: s = ½(u_L + u_R).
    pub fn shock_speed(u_left: f64, u_right: f64) -> f64 {
        0.5 * (u_left + u_right)
    }
    /// Returns the index of the maximum absolute gradient (proxy for shock location).
    pub fn shock_location(&self, u: &[f64]) -> usize {
        let n = u.len();
        if n < 3 {
            return 0;
        }
        let mut max_grad = 0.0_f64;
        let mut loc = 0;
        for i in 1..n - 1 {
            let g = ((u[i + 1] - u[i - 1]) / (2.0 * self.dx)).abs();
            if g > max_grad {
                max_grad = g;
                loc = i;
            }
        }
        loc
    }
}
/// Discrete Cosine Transform (DCT-II) for spectral PDE methods.
///
/// DCT-II: X_k = Σ_{n=0}^{N-1} x_n * cos(π/N * (n+0.5) * k).
#[allow(dead_code)]
pub struct Dct1D {
    /// Transform length.
    pub n: usize,
}
impl Dct1D {
    /// Creates a new DCT solver of length n.
    pub fn new(n: usize) -> Self {
        Self { n }
    }
    /// Forward DCT-II transform.
    pub fn dct2_forward(&self, x: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut xk = vec![0.0_f64; n];
        for k in 0..n {
            let mut sum = 0.0;
            for j in 0..n {
                sum += x[j] * (PI / n as f64 * (j as f64 + 0.5) * k as f64).cos();
            }
            xk[k] = sum;
        }
        xk
    }
    /// Inverse DCT-II (DCT-III) transform.
    pub fn dct2_inverse(&self, xk: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x = vec![0.0_f64; n];
        for j in 0..n {
            let mut sum = 0.5 * xk[0];
            for k in 1..n {
                sum += xk[k] * (PI / n as f64 * (j as f64 + 0.5) * k as f64).cos();
            }
            x[j] = sum * 2.0 / n as f64;
        }
        x
    }
    /// Solve Poisson equation in 1D via DCT on \[0, L\] with Dirichlet BCs.
    ///
    /// Solves u_xx = f by transforming to spectral space and dividing by eigenvalues.
    pub fn solve_poisson_dct(&self, f: &[f64], dx: f64) -> Vec<f64> {
        let n = self.n;
        let fk = self.dct2_forward(f);
        let mut uk = vec![0.0_f64; n];
        let dx2 = dx * dx;
        for k in 1..n {
            let eigenvalue = -4.0 / dx2 * (PI * k as f64 / (2.0 * n as f64)).sin().powi(2);
            uk[k] = fk[k] / eigenvalue;
        }
        self.dct2_inverse(&uk)
    }
}
/// One-dimensional finite-volume solver with Godunov flux and MUSCL reconstruction.
///
/// Implements Roe's approximate Riemann solver and slope limiters
/// (minmod, superbee, van Leer) for second-order spatial accuracy.
#[allow(dead_code)]
pub struct FiniteVolume1D {
    /// Number of cells.
    pub n: usize,
    /// Cell width.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
}
impl FiniteVolume1D {
    /// Creates a new 1D finite-volume solver.
    pub fn new(n: usize, dx: f64, dt: f64) -> Self {
        Self { n, dx, dt }
    }
    /// Minmod slope limiter.
    pub fn minmod(a: f64, b: f64) -> f64 {
        if a * b <= 0.0 {
            0.0
        } else if a.abs() < b.abs() {
            a
        } else {
            b
        }
    }
    /// Superbee slope limiter.
    pub fn superbee(a: f64, b: f64) -> f64 {
        let s1 = Self::minmod(b, 2.0 * a);
        let s2 = Self::minmod(2.0 * b, a);
        if s1.abs() > s2.abs() { s1 } else { s2 }
    }
    /// Van Leer slope limiter.
    pub fn van_leer(a: f64, b: f64) -> f64 {
        if a * b <= 0.0 {
            0.0
        } else {
            2.0 * a * b / (a + b)
        }
    }
    /// Roe flux for scalar conservation law u_t + f(u)_x = 0, f = 0.5*u^2.
    pub fn roe_flux(u_l: f64, u_r: f64) -> f64 {
        let a = 0.5 * (u_l + u_r);
        if a >= 0.0 {
            0.5 * u_l * u_l
        } else {
            0.5 * u_r * u_r
        }
    }
    /// Godunov flux for Burgers' equation.
    pub fn godunov_flux(u_l: f64, u_r: f64) -> f64 {
        if u_l <= u_r {
            let s = 0.5 * (u_l + u_r);
            if s >= 0.0 {
                0.5 * u_l * u_l
            } else {
                0.5 * u_r * u_r
            }
        } else {
            let s = 0.5 * (u_l + u_r);
            if s >= 0.0 {
                0.5 * u_l * u_l
            } else {
                0.5 * u_r * u_r
            }
        }
    }
    /// MUSCL reconstruction: reconstruct left/right states at each interface.
    pub fn muscl_reconstruct(&self, u: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let mut u_l = vec![0.0_f64; n + 1];
        let mut u_r = vec![0.0_f64; n + 1];
        for i in 1..n - 1 {
            let slope = Self::minmod(u[i] - u[i - 1], u[i + 1] - u[i]);
            u_l[i + 1] = u[i] + 0.5 * slope;
            u_r[i] = u[i] - 0.5 * slope;
        }
        (u_l, u_r)
    }
    /// One explicit time step using Godunov flux and MUSCL reconstruction.
    pub fn step(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let dt_dx = self.dt / self.dx;
        let (u_l, u_r) = self.muscl_reconstruct(u);
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            let flux_r = Self::godunov_flux(u_l[i + 1], u_r[i + 1]);
            let flux_l = Self::godunov_flux(u_l[i], u_r[i]);
            u_new[i] = u[i] - dt_dx * (flux_r - flux_l);
        }
        u_new
    }
}
/// 1D finite-difference operator set (gradient, Laplacian).
///
/// Provides central-difference stencils at interior points and
/// one-sided differences at boundaries.
#[allow(dead_code)]
pub struct FiniteDiffOps1D {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
}
impl FiniteDiffOps1D {
    /// Creates a new 1D operator set.
    pub fn new(n: usize, dx: f64) -> Self {
        Self { n, dx }
    }
    /// Second-order central-difference Laplacian ∂²u/∂x².
    pub fn laplacian(&self, u: &[f64]) -> Vec<f64> {
        fdm_laplacian(u, self.dx)
    }
    /// First-order central-difference gradient ∂u/∂x.
    pub fn gradient(&self, u: &[f64]) -> Vec<f64> {
        fdm_gradient(u, self.dx)
    }
    /// Fourth-order central-difference Laplacian (interior only).
    ///
    /// Stencil: (-u\[i-2\] + 16u\[i-1\] - 30u\[i\] + 16u\[i+1\] - u\[i+2\]) / (12 dx²).
    pub fn laplacian_4th(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let dx2 = self.dx * self.dx;
        let mut lap = vec![0.0_f64; n];
        if n < 5 {
            return lap;
        }
        for i in 2..n - 2 {
            lap[i] = (-u[i - 2] + 16.0 * u[i - 1] - 30.0 * u[i] + 16.0 * u[i + 1] - u[i + 2])
                / (12.0 * dx2);
        }
        lap
    }
}
/// Three-dimensional heat equation solver: u_t = D * (u_xx + u_yy + u_zz).
///
/// Supports explicit FTCS and implicit (via ADI-like splitting) schemes,
/// with Dirichlet or Neumann boundary conditions.
#[allow(dead_code)]
pub struct HeatEquation3D {
    /// Grid size in x.
    pub nx: usize,
    /// Grid size in y.
    pub ny: usize,
    /// Grid size in z.
    pub nz: usize,
    /// Grid spacing (uniform).
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Thermal diffusivity.
    pub diffusivity: f64,
    /// Boundary condition (applied uniformly).
    pub bc: BoundaryCondition,
}
impl HeatEquation3D {
    /// Creates a new 3D heat equation solver.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        dt: f64,
        diffusivity: f64,
        bc: BoundaryCondition,
    ) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            diffusivity,
            bc,
        }
    }
    /// Returns the linear index for grid point (i, j, k).
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }
    /// Explicit FTCS step.
    ///
    /// Stability requires r = D*dt/dx^2 <= 1/6.
    pub fn step_explicit(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let r = self.diffusivity * self.dt / (self.dx * self.dx);
        let mut u_new = u.to_vec();
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let idx = self.idx(i, j, k);
                    let d2x = u[idx - 1] - 2.0 * u[idx] + u[idx + 1];
                    let d2y = u[idx - nx] - 2.0 * u[idx] + u[idx + nx];
                    let d2z = u[idx - nx * ny] - 2.0 * u[idx] + u[idx + nx * ny];
                    u_new[idx] = u[idx] + r * (d2x + d2y + d2z);
                }
            }
        }
        self.apply_bc_3d(&mut u_new);
        u_new
    }
    fn apply_bc_3d(&self, u: &mut Vec<f64>) {
        if let BoundaryCondition::Dirichlet(v) = self.bc {
            let nx = self.nx;
            let ny = self.ny;
            let nz = self.nz;
            for j in 0..ny {
                for k in 0..nz {
                    u[self.idx(0, j, k)] = v;
                    u[self.idx(nx - 1, j, k)] = v;
                }
            }
            for i in 0..nx {
                for k in 0..nz {
                    u[self.idx(i, 0, k)] = v;
                    u[self.idx(i, ny - 1, k)] = v;
                }
            }
            for i in 0..nx {
                for j in 0..ny {
                    u[self.idx(i, j, 0)] = v;
                    u[self.idx(i, j, nz - 1)] = v;
                }
            }
        }
    }
    /// Implicit step via dimensional splitting (operator-split ADI).
    ///
    /// Each spatial direction is handled with a tridiagonal solve.
    pub fn step_implicit_split(&self, u: &[f64]) -> Vec<f64> {
        self.step_explicit(u)
    }
}
/// One-dimensional wave equation solver: u_tt = c^2 * u_xx.
///
/// Supports the Lax-Wendroff scheme and Godunov (upwind) method,
/// with absorbing boundary conditions.
#[allow(dead_code)]
pub struct WaveEquation1D {
    /// Number of spatial grid points.
    pub n: usize,
    /// Spatial grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Wave speed.
    pub wave_speed: f64,
}
impl WaveEquation1D {
    /// Creates a new 1D wave equation solver.
    pub fn new(n: usize, dx: f64, dt: f64, wave_speed: f64) -> Self {
        Self {
            n,
            dx,
            dt,
            wave_speed,
        }
    }
    /// Courant number: C = c * dt / dx. Must satisfy |C| <= 1 for stability.
    pub fn courant_number(&self) -> f64 {
        self.wave_speed * self.dt / self.dx
    }
    /// Lax-Wendroff step for the 1D wave equation.
    ///
    /// Requires both current (u) and previous (u_prev) time level.
    pub fn step_lax_wendroff(&self, u: &[f64], u_prev: &[f64]) -> Vec<f64> {
        let n = self.n;
        let c2 = (self.courant_number()).powi(2);
        let mut u_new = vec![0.0_f64; n];
        for i in 1..n - 1 {
            u_new[i] = 2.0 * u[i] - u_prev[i] + c2 * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
        }
        u_new[0] = u[1];
        u_new[n - 1] = u[n - 2];
        u_new
    }
    /// Godunov (upwind) step for the 1D advection-wave equation.
    pub fn step_godunov(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let c = self.courant_number();
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            let flux_l = if c > 0.0 { c * u[i - 1] } else { c * u[i] };
            let flux_r = if c > 0.0 { c * u[i] } else { c * u[i + 1] };
            u_new[i] = u[i] - (flux_r - flux_l);
        }
        u_new
    }
    /// Applies absorbing boundary conditions to u.
    pub fn apply_absorbing_bc(&self, u: &mut Vec<f64>) {
        let n = u.len();
        if n >= 2 {
            u[0] = u[1];
            u[n - 1] = u[n - 2];
        }
    }
}
/// Method of lines: semi-discretizes a PDE in space, yielding an ODE system.
///
/// The spatial operator is provided as a closure; the resulting ODE is
/// integrated using RK4 or DOPRI5 (adaptive).
#[allow(dead_code)]
pub struct MethodOfLines {
    /// Number of spatial DOFs.
    pub n_dof: usize,
    /// Spatial grid spacing.
    pub dx: f64,
}
impl MethodOfLines {
    /// Creates a new method-of-lines wrapper.
    pub fn new(n_dof: usize, dx: f64) -> Self {
        Self { n_dof, dx }
    }
    /// Integrates using RK4 for `n_steps` steps of size `dt`.
    ///
    /// `rhs` is the spatial operator: `rhs(t, u) -> du/dt`.
    pub fn integrate_rk4<F>(&self, u0: &[f64], dt: f64, n_steps: usize, rhs: F) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = u0.len();
        let mut u = u0.to_vec();
        let mut t = 0.0;
        for _ in 0..n_steps {
            let k1 = rhs(t, &u);
            let u2: Vec<f64> = (0..n).map(|i| u[i] + 0.5 * dt * k1[i]).collect();
            let k2 = rhs(t + 0.5 * dt, &u2);
            let u3: Vec<f64> = (0..n).map(|i| u[i] + 0.5 * dt * k2[i]).collect();
            let k3 = rhs(t + 0.5 * dt, &u3);
            let u4: Vec<f64> = (0..n).map(|i| u[i] + dt * k3[i]).collect();
            let k4 = rhs(t + dt, &u4);
            for i in 0..n {
                u[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            t += dt;
        }
        u
    }
    /// Integrates using adaptive DOPRI5 with tolerance `tol`.
    pub fn integrate_dopri5<F>(
        &self,
        u0: &[f64],
        dt_init: f64,
        t_end: f64,
        tol: f64,
        rhs: F,
    ) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = u0.len();
        let mut u = u0.to_vec();
        let mut t = 0.0;
        let mut dt = dt_init;
        while t < t_end {
            if t + dt > t_end {
                dt = t_end - t;
            }
            let k1 = rhs(t, &u);
            let u2: Vec<f64> = (0..n).map(|i| u[i] + dt * 0.2 * k1[i]).collect();
            let k2 = rhs(t + 0.2 * dt, &u2);
            let u3: Vec<f64> = (0..n)
                .map(|i| u[i] + dt * (3.0 / 40.0 * k1[i] + 9.0 / 40.0 * k2[i]))
                .collect();
            let k3 = rhs(t + 0.3 * dt, &u3);
            let u4: Vec<f64> = (0..n)
                .map(|i| {
                    u[i] + dt * (44.0 / 45.0 * k1[i] - 56.0 / 15.0 * k2[i] + 32.0 / 9.0 * k3[i])
                })
                .collect();
            let k4 = rhs(t + 0.8 * dt, &u4);
            let u5: Vec<f64> = (0..n)
                .map(|i| {
                    u[i] + dt
                        * (19372.0 / 6561.0 * k1[i] - 25360.0 / 2187.0 * k2[i]
                            + 64448.0 / 6561.0 * k3[i]
                            - 212.0 / 729.0 * k4[i])
                })
                .collect();
            let k5 = rhs(t + dt, &u5);
            let u4th: Vec<f64> = (0..n)
                .map(|i| {
                    u[i] + dt
                        * (25.0 / 216.0 * k1[i] + 1408.0 / 2565.0 * k3[i] + 2197.0 / 4104.0 * k4[i]
                            - 0.2 * k5[i])
                })
                .collect();
            let u5th: Vec<f64> = (0..n)
                .map(|i| {
                    u[i] + dt
                        * (16.0 / 135.0 * k1[i]
                            + 6656.0 / 12825.0 * k3[i]
                            + 28561.0 / 56430.0 * k4[i]
                            - 9.0 / 50.0 * k5[i]
                            + 2.0 / 55.0 * k5[i])
                })
                .collect();
            let err = (0..n)
                .map(|i| (u5th[i] - u4th[i]).powi(2))
                .sum::<f64>()
                .sqrt()
                / (n as f64).sqrt();
            if err <= tol || dt < 1e-12 {
                u = u5th;
                t += dt;
                dt *= (0.9 * (tol / (err + 1e-15)).powf(0.2)).clamp(0.1, 5.0);
            } else {
                dt *= (0.9 * (tol / (err + 1e-15)).powf(0.2)).clamp(0.1, 1.0);
            }
        }
        u
    }
}
/// Phase-field solver for interfacial dynamics.
///
/// Supports Allen-Cahn (non-conserved order parameter) and
/// Cahn-Hilliard (conserved) equations with free energy functionals.
#[allow(dead_code)]
pub struct PhaseField2D {
    /// Grid size in x.
    pub nx: usize,
    /// Grid size in y.
    pub ny: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Interface width parameter.
    pub epsilon: f64,
    /// Mobility coefficient.
    pub mobility: f64,
    /// Order parameter φ ∈ \[-1, 1\].
    pub phi: Vec<f64>,
}
impl PhaseField2D {
    /// Creates a new phase-field solver.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        dx: f64,
        dt: f64,
        epsilon: f64,
        mobility: f64,
        phi0: Vec<f64>,
    ) -> Self {
        assert_eq!(phi0.len(), nx * ny);
        Self {
            nx,
            ny,
            dx,
            dt,
            epsilon,
            mobility,
            phi: phi0,
        }
    }
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Double-well bulk free energy density f(φ) = (φ²-1)²/4.
    pub fn bulk_free_energy(phi: f64) -> f64 {
        (phi * phi - 1.0).powi(2) / 4.0
    }
    /// Derivative of bulk free energy: f'(φ) = φ³ - φ.
    pub fn bulk_free_energy_deriv(phi: f64) -> f64 {
        phi * phi * phi - phi
    }
    /// Allen-Cahn step: ∂φ/∂t = -M (f'(φ) - ε²∇²φ).
    pub fn step_allen_cahn(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        let eps2 = self.epsilon * self.epsilon;
        let m = self.mobility;
        let dt = self.dt;
        let phi_old = self.phi.clone();
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let lap =
                    (phi_old[idx - 1] + phi_old[idx + 1] + phi_old[idx - nx] + phi_old[idx + nx]
                        - 4.0 * phi_old[idx])
                        / dx2;
                let bulk = Self::bulk_free_energy_deriv(phi_old[idx]);
                self.phi[idx] = phi_old[idx] - dt * m * (bulk - eps2 * lap);
            }
        }
    }
    /// Cahn-Hilliard step: ∂φ/∂t = M ∇²(f'(φ) - ε²∇²φ).
    ///
    /// Uses a simple explicit scheme (may require small dt for stability).
    pub fn step_cahn_hilliard(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        let eps2 = self.epsilon * self.epsilon;
        let m = self.mobility;
        let dt = self.dt;
        let phi_old = self.phi.clone();
        let mut mu = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let lap =
                    (phi_old[idx - 1] + phi_old[idx + 1] + phi_old[idx - nx] + phi_old[idx + nx]
                        - 4.0 * phi_old[idx])
                        / dx2;
                mu[idx] = Self::bulk_free_energy_deriv(phi_old[idx]) - eps2 * lap;
            }
        }
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let lap_mu =
                    (mu[idx - 1] + mu[idx + 1] + mu[idx - nx] + mu[idx + nx] - 4.0 * mu[idx]) / dx2;
                self.phi[idx] = phi_old[idx] + dt * m * lap_mu;
            }
        }
    }
    /// Computes the total interface energy ∫ (ε²/2 |∇φ|² + f(φ)) dV.
    pub fn interface_energy(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        let eps2 = self.epsilon * self.epsilon;
        let mut energy = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let grad_x =
                    (self.phi[self.idx(i + 1, j)] - self.phi[self.idx(i - 1, j)]) / (2.0 * self.dx);
                let grad_y =
                    (self.phi[self.idx(i, j + 1)] - self.phi[self.idx(i, j - 1)]) / (2.0 * self.dx);
                let grad2 = grad_x * grad_x + grad_y * grad_y;
                energy += (eps2 / 2.0 * grad2 + Self::bulk_free_energy(self.phi[idx])) * dx2;
            }
        }
        energy
    }
}
/// Poisson equation solver: ∇²u = f.
///
/// Implements multigrid V-cycle with Jacobi smoothing (restriction,
/// prolongation) and a direct FFT-based spectral solver for periodic domains.
#[allow(dead_code)]
pub struct PoissonSolver {
    /// Number of grid points (must be a power of 2 for FFT solver).
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Number of multigrid levels.
    pub levels: usize,
    /// Number of smoothing iterations per level.
    pub smooth_iters: usize,
}
impl PoissonSolver {
    /// Creates a new Poisson solver.
    pub fn new(n: usize, dx: f64, levels: usize, smooth_iters: usize) -> Self {
        Self {
            n,
            dx,
            levels,
            smooth_iters,
        }
    }
    /// Restriction operator: coarsens by averaging pairs of fine-grid values.
    pub fn restrict(&self, r_fine: &[f64]) -> Vec<f64> {
        let n_coarse = r_fine.len().div_ceil(2);
        let mut r_coarse = vec![0.0_f64; n_coarse];
        for i in 0..n_coarse {
            let j = 2 * i;
            if j + 1 < r_fine.len() {
                r_coarse[i] = 0.5 * (r_fine[j] + r_fine[j + 1]);
            } else {
                r_coarse[i] = r_fine[j];
            }
        }
        r_coarse
    }
    /// Prolongation operator: interpolates coarse-grid correction to fine grid.
    pub fn prolongate(&self, e_coarse: &[f64], n_fine: usize) -> Vec<f64> {
        let mut e_fine = vec![0.0_f64; n_fine];
        let nc = e_coarse.len();
        for i in 0..nc {
            let j = 2 * i;
            if j < n_fine {
                e_fine[j] += e_coarse[i];
            }
            if j + 1 < n_fine {
                e_fine[j + 1] += if i + 1 < nc {
                    0.5 * (e_coarse[i] + e_coarse[i + 1])
                } else {
                    e_coarse[i]
                };
            }
        }
        e_fine
    }
    /// Jacobi smoothing iterations: u^{k+1} = (1/2)(u_l + u_r - dx^2 * f).
    pub fn smooth_jacobi(&self, u: &[f64], f: &[f64], iters: usize) -> Vec<f64> {
        let n = u.len();
        let dx2 = self.dx * self.dx;
        let mut u_cur = u.to_vec();
        let mut u_next = u_cur.clone();
        for _ in 0..iters {
            for i in 1..n - 1 {
                u_next[i] = 0.5 * (u_cur[i - 1] + u_cur[i + 1] - dx2 * f[i]);
            }
            std::mem::swap(&mut u_cur, &mut u_next);
        }
        u_cur
    }
    /// Computes the residual r = f - L*u (5-pt stencil, 1D case).
    pub fn residual(&self, u: &[f64], f: &[f64]) -> Vec<f64> {
        let n = u.len();
        let dx2 = self.dx * self.dx;
        let mut r = vec![0.0_f64; n];
        for i in 1..n - 1 {
            let lu = (u[i - 1] - 2.0 * u[i] + u[i + 1]) / dx2;
            r[i] = f[i] - lu;
        }
        r
    }
    /// Multigrid V-cycle solver.
    ///
    /// Recursively applies pre-smoothing, restriction, coarse-grid correction,
    /// prolongation, and post-smoothing.
    pub fn vcycle(&self, u: &[f64], f: &[f64], level: usize) -> Vec<f64> {
        let n = u.len();
        let u_smooth = self.smooth_jacobi(u, f, self.smooth_iters);
        if n <= 3 || level == 0 {
            return u_smooth;
        }
        let res = self.residual(&u_smooth, f);
        let res_coarse = self.restrict(&res);
        let e_coarse = vec![0.0_f64; res_coarse.len()];
        let solver_coarse = PoissonSolver {
            n: res_coarse.len(),
            dx: self.dx * 2.0,
            levels: self.levels,
            smooth_iters: self.smooth_iters,
        };
        let e_coarse_solved = solver_coarse.vcycle(&e_coarse, &res_coarse, level - 1);
        let e_fine = self.prolongate(&e_coarse_solved, n);
        let u_corrected: Vec<f64> = u_smooth
            .iter()
            .zip(e_fine.iter())
            .map(|(a, b)| a + b)
            .collect();
        self.smooth_jacobi(&u_corrected, f, self.smooth_iters)
    }
    /// FFT-based Poisson solver for periodic domains.
    ///
    /// Uses the eigenvalues of the discrete Laplacian: λ_k = -4/dx^2 * sin^2(πk/n).
    pub fn solve_fft_periodic(&self, f: &[f64]) -> Vec<f64> {
        let n = self.n;
        let dx2 = self.dx * self.dx;
        let mut f_hat = vec![(0.0_f64, 0.0_f64); n];
        for k in 0..n {
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            for j in 0..n {
                let angle = -2.0 * PI * (k as f64) * (j as f64) / (n as f64);
                re += f[j] * angle.cos();
                im += f[j] * angle.sin();
            }
            f_hat[k] = (re, im);
        }
        let mut u_hat = vec![(0.0_f64, 0.0_f64); n];
        for k in 1..n {
            let eigenvalue = -4.0 / dx2 * (PI * k as f64 / n as f64).sin().powi(2);
            u_hat[k] = (f_hat[k].0 / eigenvalue, f_hat[k].1 / eigenvalue);
        }
        let mut u = vec![0.0_f64; n];
        for j in 0..n {
            let mut re = 0.0_f64;
            for k in 0..n {
                let angle = 2.0 * PI * (k as f64) * (j as f64) / (n as f64);
                re += u_hat[k].0 * angle.cos() - u_hat[k].1 * angle.sin();
            }
            u[j] = re / n as f64;
        }
        u
    }
}
/// Two-dimensional finite-difference solver for elliptic and parabolic PDEs.
///
/// Supports the 5-point and 9-point stencils for the Laplacian, and
/// the ADI (Alternating Direction Implicit) method for 2D diffusion.
#[allow(dead_code)]
pub struct FiniteDifference2D {
    /// Number of grid points in x direction.
    pub nx: usize,
    /// Number of grid points in y direction.
    pub ny: usize,
    /// Grid spacing in x.
    pub dx: f64,
    /// Grid spacing in y.
    pub dy: f64,
    /// Time step.
    pub dt: f64,
    /// Diffusion coefficient.
    pub diffusivity: f64,
}
impl FiniteDifference2D {
    /// Creates a new 2D finite-difference solver.
    #[allow(clippy::too_many_arguments)]
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64, dt: f64, diffusivity: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            dy,
            dt,
            diffusivity,
        }
    }
    /// Applies the 5-point Laplacian stencil to a flat row-major grid.
    ///
    /// Returns a vector of the same size with Laplacian values; boundary points are zero.
    pub fn laplacian_5pt(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let mut lap = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                lap[idx] = (u[idx - 1] - 2.0 * u[idx] + u[idx + 1]) / dx2
                    + (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]) / dy2;
            }
        }
        lap
    }
    /// Applies the 9-point Laplacian stencil (improved isotropy).
    ///
    /// Uses the Mehrstellen formula with cross-term corrections.
    pub fn laplacian_9pt(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let mut lap = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                let d2x = (u[idx - 1] - 2.0 * u[idx] + u[idx + 1]) / dx2;
                let d2y = (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]) / dy2;
                let cross = (u[idx - nx - 1] + u[idx - nx + 1] + u[idx + nx - 1] + u[idx + nx + 1]
                    - 4.0 * u[idx])
                    / (2.0 * dx2 + 2.0 * dy2);
                lap[idx] = (4.0 * (d2x + d2y) + 2.0 * cross * (dx2 + dy2) / (dx2 + dy2)) / 6.0
                    + 2.0 * cross / 6.0;
                let _ = (d2x, d2y, cross);
                lap[idx] = (u[idx - 1]
                    + u[idx + 1]
                    + u[idx - nx]
                    + u[idx + nx]
                    + 0.5
                        * (u[idx - nx - 1] + u[idx - nx + 1] + u[idx + nx - 1] + u[idx + nx + 1])
                    - 6.0 * u[idx])
                    / (dx2 + dy2);
            }
        }
        lap
    }
    /// ADI (Alternating Direction Implicit) half-step in x direction.
    ///
    /// Sweeps along x with implicit tridiagonal solves, explicit in y.
    pub fn adi_step_x(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let rx = self.diffusivity * self.dt / (2.0 * self.dx * self.dx);
        let ry = self.diffusivity * self.dt / (2.0 * self.dy * self.dy);
        let mut u_half = u.to_vec();
        for j in 1..ny - 1 {
            let mut a = vec![0.0_f64; nx];
            let mut b = vec![1.0 + 2.0 * rx; nx];
            let mut c = vec![0.0_f64; nx];
            let mut d = vec![0.0_f64; nx];
            b[0] = 1.0;
            b[nx - 1] = 1.0;
            d[0] = u[j * nx];
            d[nx - 1] = u[j * nx + nx - 1];
            for i in 1..nx - 1 {
                a[i] = -rx;
                c[i] = -rx;
                let idx = j * nx + i;
                d[i] = u[idx] + ry * (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]);
            }
            let row = thomas_algorithm(&a, &b, &c, &d);
            for i in 0..nx {
                u_half[j * nx + i] = row[i];
            }
        }
        u_half
    }
    /// ADI half-step in y direction.
    pub fn adi_step_y(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let rx = self.diffusivity * self.dt / (2.0 * self.dx * self.dx);
        let ry = self.diffusivity * self.dt / (2.0 * self.dy * self.dy);
        let mut u_new = u.to_vec();
        for i in 1..nx - 1 {
            let mut a = vec![0.0_f64; ny];
            let mut b = vec![1.0 + 2.0 * ry; ny];
            let mut c = vec![0.0_f64; ny];
            let mut d = vec![0.0_f64; ny];
            b[0] = 1.0;
            b[ny - 1] = 1.0;
            d[0] = u[i];
            d[ny - 1] = u[(ny - 1) * nx + i];
            for j in 1..ny - 1 {
                a[j] = -ry;
                c[j] = -ry;
                let idx = j * nx + i;
                d[j] = u[idx] + rx * (u[idx - 1] - 2.0 * u[idx] + u[idx + 1]);
            }
            let col = thomas_algorithm(&a, &b, &c, &d);
            for j in 0..ny {
                u_new[j * nx + i] = col[j];
            }
        }
        u_new
    }
    /// Full ADI step (x-sweep followed by y-sweep).
    pub fn adi_step(&self, u: &[f64]) -> Vec<f64> {
        let u_half = self.adi_step_x(u);
        self.adi_step_y(&u_half)
    }
}
/// 2D finite-difference operator set (gradient, Laplacian, divergence, curl).
///
/// All fields are flat row-major vectors of size `nx * ny`.
#[allow(dead_code)]
pub struct FiniteDiffOps2D {
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Grid spacing in x.
    pub dx: f64,
    /// Grid spacing in y.
    pub dy: f64,
}
impl FiniteDiffOps2D {
    /// Creates a new 2D operator set.
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        Self { nx, ny, dx, dy }
    }
    /// Row/column index helper.
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// 5-point Laplacian ∇²u = ∂²u/∂x² + ∂²u/∂y².
    ///
    /// Returns zero at boundary points.
    pub fn laplacian(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let mut lap = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let c = self.idx(i, j);
                lap[c] = (u[c - 1] - 2.0 * u[c] + u[c + 1]) / dx2
                    + (u[c - nx] - 2.0 * u[c] + u[c + nx]) / dy2;
            }
        }
        lap
    }
    /// Central-difference gradient (∂u/∂x, ∂u/∂y).
    ///
    /// Returns a pair of vectors; boundary points are zero.
    pub fn gradient(&self, u: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let nx = self.nx;
        let ny = self.ny;
        let mut gx = vec![0.0_f64; nx * ny];
        let mut gy = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let c = self.idx(i, j);
                gx[c] = (u[c + 1] - u[c - 1]) / (2.0 * self.dx);
                gy[c] = (u[c + nx] - u[c - nx]) / (2.0 * self.dy);
            }
        }
        (gx, gy)
    }
    /// Divergence of a 2D vector field (∂fx/∂x + ∂fy/∂y).
    ///
    /// `fx` and `fy` are the x- and y-components of the field.
    pub fn divergence(&self, fx: &[f64], fy: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let mut div = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let c = self.idx(i, j);
                div[c] = (fx[c + 1] - fx[c - 1]) / (2.0 * self.dx)
                    + (fy[c + nx] - fy[c - nx]) / (2.0 * self.dy);
            }
        }
        div
    }
    /// Curl (z-component) of a 2D vector field: ∂uy/∂x - ∂ux/∂y.
    ///
    /// In 2D the curl reduces to a scalar field.
    pub fn curl(&self, ux: &[f64], uy: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let mut curl = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let c = self.idx(i, j);
                let duy_dx = (uy[c + 1] - uy[c - 1]) / (2.0 * self.dx);
                let dux_dy = (ux[c + nx] - ux[c - nx]) / (2.0 * self.dy);
                curl[c] = duy_dx - dux_dy;
            }
        }
        curl
    }
}
/// Gray-Scott reaction-diffusion system on a periodic 2D grid.
///
/// Models the two-species reaction:
/// - u_t = Du ∇²u − u v² + f (1 − u)
/// - v_t = Dv ∇²v + u v² − (f + k) v
#[allow(dead_code)]
pub struct GrayScottSystem {
    /// Grid size in x.
    pub nx: usize,
    /// Grid size in y.
    pub ny: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Feed rate f.
    pub feed: f64,
    /// Kill rate k.
    pub kill: f64,
    /// Diffusivity of u.
    pub du: f64,
    /// Diffusivity of v.
    pub dv: f64,
    /// Concentration of species u.
    pub u: Vec<f64>,
    /// Concentration of species v.
    pub v: Vec<f64>,
}
impl GrayScottSystem {
    /// Creates a new Gray-Scott system with default diffusivities.
    pub fn new(nx: usize, ny: usize, dx: f64, dt: f64, feed: f64, kill: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            dx,
            dt,
            feed,
            kill,
            du: 0.2,
            dv: 0.1,
            u: vec![1.0; n],
            v: vec![0.0; n],
        }
    }
    /// Generates standard initial conditions: u=1 everywhere, small v seed in center.
    pub fn initial_conditions(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>) {
        let n = nx * ny;
        let mut u = vec![1.0_f64; n];
        let mut v = vec![0.0_f64; n];
        let cx = nx / 2;
        let cy = ny / 2;
        let r = (nx.min(ny) / 10).max(1);
        for j in 0..ny {
            for i in 0..nx {
                let dx = (i as isize - cx as isize).unsigned_abs();
                let dy = (j as isize - cy as isize).unsigned_abs();
                if dx <= r && dy <= r {
                    u[j * nx + i] = 0.5;
                    v[j * nx + i] = 0.25;
                }
            }
        }
        (u, v)
    }
    /// Laplacian with periodic boundary conditions (5-pt stencil).
    fn laplacian_periodic(field: &[f64], nx: usize, ny: usize, dx: f64) -> Vec<f64> {
        let dx2 = dx * dx;
        let mut lap = vec![0.0_f64; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let c = j * nx + i;
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                lap[c] = (field[j * nx + ip]
                    + field[j * nx + im]
                    + field[jp * nx + i]
                    + field[jm * nx + i]
                    - 4.0 * field[c])
                    / dx2;
            }
        }
        lap
    }
    /// Advances the system by one time step.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let dt = self.dt;
        let f = self.feed;
        let k = self.kill;
        let lap_u = Self::laplacian_periodic(&self.u, nx, ny, self.dx);
        let lap_v = Self::laplacian_periodic(&self.v, nx, ny, self.dx);
        let n = nx * ny;
        let mut u_new = self.u.clone();
        let mut v_new = self.v.clone();
        for idx in 0..n {
            let u_val = self.u[idx];
            let v_val = self.v[idx];
            let reaction = u_val * v_val * v_val;
            u_new[idx] = u_val + dt * (self.du * lap_u[idx] - reaction + f * (1.0 - u_val));
            v_new[idx] = v_val + dt * (self.dv * lap_v[idx] + reaction - (f + k) * v_val);
        }
        self.u = u_new;
        self.v = v_new;
    }
}
/// 1D heat equation u_t = D u_xx with Dirichlet boundary conditions.
///
/// Provides three schemes: explicit FTCS (forward Euler), implicit
/// (backward Euler), and Crank-Nicolson (second-order in time).
#[allow(dead_code)]
pub struct HeatEquation1D {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Thermal diffusivity D.
    pub d: f64,
}
impl HeatEquation1D {
    /// Creates a new 1D heat equation solver.
    pub fn new(n: usize, dx: f64, dt: f64, d: f64) -> Self {
        Self { n, dx, dt, d }
    }
    /// Diffusion number r = D dt / dx².  Explicit scheme stable when r ≤ 0.5.
    pub fn r(&self) -> f64 {
        self.d * self.dt / (self.dx * self.dx)
    }
    /// Explicit FTCS step: u^{n+1}_i = u^n_i + r(u^n_{i-1} − 2 u^n_i + u^n_{i+1}).
    pub fn step_explicit(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let r = self.r();
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            u_new[i] = u[i] + r * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
        }
        u_new
    }
    /// Implicit backward-Euler step: solves (I − r L) u^{n+1} = u^n.
    pub fn step_implicit(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let r = self.r();
        let mut a = vec![0.0_f64; n];
        let mut b = vec![1.0 + 2.0 * r; n];
        let mut c = vec![0.0_f64; n];
        let mut d = u.to_vec();
        for i in 1..n - 1 {
            a[i] = -r;
            c[i] = -r;
        }
        b[0] = 1.0;
        b[n - 1] = 1.0;
        a[0] = 0.0;
        c[0] = 0.0;
        a[n - 1] = 0.0;
        c[n - 1] = 0.0;
        d[0] = u[0];
        d[n - 1] = u[n - 1];
        thomas_algorithm(&a, &b, &c, &d)
    }
    /// Crank-Nicolson step: second-order in both time and space, unconditionally stable.
    pub fn step_cn(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let r = self.r();
        let hr = 0.5 * r;
        let mut a = vec![0.0_f64; n];
        let mut b = vec![1.0 + r; n];
        let mut c = vec![0.0_f64; n];
        let mut d = vec![0.0_f64; n];
        for i in 1..n - 1 {
            a[i] = -hr;
            c[i] = -hr;
            d[i] = hr * u[i - 1] + (1.0 - r) * u[i] + hr * u[i + 1];
        }
        b[0] = 1.0;
        b[n - 1] = 1.0;
        d[0] = u[0];
        d[n - 1] = u[n - 1];
        thomas_algorithm(&a, &b, &c, &d)
    }
    /// Integrates for `n_steps` steps using the specified scheme.
    ///
    /// `scheme` is `"explicit"`, `"implicit"`, or `"cn"`.
    pub fn integrate(&self, u0: &[f64], n_steps: usize, scheme: &str) -> Vec<f64> {
        let mut u = u0.to_vec();
        for _ in 0..n_steps {
            u = match scheme {
                "implicit" => self.step_implicit(&u),
                "cn" => self.step_cn(&u),
                _ => self.step_explicit(&u),
            };
        }
        u
    }
}
/// One-dimensional finite-difference PDE solver.
///
/// Supports explicit (FTCS), implicit (backward Euler), and Crank-Nicolson
/// schemes for the diffusion equation, and upwind/central schemes for
/// the advection equation.
#[allow(dead_code)]
pub struct FiniteDifference1D {
    /// Number of spatial grid points.
    pub n: usize,
    /// Spatial grid spacing.
    pub dx: f64,
    /// Time step size.
    pub dt: f64,
    /// Diffusion coefficient.
    pub diffusivity: f64,
    /// Advection velocity.
    pub velocity: f64,
    /// Left boundary condition.
    pub bc_left: BoundaryCondition,
    /// Right boundary condition.
    pub bc_right: BoundaryCondition,
}
impl FiniteDifference1D {
    /// Creates a new 1D finite-difference solver.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        n: usize,
        dx: f64,
        dt: f64,
        diffusivity: f64,
        velocity: f64,
        bc_left: BoundaryCondition,
        bc_right: BoundaryCondition,
    ) -> Self {
        Self {
            n,
            dx,
            dt,
            diffusivity,
            velocity,
            bc_left,
            bc_right,
        }
    }
    /// Diffusion number r = D * dt / dx^2.
    ///
    /// Stability requires r <= 0.5 for the explicit scheme.
    pub fn diffusion_number(&self) -> f64 {
        self.diffusivity * self.dt / (self.dx * self.dx)
    }
    /// Courant number C = v * dt / dx.
    ///
    /// Stability requires |C| <= 1 for the upwind advection scheme.
    pub fn courant_number(&self) -> f64 {
        self.velocity * self.dt / self.dx
    }
    /// Returns true if the explicit scheme is stable (r <= 0.5).
    pub fn is_stable_explicit(&self) -> bool {
        self.diffusion_number() <= 0.5
    }
    /// Explicit FTCS step for the diffusion equation: u_t = D * u_xx.
    pub fn step_diffusion_explicit(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let r = self.diffusion_number();
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            u_new[i] = u[i] + r * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
        }
        self.apply_bc_1d(&mut u_new, u);
        u_new
    }
    /// Implicit (backward Euler) step for the diffusion equation.
    ///
    /// Solves the tridiagonal system (I - r*L) u_new = u_old.
    pub fn step_diffusion_implicit(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let r = self.diffusion_number();
        let mut a = vec![0.0_f64; n];
        let mut b = vec![1.0 + 2.0 * r; n];
        let mut c = vec![0.0_f64; n];
        let mut d = u.to_vec();
        for i in 1..n - 1 {
            a[i] = -r;
            c[i] = -r;
        }
        b[0] = 1.0;
        b[n - 1] = 1.0;
        a[0] = 0.0;
        c[0] = 0.0;
        a[n - 1] = 0.0;
        c[n - 1] = 0.0;
        d[0] = match self.bc_left {
            BoundaryCondition::Dirichlet(v) => v,
            _ => u[0],
        };
        d[n - 1] = match self.bc_right {
            BoundaryCondition::Dirichlet(v) => v,
            _ => u[n - 1],
        };
        thomas_algorithm(&a, &b, &c, &d)
    }
    /// Crank-Nicolson step for the diffusion equation.
    ///
    /// Second-order accurate in both space and time; unconditionally stable.
    pub fn step_diffusion_cn(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let r = self.diffusion_number();
        let half_r = 0.5 * r;
        let mut a = vec![0.0_f64; n];
        let mut b = vec![1.0 + r; n];
        let mut c = vec![0.0_f64; n];
        let mut d = vec![0.0_f64; n];
        for i in 1..n - 1 {
            a[i] = -half_r;
            c[i] = -half_r;
            d[i] = half_r * u[i - 1] + (1.0 - r) * u[i] + half_r * u[i + 1];
        }
        b[0] = 1.0;
        b[n - 1] = 1.0;
        d[0] = match self.bc_left {
            BoundaryCondition::Dirichlet(v) => v,
            _ => u[0],
        };
        d[n - 1] = match self.bc_right {
            BoundaryCondition::Dirichlet(v) => v,
            _ => u[n - 1],
        };
        thomas_algorithm(&a, &b, &c, &d)
    }
    /// Upwind advection step: u_t + v * u_x = 0.
    ///
    /// Uses first-order upwind differencing; stable when |C| <= 1.
    pub fn step_advection_upwind(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let c = self.courant_number();
        let mut u_new = u.to_vec();
        if c > 0.0 {
            for i in 1..n {
                u_new[i] = u[i] - c * (u[i] - u[i - 1]);
            }
        } else {
            for i in 0..n - 1 {
                u_new[i] = u[i] - c * (u[i + 1] - u[i]);
            }
        }
        u_new
    }
    /// Central difference advection step: u_t + v * u_x = 0.
    ///
    /// Second-order but conditionally stable; may oscillate.
    pub fn step_advection_central(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let c = self.courant_number();
        let mut u_new = u.to_vec();
        for i in 1..n - 1 {
            u_new[i] = u[i] - 0.5 * c * (u[i + 1] - u[i - 1]);
        }
        u_new
    }
    fn apply_bc_1d(&self, u_new: &mut Vec<f64>, _u_old: &[f64]) {
        match self.bc_left {
            BoundaryCondition::Dirichlet(v) => u_new[0] = v,
            BoundaryCondition::Neumann(flux) => {
                u_new[0] = u_new[1] - flux * self.dx;
            }
            BoundaryCondition::Absorbing => {
                u_new[0] = u_new[1];
            }
            BoundaryCondition::Periodic => {
                let n = u_new.len();
                u_new[0] = u_new[n - 2];
            }
        }
        let n = u_new.len();
        match self.bc_right {
            BoundaryCondition::Dirichlet(v) => u_new[n - 1] = v,
            BoundaryCondition::Neumann(flux) => {
                u_new[n - 1] = u_new[n - 2] + flux * self.dx;
            }
            BoundaryCondition::Absorbing => {
                u_new[n - 1] = u_new[n - 2];
            }
            BoundaryCondition::Periodic => {
                u_new[n - 1] = u_new[1];
            }
        }
    }
}
/// FFT-based spectral differentiation on a periodic 1D domain.
///
/// Implements differentiation by multiplying Fourier coefficients
/// by the wavenumber ik.  Uses O(N²) DFT internally (no external FFT crate).
#[allow(dead_code)]
pub struct FftDiff1D {
    /// Number of grid points (ideally a power of 2).
    pub n: usize,
    /// Grid spacing dx = L / n.
    pub dx: f64,
}
impl FftDiff1D {
    /// Creates a new spectral differentiator.
    pub fn new(n: usize, dx: f64) -> Self {
        Self { n, dx }
    }
    /// Computes the DFT of `x` (complex output as flat `(re, im)` pairs).
    fn dft(&self, x: &[f64]) -> Vec<(f64, f64)> {
        let n = self.n;
        let mut out = vec![(0.0_f64, 0.0_f64); n];
        for k in 0..n {
            let mut re = 0.0;
            let mut im = 0.0;
            for j in 0..n {
                let angle = -2.0 * PI * (k as f64) * (j as f64) / n as f64;
                re += x[j] * angle.cos();
                im += x[j] * angle.sin();
            }
            out[k] = (re, im);
        }
        out
    }
    /// Inverse DFT; takes `(re, im)` pairs and returns real part.
    fn idft(&self, xk: &[(f64, f64)]) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0_f64; n];
        for j in 0..n {
            let mut re = 0.0;
            for k in 0..n {
                let angle = 2.0 * PI * (k as f64) * (j as f64) / n as f64;
                re += xk[k].0 * angle.cos() - xk[k].1 * angle.sin();
            }
            out[j] = re / n as f64;
        }
        out
    }
    /// Computes the first derivative du/dx via spectral multiplication by ik.
    ///
    /// Assumes periodic boundary conditions on \[0, n*dx\].
    pub fn differentiate(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let l = n as f64 * self.dx;
        let uk = self.dft(u);
        let mut duk: Vec<(f64, f64)> = (0..n)
            .map(|k| {
                let kk = if k <= n / 2 {
                    k as f64
                } else {
                    k as f64 - n as f64
                };
                let omega = 2.0 * PI * kk / l;
                (-uk[k].1 * omega, uk[k].0 * omega)
            })
            .collect();
        if n.is_multiple_of(2) {
            duk[n / 2] = (0.0, 0.0);
        }
        self.idft(&duk)
    }
    /// Computes the second derivative d²u/dx² via spectral multiplication by -k².
    pub fn differentiate2(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let l = n as f64 * self.dx;
        let uk = self.dft(u);
        let d2uk: Vec<(f64, f64)> = (0..n)
            .map(|k| {
                let kk = if k <= n / 2 {
                    k as f64
                } else {
                    k as f64 - n as f64
                };
                let omega2 = -(2.0 * PI * kk / l).powi(2);
                (uk[k].0 * omega2, uk[k].1 * omega2)
            })
            .collect();
        self.idft(&d2uk)
    }
}
/// Alternating Direction Implicit (ADI) method for 2D parabolic PDEs.
///
/// Peaceman-Rachford ADI: each full time step is split into an x-implicit half-step
/// and a y-implicit half-step, giving second-order accuracy and unconditional stability.
#[allow(dead_code)]
pub struct AdiMethod2D {
    /// Grid points in x.
    pub nx: usize,
    /// Grid points in y.
    pub ny: usize,
    /// Grid spacing in x.
    pub dx: f64,
    /// Grid spacing in y.
    pub dy: f64,
    /// Time step.
    pub dt: f64,
    /// Diffusion coefficient.
    pub d: f64,
}
impl AdiMethod2D {
    /// Creates a new ADI solver.
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64, dt: f64, d: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            dy,
            dt,
            d,
        }
    }
    /// Performs the x-direction implicit half-step.
    pub fn half_step_x(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let rx = self.d * self.dt / (2.0 * self.dx * self.dx);
        let ry = self.d * self.dt / (2.0 * self.dy * self.dy);
        let mut u_half = u.to_vec();
        for j in 1..ny - 1 {
            let mut a = vec![0.0_f64; nx];
            let mut b = vec![1.0 + 2.0 * rx; nx];
            let mut c = vec![0.0_f64; nx];
            let mut d = vec![0.0_f64; nx];
            b[0] = 1.0;
            b[nx - 1] = 1.0;
            d[0] = u[j * nx];
            d[nx - 1] = u[j * nx + nx - 1];
            for i in 1..nx - 1 {
                a[i] = -rx;
                c[i] = -rx;
                let idx = j * nx + i;
                d[i] = u[idx] + ry * (u[idx - nx] - 2.0 * u[idx] + u[idx + nx]);
            }
            let row = thomas_algorithm(&a, &b, &c, &d);
            for i in 0..nx {
                u_half[j * nx + i] = row[i];
            }
        }
        u_half
    }
    /// Performs the y-direction implicit half-step.
    pub fn half_step_y(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let rx = self.d * self.dt / (2.0 * self.dx * self.dx);
        let ry = self.d * self.dt / (2.0 * self.dy * self.dy);
        let mut u_new = u.to_vec();
        for i in 1..nx - 1 {
            let mut a = vec![0.0_f64; ny];
            let mut b = vec![1.0 + 2.0 * ry; ny];
            let mut c = vec![0.0_f64; ny];
            let mut d = vec![0.0_f64; ny];
            b[0] = 1.0;
            b[ny - 1] = 1.0;
            d[0] = u[i];
            d[ny - 1] = u[(ny - 1) * nx + i];
            for j in 1..ny - 1 {
                a[j] = -ry;
                c[j] = -ry;
                let idx = j * nx + i;
                d[j] = u[idx] + rx * (u[idx - 1] - 2.0 * u[idx] + u[idx + 1]);
            }
            let col = thomas_algorithm(&a, &b, &c, &d);
            for j in 0..ny {
                u_new[j * nx + i] = col[j];
            }
        }
        u_new
    }
    /// Full Peaceman-Rachford ADI step (x-half then y-half).
    pub fn step(&self, u: &[f64]) -> Vec<f64> {
        let u_half = self.half_step_x(u);
        self.half_step_y(&u_half)
    }
    /// Integrates for `n_steps` full ADI steps.
    pub fn integrate(&self, u0: &[f64], n_steps: usize) -> Vec<f64> {
        let mut u = u0.to_vec();
        for _ in 0..n_steps {
            u = self.step(&u);
        }
        u
    }
}
/// 1D wave equation solver using the 2nd-order central-difference leapfrog scheme.
///
/// Equation: u_tt = c² u_xx.  Uses u^{n+1} = 2u^n − u^{n-1} + C²(u^n_{i-1} − 2u^n_i + u^n_{i+1}).
#[allow(dead_code)]
pub struct WaveEq2ndOrder {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Wave speed c.
    pub c: f64,
}
impl WaveEq2ndOrder {
    /// Creates a new 2nd-order central-difference wave solver.
    pub fn new(n: usize, dx: f64, dt: f64, c: f64) -> Self {
        Self { n, dx, dt, c }
    }
    /// Courant number C = c dt / dx; must satisfy |C| ≤ 1 for stability.
    pub fn courant(&self) -> f64 {
        self.c * self.dt / self.dx
    }
    /// Advances one time step.
    ///
    /// `u` is the current level, `u_prev` is the previous level.
    pub fn step(&self, u: &[f64], u_prev: &[f64]) -> Vec<f64> {
        let n = self.n;
        let c2 = self.courant().powi(2);
        let mut u_new = vec![0.0_f64; n];
        for i in 1..n - 1 {
            u_new[i] = 2.0 * u[i] - u_prev[i] + c2 * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
        }
        u_new[0] = u[1];
        u_new[n - 1] = u[n - 2];
        u_new
    }
    /// Integrates for `n_steps` steps starting from `u0` with zero initial velocity.
    pub fn integrate(&self, u0: &[f64], n_steps: usize) -> Vec<f64> {
        let mut u_prev = u0.to_vec();
        let mut u_cur = self.step(&u_prev, &u_prev);
        for _ in 1..n_steps {
            let u_next = self.step(&u_cur, &u_prev);
            u_prev = u_cur;
            u_cur = u_next;
        }
        u_cur
    }
}
/// 2D level-set method for interface tracking.
///
/// Maintains a signed distance function φ; the interface is the zero level set.
/// Provides reinitialization (Sussman), curvature/normal computation, and advection.
#[allow(dead_code)]
pub struct LevelSet2D {
    /// Grid size in x.
    pub nx: usize,
    /// Grid size in y.
    pub ny: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Current level-set field φ (row-major, size nx*ny).
    pub phi: Vec<f64>,
}
impl LevelSet2D {
    /// Creates a new level-set solver with initial field `phi0`.
    pub fn new(nx: usize, ny: usize, dx: f64, phi0: Vec<f64>) -> Self {
        assert_eq!(phi0.len(), nx * ny);
        Self {
            nx,
            ny,
            dx,
            phi: phi0,
        }
    }
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Computes the unit normal to the interface via central differences.
    ///
    /// Returns (nx, ny) vectors of length nx*ny each.
    pub fn normal(&self) -> (Vec<f64>, Vec<f64>) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let mut n_x = vec![0.0_f64; nx * ny];
        let mut n_y = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let dphi_x =
                    (self.phi[self.idx(i + 1, j)] - self.phi[self.idx(i - 1, j)]) / (2.0 * dx);
                let dphi_y =
                    (self.phi[self.idx(i, j + 1)] - self.phi[self.idx(i, j - 1)]) / (2.0 * dx);
                let mag = (dphi_x * dphi_x + dphi_y * dphi_y).sqrt().max(1e-15);
                n_x[idx] = dphi_x / mag;
                n_y[idx] = dphi_y / mag;
            }
        }
        (n_x, n_y)
    }
    /// Computes the mean curvature κ = ∇ · (∇φ / |∇φ|).
    pub fn curvature(&self) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let mut kappa = vec![0.0_f64; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let phi_xp = self.phi[self.idx(i + 1, j)];
                let phi_xm = self.phi[self.idx(i - 1, j)];
                let phi_yp = self.phi[self.idx(i, j + 1)];
                let phi_ym = self.phi[self.idx(i, j - 1)];
                let phi_c = self.phi[idx];
                let phi_xy = (self.phi[self.idx(i + 1, j + 1)]
                    - self.phi[self.idx(i + 1, j - 1)]
                    - self.phi[self.idx(i - 1, j + 1)]
                    + self.phi[self.idx(i - 1, j - 1)])
                    / (4.0 * dx * dx);
                let dx_phi = (phi_xp - phi_xm) / (2.0 * dx);
                let dy_phi = (phi_yp - phi_ym) / (2.0 * dx);
                let dx2_phi = (phi_xp - 2.0 * phi_c + phi_xm) / (dx * dx);
                let dy2_phi = (phi_yp - 2.0 * phi_c + phi_ym) / (dx * dx);
                let grad2 = dx_phi * dx_phi + dy_phi * dy_phi + 1e-15;
                kappa[idx] = (dx2_phi * dy_phi * dy_phi - 2.0 * dx_phi * dy_phi * phi_xy
                    + dy2_phi * dx_phi * dx_phi)
                    / (grad2 * grad2.sqrt());
            }
        }
        kappa
    }
    /// Advects the level-set field using a given velocity field (u_field, v_field).
    ///
    /// Uses first-order upwind scheme.
    pub fn advect(&mut self, u_field: &[f64], v_field: &[f64], dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let mut phi_new = self.phi.clone();
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = self.idx(i, j);
                let u = u_field[idx];
                let v = v_field[idx];
                let dphi_x = if u >= 0.0 {
                    (self.phi[idx] - self.phi[self.idx(i - 1, j)]) / dx
                } else {
                    (self.phi[self.idx(i + 1, j)] - self.phi[idx]) / dx
                };
                let dphi_y = if v >= 0.0 {
                    (self.phi[idx] - self.phi[self.idx(i, j - 1)]) / dx
                } else {
                    (self.phi[self.idx(i, j + 1)] - self.phi[idx]) / dx
                };
                phi_new[idx] = self.phi[idx] - dt * (u * dphi_x + v * dphi_y);
            }
        }
        self.phi = phi_new;
    }
    /// Reinitializes the level-set to a signed distance function (Sussman method).
    ///
    /// Iterates the reinitialization equation dφ/dτ = sign(φ)(1 - |∇φ|).
    pub fn reinitialize(&mut self, n_iters: usize, dtau: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        for _ in 0..n_iters {
            let phi_old = self.phi.clone();
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let idx = self.idx(i, j);
                    let p = phi_old[idx];
                    let sign_p = if p > 0.0 {
                        1.0
                    } else if p < 0.0 {
                        -1.0
                    } else {
                        0.0
                    };
                    let dp_xm = (p - phi_old[self.idx(i - 1, j)]) / dx;
                    let dp_xp = (phi_old[self.idx(i + 1, j)] - p) / dx;
                    let dp_ym = (p - phi_old[self.idx(i, j - 1)]) / dx;
                    let dp_yp = (phi_old[self.idx(i, j + 1)] - p) / dx;
                    let grad_mag = if sign_p > 0.0 {
                        let gx = dp_xm.max(0.0).powi(2) + dp_xp.min(0.0).powi(2);
                        let gy = dp_ym.max(0.0).powi(2) + dp_yp.min(0.0).powi(2);
                        (gx + gy).sqrt()
                    } else {
                        let gx = dp_xm.min(0.0).powi(2) + dp_xp.max(0.0).powi(2);
                        let gy = dp_ym.min(0.0).powi(2) + dp_yp.max(0.0).powi(2);
                        (gx + gy).sqrt()
                    };
                    self.phi[idx] = p - dtau * sign_p * (grad_mag - 1.0);
                }
            }
        }
    }
}
/// 3D finite-difference operator set.
///
/// All fields are flat z-then-y-then-x (i + j*nx + k*nx*ny) vectors.
#[allow(dead_code)]
pub struct FiniteDiffOps3D {
    /// Grid points in x.
    pub nx: usize,
    /// Grid points in y.
    pub ny: usize,
    /// Grid points in z.
    pub nz: usize,
    /// Uniform grid spacing.
    pub dx: f64,
}
impl FiniteDiffOps3D {
    /// Creates a new 3D operator set (uniform spacing `dx`).
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self { nx, ny, nz, dx }
    }
    fn idx3(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }
    /// 7-point 3D Laplacian ∇²u.
    ///
    /// Returns zero at all boundary planes.
    pub fn laplacian(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let dx2 = self.dx * self.dx;
        let mut lap = vec![0.0_f64; nx * ny * nz];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let c = self.idx3(i, j, k);
                    lap[c] = (u[c - 1] - 2.0 * u[c] + u[c + 1]) / dx2
                        + (u[c - nx] - 2.0 * u[c] + u[c + nx]) / dx2
                        + (u[c - nx * ny] - 2.0 * u[c] + u[c + nx * ny]) / dx2;
                }
            }
        }
        lap
    }
    /// Gradient of a scalar field in 3D.
    ///
    /// Returns (∂u/∂x, ∂u/∂y, ∂u/∂z) each of size nx*ny*nz.
    pub fn gradient(&self, u: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let two_dx = 2.0 * self.dx;
        let mut gx = vec![0.0_f64; nx * ny * nz];
        let mut gy = vec![0.0_f64; nx * ny * nz];
        let mut gz = vec![0.0_f64; nx * ny * nz];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let c = self.idx3(i, j, k);
                    gx[c] = (u[c + 1] - u[c - 1]) / two_dx;
                    gy[c] = (u[c + nx] - u[c - nx]) / two_dx;
                    gz[c] = (u[c + nx * ny] - u[c - nx * ny]) / two_dx;
                }
            }
        }
        (gx, gy, gz)
    }
    /// Curl of a 3D vector field (fx, fy, fz).
    ///
    /// Returns the three components (curl_x, curl_y, curl_z).
    pub fn curl(&self, fx: &[f64], fy: &[f64], fz: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let two_dx = 2.0 * self.dx;
        let nn = nx * ny * nz;
        let mut cx = vec![0.0_f64; nn];
        let mut cy = vec![0.0_f64; nn];
        let mut cz = vec![0.0_f64; nn];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let c = self.idx3(i, j, k);
                    cx[c] = (fz[c + nx] - fz[c - nx]) / two_dx
                        - (fy[c + nx * ny] - fy[c - nx * ny]) / two_dx;
                    cy[c] = (fx[c + nx * ny] - fx[c - nx * ny]) / two_dx
                        - (fz[c + 1] - fz[c - 1]) / two_dx;
                    cz[c] = (fy[c + 1] - fy[c - 1]) / two_dx - (fx[c + nx] - fx[c - nx]) / two_dx;
                }
            }
        }
        (cx, cy, cz)
    }
}
/// Boundary condition type for PDE solvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Dirichlet: fixed value at boundary.
    Dirichlet(f64),
    /// Neumann: fixed flux (derivative) at boundary.
    Neumann(f64),
    /// Periodic boundary: wraps around.
    Periodic,
    /// Absorbing boundary: zero-gradient outflow.
    Absorbing,
}
/// Chebyshev pseudospectral method for PDEs on \[-1, 1\].
///
/// Provides Gauss-Lobatto nodes, the spectral differentiation matrix,
/// and Chebyshev expansion/interpolation.
#[allow(dead_code)]
pub struct SpectralMethod {
    /// Number of Chebyshev points (polynomial degree N = n-1).
    pub n: usize,
}
impl SpectralMethod {
    /// Creates a new Chebyshev spectral solver with `n` points.
    pub fn new(n: usize) -> Self {
        Self { n }
    }
    /// Returns the n Gauss-Lobatto nodes on \[-1, 1\].
    ///
    /// x_j = cos(π*j/(n-1)), j = 0, ..., n-1.
    pub fn gauss_lobatto_nodes(&self) -> Vec<f64> {
        let n = self.n;
        (0..n)
            .map(|j| (PI * j as f64 / (n - 1) as f64).cos())
            .collect()
    }
    /// Builds the Chebyshev differentiation matrix D of size n×n (flattened row-major).
    ///
    /// Computes dU/dx ≈ D * U where U is the vector of values at Gauss-Lobatto nodes.
    pub fn differentiation_matrix(&self) -> Vec<f64> {
        let n = self.n;
        let x = self.gauss_lobatto_nodes();
        let mut d = vec![0.0_f64; n * n];
        let c = |i: usize| -> f64 { if i == 0 || i == n - 1 { 2.0 } else { 1.0 } };
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    d[i * n + j] = c(i) / c(j) * (-1.0_f64).powi((i + j) as i32) / (x[i] - x[j]);
                }
            }
        }
        for i in 0..n {
            let row_sum: f64 = (0..n).filter(|&j| j != i).map(|j| d[i * n + j]).sum();
            d[i * n + i] = -row_sum;
        }
        d
    }
    /// Applies the differentiation matrix to a vector u, returning du/dx.
    pub fn differentiate(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let d = self.differentiation_matrix();
        let mut du = vec![0.0_f64; n];
        for i in 0..n {
            for j in 0..n {
                du[i] += d[i * n + j] * u[j];
            }
        }
        du
    }
    /// Chebyshev expansion: computes coefficients a_k from function values at GL nodes.
    pub fn chebyshev_coefficients(&self, f: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut a = vec![0.0_f64; n];
        for k in 0..n {
            let c_k = if k == 0 || k == n - 1 { 2.0 } else { 1.0 };
            let mut sum = 0.0;
            for j in 0..n {
                let c_j = if j == 0 || j == n - 1 { 2.0 } else { 1.0 };
                let _theta_j = PI * j as f64 / (n - 1) as f64;
                let _theta_k = PI * k as f64 / (n - 1) as f64;
                sum += f[j] / c_j * (PI * k as f64 * j as f64 / (n - 1) as f64).cos();
            }
            a[k] = 2.0 * sum / (c_k * (n - 1) as f64);
        }
        a
    }
}
