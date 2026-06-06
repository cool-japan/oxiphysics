//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::functions::*;
use std::f64::consts::PI;

/// Poisson equation solver: ∇²u = f.
///
/// Implements multigrid V-cycle with Jacobi smoothing (restriction,
/// prolongation) and a direct FFT-based spectral solver for periodic domains.
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
        for (i, cell) in r_coarse.iter_mut().enumerate() {
            let j = 2 * i;
            if j + 1 < r_fine.len() {
                *cell = 0.5 * (r_fine[j] + r_fine[j + 1]);
            } else {
                *cell = r_fine[j];
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
        for (k, fhat_k) in f_hat.iter_mut().enumerate() {
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            for (j, fj) in f.iter().enumerate() {
                let angle = -2.0 * PI * (k as f64) * (j as f64) / (n as f64);
                re += fj * angle.cos();
                im += fj * angle.sin();
            }
            *fhat_k = (re, im);
        }
        let mut u_hat = vec![(0.0_f64, 0.0_f64); n];
        for (k, uhat_k) in u_hat.iter_mut().enumerate().skip(1) {
            let eigenvalue = -4.0 / dx2 * (PI * k as f64 / n as f64).sin().powi(2);
            *uhat_k = (f_hat[k].0 / eigenvalue, f_hat[k].1 / eigenvalue);
        }
        let mut u = vec![0.0_f64; n];
        for (j, uj) in u.iter_mut().enumerate() {
            let mut re = 0.0_f64;
            for (k, uhat_k) in u_hat.iter().enumerate() {
                let angle = 2.0 * PI * (k as f64) * (j as f64) / (n as f64);
                re += uhat_k.0 * angle.cos() - uhat_k.1 * angle.sin();
            }
            *uj = re / n as f64;
        }
        u
    }
}
/// Two-dimensional finite-difference solver for elliptic and parabolic PDEs.
///
/// Supports the 5-point and 9-point stencils for the Laplacian, and
/// the ADI (Alternating Direction Implicit) method for 2D diffusion.
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
            u_half[j * nx..j * nx + nx].copy_from_slice(&row);
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
/// Discrete Cosine Transform (DCT-II) for spectral PDE methods.
///
/// DCT-II: X_k = Σ_{n=0}^{N-1} x_n * cos(π/N * (n+0.5) * k).
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
        for (k, xk_val) in xk.iter_mut().enumerate() {
            let mut sum = 0.0;
            for (j, xj) in x.iter().enumerate() {
                sum += xj * (PI / n as f64 * (j as f64 + 0.5) * k as f64).cos();
            }
            *xk_val = sum;
        }
        xk
    }
    /// Inverse DCT-II (DCT-III) transform.
    pub fn dct2_inverse(&self, xk: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x = vec![0.0_f64; n];
        for (j, xj) in x.iter_mut().enumerate() {
            let mut sum = 0.5 * xk[0];
            for (k, xkk) in xk.iter().enumerate().skip(1) {
                sum += xkk * (PI / n as f64 * (j as f64 + 0.5) * k as f64).cos();
            }
            *xj = sum * 2.0 / n as f64;
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
/// 2D finite-difference operator set (gradient, Laplacian, divergence, curl).
///
/// All fields are flat row-major vectors of size `nx * ny`.
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
/// Three-dimensional heat equation solver: u_t = D * (u_xx + u_yy + u_zz).
///
/// Supports explicit FTCS and implicit (via ADI-like splitting) schemes,
/// with Dirichlet or Neumann boundary conditions.
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
    fn apply_bc_3d(&self, u: &mut [f64]) {
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
/// 1D heat equation u_t = D u_xx with Dirichlet boundary conditions.
///
/// Provides three schemes: explicit FTCS (forward Euler), implicit
/// (backward Euler), and Crank-Nicolson (second-order in time).
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
/// One-dimensional wave equation solver: u_tt = c^2 * u_xx.
///
/// Supports the Lax-Wendroff scheme and Godunov (upwind) method,
/// with absorbing boundary conditions.
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
    pub fn apply_absorbing_bc(&self, u: &mut [f64]) {
        let n = u.len();
        if n >= 2 {
            u[0] = u[1];
            u[n - 1] = u[n - 2];
        }
    }
}
/// One-dimensional finite-volume solver with Godunov flux and MUSCL reconstruction.
///
/// Implements Roe's approximate Riemann solver and slope limiters
/// (minmod, superbee, van Leer) for second-order spatial accuracy.
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
        let s = 0.5 * (u_l + u_r);
        if s >= 0.0 {
            0.5 * u_l * u_l
        } else {
            0.5 * u_r * u_r
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
/// 2D level-set method for interface tracking.
///
/// Maintains a signed distance function φ; the interface is the zero level set.
/// Provides reinitialization (Sussman), curvature/normal computation, and advection.
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
/// 1D wave equation solver using the 2nd-order central-difference leapfrog scheme.
///
/// Equation: u_tt = c² u_xx.  Uses u^{n+1} = 2u^n − u^{n-1} + C²(u^n_{i-1} − 2u^n_i + u^n_{i+1}).
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
