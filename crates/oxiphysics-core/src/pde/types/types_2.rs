//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::functions::*;
use std::f64::consts::PI;

use super::types_impl::BoundaryCondition;

/// Alternating Direction Implicit (ADI) method for 2D parabolic PDEs.
///
/// Peaceman-Rachford ADI: each full time step is split into an x-implicit half-step
/// and a y-implicit half-step, giving second-order accuracy and unconditional stability.
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
            u_half[j * nx..j * nx + nx].copy_from_slice(&row);
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
/// Solver for the viscous Burgers equation: u_t + u*u_x = ν*u_xx.
///
/// Uses upwind differencing with entropy fix for the convective term
/// and central differences for the diffusive term.
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
/// Chebyshev pseudospectral method for PDEs on \[-1, 1\].
///
/// Provides Gauss-Lobatto nodes, the spectral differentiation matrix,
/// and Chebyshev expansion/interpolation.
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
        for (k, ak) in a.iter_mut().enumerate() {
            let c_k = if k == 0 || k == n - 1 { 2.0 } else { 1.0 };
            let mut sum = 0.0;
            for (j, &fj) in f.iter().enumerate() {
                let c_j = if j == 0 || j == n - 1 { 2.0 } else { 1.0 };
                sum += fj / c_j * (PI * k as f64 * j as f64 / (n - 1) as f64).cos();
            }
            *ak = 2.0 * sum / (c_k * (n - 1) as f64);
        }
        a
    }
}
/// Gray-Scott reaction-diffusion system on a periodic 2D grid.
///
/// Models the two-species reaction:
/// - u_t = Du ∇²u − u v² + f (1 − u)
/// - v_t = Dv ∇²v + u v² − (f + k) v
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
/// 1D finite-difference operator set (gradient, Laplacian).
///
/// Provides central-difference stencils at interior points and
/// one-sided differences at boundaries.
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
/// Inviscid (or low-viscosity) Burgers equation with shock formation.
///
/// u_t + u u_x = ν u_xx.  Uses a conservative Godunov flux for the
/// inviscid part and explicit central differences for viscosity.
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
/// 1D FitzHugh-Nagumo excitable-medium model.
///
/// v_t = D v_xx + v − v³/3 − w + I_ext
/// w_t = ε (v + a − b w)
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
/// Method of lines: semi-discretizes a PDE in space, yielding an ODE system.
///
/// The spatial operator is provided as a closure; the resulting ODE is
/// integrated using RK4 or DOPRI5 (adaptive).
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
/// One-dimensional finite-difference PDE solver.
///
/// Supports explicit (FTCS), implicit (backward Euler), and Crank-Nicolson
/// schemes for the diffusion equation, and upwind/central schemes for
/// the advection equation.
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
    fn apply_bc_1d(&self, u_new: &mut [f64], _u_old: &[f64]) {
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
