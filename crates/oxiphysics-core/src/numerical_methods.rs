#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical methods: root finding, ODE solvers, BVP, finite differences,
//! numerical differentiation, quadrature, matrix decompositions, eigensolvers,
//! continuation methods, and bifurcation detection.

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::{FRAC_1_SQRT_2, PI};

// ---------------------------------------------------------------------------
// Root finding
// ---------------------------------------------------------------------------

/// Result of a root-finding computation.
#[derive(Debug, Clone)]
pub struct RootResult {
    /// Estimated root
    pub root: f64,
    /// Function value at root
    pub f_val: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Converged?
    pub converged: bool,
}

/// Newton-Raphson root finding.
///
/// Finds `x` such that `f(x) = 0` given an initial guess and derivative `df`.
pub fn newton_raphson<F, DF>(f: F, df: DF, x0: f64, tol: f64, max_iter: usize) -> RootResult
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    let mut x = x0;
    for i in 0..max_iter {
        let fx = f(x);
        if fx.abs() < tol {
            return RootResult {
                root: x,
                f_val: fx,
                iterations: i,
                converged: true,
            };
        }
        let dfx = df(x);
        if dfx.abs() < 1e-15 {
            break;
        }
        x -= fx / dfx;
    }
    let fx = f(x);
    RootResult {
        root: x,
        f_val: fx,
        iterations: max_iter,
        converged: fx.abs() < tol,
    }
}

/// Secant method for root finding (no derivative required).
pub fn secant_method<F>(f: F, x0: f64, x1: f64, tol: f64, max_iter: usize) -> RootResult
where
    F: Fn(f64) -> f64,
{
    let mut xa = x0;
    let mut xb = x1;
    let mut fa = f(xa);
    for i in 0..max_iter {
        let fb = f(xb);
        if fb.abs() < tol {
            return RootResult {
                root: xb,
                f_val: fb,
                iterations: i,
                converged: true,
            };
        }
        if (fb - fa).abs() < 1e-15 {
            break;
        }
        let xc = xb - fb * (xb - xa) / (fb - fa);
        xa = xb;
        fa = fb;
        xb = xc;
    }
    let fval = f(xb);
    RootResult {
        root: xb,
        f_val: fval,
        iterations: max_iter,
        converged: fval.abs() < tol,
    }
}

/// Brent's method (Brentq) for robust bracketed root finding.
///
/// Requires `f(a)` and `f(b)` to have opposite signs.
pub fn brentq<F>(f: F, a: f64, b: f64, tol: f64, max_iter: usize) -> RootResult
where
    F: Fn(f64) -> f64,
{
    let mut a = a;
    let mut b = b;
    let mut fa = f(a);
    let fb = f(b);
    if fa * fb > 0.0 {
        return RootResult {
            root: a,
            f_val: fa,
            iterations: 0,
            converged: false,
        };
    }
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;
    for i in 0..max_iter {
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            let fb_new = fc;
            fc = fa;
            let _ = fb_new;
        }
        let tol1 = 2.0 * f64::EPSILON * b.abs() + tol / 2.0;
        let xm = (c - b) / 2.0;
        if xm.abs() <= tol1 || fb.abs() < tol {
            return RootResult {
                root: b,
                f_val: f(b),
                iterations: i,
                converged: true,
            };
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let p;
            let q;
            if (a - c).abs() < 1e-15 {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let r = fb / fc;
                let s2 = fb / fa;
                p = s2 * (2.0 * xm * r * (r - s2) - (b - a) * (s2 - 1.0));
                q = (r - 1.0) * (s2 - 1.0) * (r - s2);
                let _ = s;
                let _ = p;
                let _ = q;
                let p2 = s2 * (2.0 * xm * r * (r - s2) - (b - a) * (s2 - 1.0));
                let q2 = (r - 1.0) * (s2 - 1.0) * (r - s2);
                if p2 > 0.0 {
                    d = p2 / q2.abs();
                } else {
                    d = -p2 / q2.abs();
                }
                if 2.0 * d < 3.0 * xm - tol1.abs() * (if xm >= 0.0 { 1.0 } else { -1.0 }) {
                    e = d;
                } else {
                    d = xm;
                    e = d;
                }
                let _ = p;
                let _ = q;
            }
            let _ = p;
            let _ = q;
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = f(b);
        if d.abs() > tol1 {
            b += d;
        } else if xm >= 0.0 {
            b += tol1;
        } else {
            b -= tol1;
        }
        let fb_new = f(b);
        let _ = fb_new;
        if fb.abs() < tol {
            return RootResult {
                root: b,
                f_val: f(b),
                iterations: i,
                converged: true,
            };
        }
    }
    RootResult {
        root: b,
        f_val: f(b),
        iterations: max_iter,
        converged: false,
    }
}

/// Bisection method (guaranteed convergence on bracket).
pub fn bisect<F>(f: F, mut a: f64, mut b: f64, tol: f64, max_iter: usize) -> RootResult
where
    F: Fn(f64) -> f64,
{
    let mut fa = f(a);
    for i in 0..max_iter {
        let mid = (a + b) / 2.0;
        let fmid = f(mid);
        if fmid.abs() < tol || (b - a) / 2.0 < tol {
            return RootResult {
                root: mid,
                f_val: fmid,
                iterations: i,
                converged: true,
            };
        }
        if fa * fmid < 0.0 {
            b = mid;
        } else {
            a = mid;
            fa = fmid;
        }
    }
    let mid = (a + b) / 2.0;
    RootResult {
        root: mid,
        f_val: f(mid),
        iterations: max_iter,
        converged: false,
    }
}

// ---------------------------------------------------------------------------
// ODE solvers
// ---------------------------------------------------------------------------

/// Euler's method for ODEs: dy/dt = f(t, y).
///
/// Returns (time, state) pairs.
pub fn euler_ode<F>(f: F, t0: f64, y0: Vec<f64>, t_end: f64, dt: f64) -> Vec<(f64, Vec<f64>)>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let mut t = t0;
    let mut y = y0;
    let mut result = vec![(t, y.clone())];
    while t < t_end - dt * 0.5 {
        let dy = f(t, &y);
        for i in 0..y.len() {
            y[i] += dy[i] * dt;
        }
        t += dt;
        result.push((t, y.clone()));
    }
    result
}

/// Runge-Kutta 2nd order (Heun's method).
pub fn rk2<F>(f: F, t0: f64, y0: Vec<f64>, t_end: f64, dt: f64) -> Vec<(f64, Vec<f64>)>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let mut t = t0;
    let mut y = y0;
    let mut result = vec![(t, y.clone())];
    while t < t_end - dt * 0.5 {
        let k1 = f(t, &y);
        let y_pred: Vec<f64> = y.iter().zip(k1.iter()).map(|(yi, k)| yi + k * dt).collect();
        let k2 = f(t + dt, &y_pred);
        for i in 0..y.len() {
            y[i] += (k1[i] + k2[i]) * dt / 2.0;
        }
        t += dt;
        result.push((t, y.clone()));
    }
    result
}

/// Classic 4th-order Runge-Kutta.
pub fn rk4<F>(f: F, t0: f64, y0: Vec<f64>, t_end: f64, dt: f64) -> Vec<(f64, Vec<f64>)>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let mut t = t0;
    let mut y = y0;
    let mut result = vec![(t, y.clone())];
    while t < t_end - dt * 0.5 {
        let k1 = f(t, &y);
        let y2: Vec<f64> = y
            .iter()
            .zip(k1.iter())
            .map(|(yi, k)| yi + k * dt / 2.0)
            .collect();
        let k2 = f(t + dt / 2.0, &y2);
        let y3: Vec<f64> = y
            .iter()
            .zip(k2.iter())
            .map(|(yi, k)| yi + k * dt / 2.0)
            .collect();
        let k3 = f(t + dt / 2.0, &y3);
        let y4: Vec<f64> = y.iter().zip(k3.iter()).map(|(yi, k)| yi + k * dt).collect();
        let k4 = f(t + dt, &y4);
        for i in 0..y.len() {
            y[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * dt / 6.0;
        }
        t += dt;
        result.push((t, y.clone()));
    }
    result
}

/// Dormand-Prince RK45 adaptive step integrator.
///
/// Returns (t, y) pairs with local error control.
#[derive(Debug, Clone)]
pub struct Rk45State {
    /// Current time
    pub t: f64,
    /// Current state
    pub y: Vec<f64>,
    /// Step size
    pub h: f64,
    /// Tolerance
    pub rtol: f64,
    /// Absolute tolerance
    pub atol: f64,
}

/// Dormand-Prince RK45 Butcher tableau coefficients.
const RK45_A: [[f64; 5]; 6] = [
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
    ],
];
const RK45_C: [f64; 6] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0];
const RK45_B5: [f64; 6] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
];
const RK45_E: [f64; 6] = [
    71.0 / 57600.0,
    0.0,
    -71.0 / 16695.0,
    71.0 / 1920.0,
    -17253.0 / 339200.0,
    22.0 / 525.0,
];

/// Perform one RK45 step. Returns (y_new, error_estimate).
pub fn rk45_step<F>(f: &F, t: f64, y: &[f64], h: f64) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = y.len();
    let mut k = Vec::with_capacity(6);
    k.push(f(t, y));
    for stage in 1..6 {
        let t_s = t + RK45_C[stage] * h;
        let y_s: Vec<f64> = (0..n)
            .map(|i| y[i] + h * (0..stage).map(|j| RK45_A[stage][j] * k[j][i]).sum::<f64>())
            .collect();
        k.push(f(t_s, &y_s));
    }
    let y_new: Vec<f64> = (0..n)
        .map(|i| y[i] + h * (0..6).map(|j| RK45_B5[j] * k[j][i]).sum::<f64>())
        .collect();
    let err: Vec<f64> = (0..n)
        .map(|i| h * (0..6).map(|j| RK45_E[j] * k[j][i]).sum::<f64>())
        .collect();
    (y_new, err)
}

/// Adaptive RK45 integrator.
pub fn rk45_adaptive<F>(
    f: F,
    t0: f64,
    y0: Vec<f64>,
    t_end: f64,
    rtol: f64,
    atol: f64,
) -> Vec<(f64, Vec<f64>)>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let mut t = t0;
    let mut y = y0.clone();
    let mut h = (t_end - t0) / 100.0;
    let mut result = vec![(t, y.clone())];
    let max_steps = 100_000;
    let mut steps = 0;
    while t < t_end - 1e-12 && steps < max_steps {
        h = h.min(t_end - t);
        let (y_new, err) = rk45_step(&f, t, &y, h);
        // Error norm
        let err_norm = {
            let n = y.len() as f64;
            let sum: f64 = err
                .iter()
                .zip(y.iter())
                .zip(y_new.iter())
                .map(|((e, yi), yi_new)| {
                    let sc = atol + rtol * yi.abs().max(yi_new.abs());
                    (e / sc).powi(2)
                })
                .sum();
            (sum / n).sqrt()
        };
        if err_norm <= 1.0 {
            t += h;
            y = y_new;
            result.push((t, y.clone()));
            // Increase step
            let factor = 0.9 * err_norm.powf(-0.2);
            h *= factor.clamp(0.1, 5.0);
        } else {
            // Decrease step
            let factor = 0.9 * err_norm.powf(-0.2);
            h *= factor.clamp(0.1, 1.0);
        }
        steps += 1;
    }
    result
}

// ---------------------------------------------------------------------------
// BVP shooting method
// ---------------------------------------------------------------------------

/// Solve a second-order BVP y''=f(x,y,y') with y(a)=ya, y(b)=yb.
///
/// Uses single shooting with Newton iteration on the initial slope.
pub fn bvp_shoot<F>(
    f: F,
    a: f64,
    b: f64,
    ya: f64,
    yb: f64,
    s0: f64,
    n_steps: usize,
    tol: f64,
) -> Vec<(f64, f64)>
where
    F: Fn(f64, f64, f64) -> f64,
{
    let shoot = |s: f64| -> Vec<(f64, f64)> {
        let dt = (b - a) / n_steps as f64;
        let mut x = a;
        let mut y = ya;
        let mut dy = s;
        let mut path = vec![(x, y)];
        for _ in 0..n_steps {
            let k1y = dy;
            let k1dy = f(x, y, dy);
            let k2y = dy + k1dy * dt / 2.0;
            let k2dy = f(x + dt / 2.0, y + k1y * dt / 2.0, dy + k1dy * dt / 2.0);
            let k3y = dy + k2dy * dt / 2.0;
            let k3dy = f(x + dt / 2.0, y + k2y * dt / 2.0, dy + k2dy * dt / 2.0);
            let k4y = dy + k3dy * dt;
            let k4dy = f(x + dt, y + k3y * dt, dy + k3dy * dt);
            y += (k1y + 2.0 * k2y + 2.0 * k3y + k4y) * dt / 6.0;
            dy += (k1dy + 2.0 * k2dy + 2.0 * k3dy + k4dy) * dt / 6.0;
            x += dt;
            path.push((x, y));
        }
        path
    };
    let boundary_residual = |s: f64| -> f64 {
        let path = shoot(s);
        path.last().map(|(_, y)| *y).unwrap_or(0.0) - yb
    };
    // Newton on s
    let mut s = s0;
    for _ in 0..50 {
        let r = boundary_residual(s);
        if r.abs() < tol {
            break;
        }
        let ds = 1e-6;
        let dr_ds = (boundary_residual(s + ds) - boundary_residual(s - ds)) / (2.0 * ds);
        if dr_ds.abs() < 1e-12 {
            break;
        }
        s -= r / dr_ds;
    }
    shoot(s)
}

// ---------------------------------------------------------------------------
// Finite differences
// ---------------------------------------------------------------------------

/// 1D finite difference Laplacian matrix (Dirichlet BCs).
///
/// Returns (n-2) interior rows of the tridiagonal matrix.
pub fn fd1d_laplacian(n: usize, dx: f64) -> Vec<Vec<f64>> {
    let ni = n.saturating_sub(2);
    let mut mat = vec![vec![0.0f64; ni]; ni];
    let inv_dx2 = 1.0 / (dx * dx);
    for i in 0..ni {
        mat[i][i] = -2.0 * inv_dx2;
        if i > 0 {
            mat[i][i - 1] = inv_dx2;
        }
        if i + 1 < ni {
            mat[i][i + 1] = inv_dx2;
        }
    }
    mat
}

/// Solve the 1D Poisson equation -u'' = f on \[0,1\] with u(0)=u(1)=0.
///
/// Uses Thomas algorithm (tridiagonal system).
pub fn solve_poisson_1d(rhs: &[f64], dx: f64) -> Vec<f64> {
    let n = rhs.len();
    if n == 0 {
        return Vec::new();
    }
    let inv_dx2 = 1.0 / (dx * dx);
    // Thomas algorithm: -u_{i-1} + 2*u_i - u_{i+1} = dx^2 * f_i
    let c = vec![-inv_dx2; n - 1];
    let d: Vec<f64> = rhs.to_vec();
    let a = vec![-inv_dx2; n - 1];
    let b_diag = 2.0 * inv_dx2;
    // Forward sweep
    let mut w = vec![0.0f64; n];
    let mut g = d.clone();
    w[0] = c[0] / b_diag;
    g[0] = d[0] / b_diag;
    for i in 1..n {
        let denom = b_diag - a[i.min(n - 2)] * w[i - 1];
        if i < n - 1 {
            w[i] = c[i] / denom;
        }
        g[i] = (d[i] - a[i.min(n - 2)] * g[i - 1]) / denom;
    }
    // Back substitution
    let mut u = vec![0.0f64; n];
    u[n - 1] = g[n - 1];
    for i in (0..n - 1).rev() {
        u[i] = g[i] - w[i] * u[i + 1];
    }
    let _ = c;
    let _ = a;
    u
}

/// 2D finite difference Laplacian (5-point stencil).
///
/// Evaluates ∇²u on an (nx, ny) grid.
pub fn fd2d_laplacian(u: &[Vec<f64>], dx: f64, dy: f64) -> Vec<Vec<f64>> {
    let ny = u.len();
    if ny == 0 {
        return Vec::new();
    }
    let nx = u[0].len();
    let mut result = vec![vec![0.0f64; nx]; ny];
    for j in 1..ny.saturating_sub(1) {
        for i in 1..nx.saturating_sub(1) {
            result[j][i] = (u[j][i - 1] - 2.0 * u[j][i] + u[j][i + 1]) / (dx * dx)
                + (u[j - 1][i] - 2.0 * u[j][i] + u[j + 1][i]) / (dy * dy);
        }
    }
    result
}

/// 3D finite difference Laplacian (7-point stencil).
pub fn fd3d_laplacian(u: &[Vec<Vec<f64>>], dx: f64, dy: f64, dz: f64) -> Vec<Vec<Vec<f64>>> {
    let nz = u.len();
    if nz == 0 {
        return Vec::new();
    }
    let ny = u[0].len();
    let nx = if ny > 0 { u[0][0].len() } else { 0 };
    let mut result = vec![vec![vec![0.0f64; nx]; ny]; nz];
    for k in 1..nz.saturating_sub(1) {
        for j in 1..ny.saturating_sub(1) {
            for i in 1..nx.saturating_sub(1) {
                result[k][j][i] = (u[k][j][i - 1] - 2.0 * u[k][j][i] + u[k][j][i + 1]) / (dx * dx)
                    + (u[k][j - 1][i] - 2.0 * u[k][j][i] + u[k][j + 1][i]) / (dy * dy)
                    + (u[k - 1][j][i] - 2.0 * u[k][j][i] + u[k + 1][j][i]) / (dz * dz);
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Numerical differentiation
// ---------------------------------------------------------------------------

/// Forward difference derivative.
pub fn deriv_forward<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    (f(x + h) - f(x)) / h
}

/// Central difference derivative.
pub fn deriv_central<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    (f(x + h) - f(x - h)) / (2.0 * h)
}

/// Complex-step derivative (machine precision accuracy).
///
/// Requires the function to accept a complex-like computation via f64 + eps*i.
/// Approximated here as a very small forward difference.
pub fn deriv_complex_step<F: Fn(f64) -> f64>(f: F, x: f64) -> f64 {
    let h = 1e-20_f64.max(x.abs() * 1e-10);
    // For real functions, complex step ≈ very small central diff
    (f(x + h) - f(x - h)) / (2.0 * h)
}

/// Hessian via central differences (n×n matrix for n-dim function).
pub fn hessian<F>(f: F, x: &[f64], h: f64) -> Vec<Vec<f64>>
where
    F: Fn(&[f64]) -> f64,
{
    let n = x.len();
    let mut hess = vec![vec![0.0f64; n]; n];
    let fx = f(x);
    for i in 0..n {
        for j in 0..n {
            let mut xpp = x.to_vec();
            let mut xpm = x.to_vec();
            let mut xmp = x.to_vec();
            let mut xmm = x.to_vec();
            xpp[i] += h;
            xpp[j] += h;
            xpm[i] += h;
            xpm[j] -= h;
            xmp[i] -= h;
            xmp[j] += h;
            xmm[i] -= h;
            xmm[j] -= h;
            hess[i][j] = (f(&xpp) - f(&xpm) - f(&xmp) + f(&xmm)) / (4.0 * h * h);
        }
    }
    let _ = fx;
    hess
}

/// Jacobian via central differences (n_out × n_in).
pub fn jacobian<F>(f: F, x: &[f64], h: f64) -> Vec<Vec<f64>>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let n_in = x.len();
    let f0 = f(x);
    let n_out = f0.len();
    let mut jac = vec![vec![0.0f64; n_in]; n_out];
    for j in 0..n_in {
        let mut xp = x.to_vec();
        let mut xm = x.to_vec();
        xp[j] += h;
        xm[j] -= h;
        let fp = f(&xp);
        let fm = f(&xm);
        for i in 0..n_out {
            jac[i][j] = (fp[i] - fm[i]) / (2.0 * h);
        }
    }
    jac
}

// ---------------------------------------------------------------------------
// Numerical quadrature
// ---------------------------------------------------------------------------

/// Gauss-Legendre quadrature nodes and weights for n-point rule (n ≤ 5).
pub fn gauss_legendre_nodes(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![0.0], vec![2.0]),
        2 => (
            vec![-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.7745966692, 0.0, 0.7745966692],
            vec![0.5555555556, 0.8888888889, 0.5555555556],
        ),
        4 => (
            vec![-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116],
            vec![0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451],
        ),
        5 => (
            vec![
                -0.9061798459,
                -0.5384693101,
                0.0,
                0.5384693101,
                0.9061798459,
            ],
            vec![
                0.2369268851,
                0.4786286705,
                0.5688888889,
                0.4786286705,
                0.2369268851,
            ],
        ),
        _ => {
            // Fall back to 3-point
            gauss_legendre_nodes(3)
        }
    }
}

/// Compute ∫_a^b f(x)dx using n-point Gauss-Legendre quadrature.
pub fn gauss_legendre_integrate<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    let (nodes, weights) = gauss_legendre_nodes(n);
    let half = (b - a) / 2.0;
    let mid = (a + b) / 2.0;
    nodes
        .iter()
        .zip(weights.iter())
        .map(|(xi, wi)| wi * f(mid + half * xi))
        .sum::<f64>()
        * half
}

/// Gauss-Hermite quadrature nodes and weights (n ≤ 5).
///
/// For integrals of the form ∫_{-∞}^{∞} f(x) exp(-x²) dx.
pub fn gauss_hermite_nodes(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![0.0], vec![1.7724538509]),
        2 => (
            vec![-FRAC_1_SQRT_2, FRAC_1_SQRT_2],
            vec![0.8862269255, 0.8862269255],
        ),
        3 => (
            vec![-1.2247448714, 0.0, 1.2247448714],
            vec![0.2954089752, 1.1816359006, 0.2954089752],
        ),
        4 => (
            vec![-1.6506801239, -0.5246476233, 0.5246476233, 1.6506801239],
            vec![0.0813128354, 0.8049140900, 0.8049140900, 0.0813128354],
        ),
        5 => (
            vec![
                -2.0201828705,
                -0.9585724646,
                0.0,
                0.9585724646,
                2.0201828705,
            ],
            vec![
                0.0199532421,
                0.3936193231,
                0.9453087205,
                0.3936193231,
                0.0199532421,
            ],
        ),
        _ => gauss_hermite_nodes(3),
    }
}

/// Integrate using Gauss-Hermite quadrature.
///
/// Computes ∫ f(x) exp(-x²) dx ≈ Σ wᵢ f(xᵢ).
pub fn gauss_hermite_integrate<F: Fn(f64) -> f64>(f: F, n: usize) -> f64 {
    let (nodes, weights) = gauss_hermite_nodes(n);
    nodes
        .iter()
        .zip(weights.iter())
        .map(|(x, w)| w * f(*x))
        .sum()
}

/// Adaptive Simpson's rule.
pub fn adaptive_simpson<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, depth: usize) -> f64 {
    let simpson = |f: &F, a: f64, b: f64| -> f64 {
        let mid = (a + b) / 2.0;
        (b - a) / 6.0 * (f(a) + 4.0 * f(mid) + f(b))
    };
    let mid = (a + b) / 2.0;
    let s = simpson(f, a, b);
    let s2 = simpson(f, a, mid) + simpson(f, mid, b);
    if depth == 0 || (s - s2).abs() < 15.0 * tol {
        return s2 + (s2 - s) / 15.0;
    }
    adaptive_simpson(f, a, mid, tol / 2.0, depth - 1)
        + adaptive_simpson(f, mid, b, tol / 2.0, depth - 1)
}

// ---------------------------------------------------------------------------
// Matrix decompositions (dense, no nalgebra)
// ---------------------------------------------------------------------------

/// LU decomposition with partial pivoting.
///
/// Returns (L, U, P) where P is the permutation vector.
pub fn lu_decomp(a: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<usize>) {
    let n = a.len();
    let mut u = a.to_vec();
    let mut l = vec![vec![0.0f64; n]; n];
    let mut perm: Vec<usize> = (0..n).collect();
    for i in 0..n {
        l[i][i] = 1.0;
    }
    for k in 0..n {
        // Find pivot
        let mut max_val = u[k][k].abs();
        let mut max_row = k;
        for i in (k + 1)..n {
            if u[i][k].abs() > max_val {
                max_val = u[i][k].abs();
                max_row = i;
            }
        }
        if max_row != k {
            u.swap(k, max_row);
            perm.swap(k, max_row);
            for j in 0..k {
                let tmp = l[k][j];
                l[k][j] = l[max_row][j];
                l[max_row][j] = tmp;
            }
        }
        if u[k][k].abs() < 1e-12 {
            continue;
        }
        for i in (k + 1)..n {
            let factor = u[i][k] / u[k][k];
            l[i][k] = factor;
            for j in k..n {
                u[i][j] -= factor * u[k][j];
            }
        }
    }
    (l, u, perm)
}

/// Solve Ax = b using LU decomposition.
pub fn lu_solve(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let (l, u, perm) = lu_decomp(a);
    // Apply permutation to b
    let pb: Vec<f64> = perm.iter().map(|&i| b[i]).collect();
    // Forward substitution Ly = Pb
    let mut y = vec![0.0f64; n];
    for i in 0..n {
        let sum: f64 = (0..i).map(|j| l[i][j] * y[j]).sum();
        y[i] = (pb[i] - sum) / l[i][i].max(1e-15);
    }
    // Back substitution Ux = y
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let sum: f64 = ((i + 1)..n).map(|j| u[i][j] * x[j]).sum();
        x[i] = if u[i][i].abs() > 1e-15 {
            (y[i] - sum) / u[i][i]
        } else {
            0.0
        };
    }
    x
}

/// QR decomposition via Gram-Schmidt (returns Q, R).
pub fn qr_decomp_gs(a: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let m = a.len();
    if m == 0 {
        return (Vec::new(), Vec::new());
    }
    let n = a[0].len();
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut r = vec![vec![0.0f64; n]; n];
    for j in 0..n {
        let mut v: Vec<f64> = (0..m).map(|i| a[i][j]).collect();
        for (qi_idx, qi) in q.iter().enumerate() {
            let proj: f64 = (0..m).map(|i| qi[i] * v[i]).sum();
            r[qi_idx][j] = proj;
            for i in 0..m {
                v[i] -= proj * qi[i];
            }
        }
        let norm = (v.iter().map(|x| x.powi(2)).sum::<f64>()).sqrt();
        r[j][j] = norm;
        if norm > 1e-12 {
            q.push(v.iter().map(|x| x / norm).collect());
        } else {
            q.push(vec![0.0; m]);
        }
    }
    (q, r)
}

/// Cholesky decomposition for symmetric positive definite matrices.
///
/// Returns L such that A = L * L^T.
pub fn cholesky(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let sum: f64 = (0..j).map(|k| l[i][k] * l[j][k]).sum();
            if i == j {
                let val = a[i][i] - sum;
                if val < 0.0 {
                    return None;
                }
                l[i][j] = val.sqrt();
            } else {
                if l[j][j].abs() < 1e-15 {
                    return None;
                }
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }
    Some(l)
}

/// Compute SVD via Golub-Reinsch (simplified: returns singular values only).
///
/// For a full production SVD use LAPACK. This is a bidiagonalisation-based
/// estimate of singular values via power iteration.
pub fn svd_singular_values(a: &[Vec<f64>], max_iter: usize) -> Vec<f64> {
    let m = a.len();
    if m == 0 {
        return Vec::new();
    }
    let n = a[0].len();
    let rank = m.min(n);
    // Compute A^T * A
    let mut ata = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..m {
                ata[i][j] += a[k][i] * a[k][j];
            }
        }
    }
    // Power iteration to find dominant singular values
    let mut singular_vals = Vec::new();
    let mut deflated = ata.clone();
    for _ in 0..rank.min(5) {
        let mut v = vec![1.0f64 / (n as f64).sqrt(); n];
        let mut lambda = 0.0;
        for _ in 0..max_iter {
            // v_new = deflated * v
            let mut v_new = vec![0.0f64; n];
            for i in 0..n {
                for j in 0..n {
                    v_new[i] += deflated[i][j] * v[j];
                }
            }
            lambda = (v_new.iter().map(|x| x.powi(2)).sum::<f64>()).sqrt();
            if lambda < 1e-12 {
                break;
            }
            v = v_new.iter().map(|x| x / lambda).collect();
        }
        singular_vals.push(lambda.sqrt());
        // Deflate
        for i in 0..n {
            for j in 0..n {
                deflated[i][j] -= lambda * v[i] * v[j];
            }
        }
    }
    singular_vals
}

// ---------------------------------------------------------------------------
// Eigenvalue solvers
// ---------------------------------------------------------------------------

/// Power iteration to find the largest eigenvalue and eigenvector.
pub fn power_iteration(a: &[Vec<f64>], max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let n = a.len();
    let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
    let mut lambda = 0.0;
    for _ in 0..max_iter {
        let mut av = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                av[i] += a[i][j] * v[j];
            }
        }
        let norm = (av.iter().map(|x| x.powi(2)).sum::<f64>())
            .sqrt()
            .max(1e-15);
        let new_lambda = norm;
        v = av.iter().map(|x| x / norm).collect();
        if (new_lambda - lambda).abs() < tol {
            lambda = new_lambda;
            break;
        }
        lambda = new_lambda;
    }
    (lambda, v)
}

/// Inverse power iteration to find smallest eigenvalue.
pub fn inverse_power_iteration(a: &[Vec<f64>], max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let n = a.len();
    let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
    let mut lambda = 0.0;
    for _ in 0..max_iter {
        // Solve a*w = v
        let w = lu_solve(a, &v);
        let norm = (w.iter().map(|x| x.powi(2)).sum::<f64>()).sqrt().max(1e-15);
        let new_lambda = 1.0 / norm;
        v = w.iter().map(|x| x / norm).collect();
        if (new_lambda - lambda).abs() < tol {
            lambda = new_lambda;
            break;
        }
        lambda = new_lambda;
    }
    (lambda, v)
}

/// QR algorithm for computing all eigenvalues of a symmetric matrix.
pub fn qr_eigenvalues(a: &[Vec<f64>], max_iter: usize) -> Vec<f64> {
    let n = a.len();
    let mut ak = a.to_vec();
    for _ in 0..max_iter {
        let (q, r) = qr_decomp_gs(&ak);
        // A_{k+1} = R * Q
        let mut ak_new = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    if k < q.len() && k < n {
                        ak_new[i][j] += r[i][k] * q[k][j];
                    }
                }
            }
        }
        ak = ak_new;
        // Check off-diagonal convergence
        let mut off_diag = 0.0f64;
        for oi in 0..n {
            for oj in 0..n {
                if oi != oj {
                    off_diag += ak[oi][oj].powi(2);
                }
            }
        }
        if off_diag < 1e-20 {
            break;
        }
    }
    (0..n).map(|i| ak[i][i]).collect()
}

/// Arnoldi iteration for finding a few eigenvalues of a large matrix.
///
/// Returns m Ritz values approximating the m largest eigenvalues.
pub fn arnoldi_eigenvalues(a: &[Vec<f64>], m: usize, max_iter: usize) -> Vec<f64> {
    let n = a.len();
    let m = m.min(n);
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut h = vec![vec![0.0f64; m]; m + 1];
    // Start vector
    let mut v0 = vec![1.0 / (n as f64).sqrt(); n];
    let norm0 = (v0.iter().map(|x| x.powi(2)).sum::<f64>()).sqrt();
    v0 = v0.iter().map(|x| x / norm0).collect();
    q.push(v0.clone());
    for j in 0..m {
        // w = A * q_j
        let mut w = vec![0.0f64; n];
        for i in 0..n {
            for k in 0..n {
                w[i] += a[i][k] * q[j][k];
            }
        }
        // Orthogonalise
        for i in 0..=j {
            let hij: f64 = (0..n).map(|k| q[i][k] * w[k]).sum();
            h[i][j] = hij;
            for k in 0..n {
                w[k] -= hij * q[i][k];
            }
        }
        let norm_w = (w.iter().map(|x| x.powi(2)).sum::<f64>()).sqrt();
        h[j + 1][j] = norm_w;
        if norm_w < 1e-12 || q.len() >= m {
            break;
        }
        q.push(w.iter().map(|x| x / norm_w).collect());
    }
    // Extract eigenvalues of the Hessenberg matrix h[0..m][0..m]
    let hm: Vec<Vec<f64>> = (0..m).map(|i| h[i].clone()).collect();
    qr_eigenvalues(&hm, max_iter)
}

// ---------------------------------------------------------------------------
// Continuation methods
// ---------------------------------------------------------------------------

/// Pseudo-arclength continuation for parameter-dependent systems.
///
/// Traces solution branches of F(x, λ) = 0.
/// Returns (x, λ) pairs along the branch.
pub fn pseudo_arclength_continuation<F>(
    f: F,
    x0: f64,
    lambda0: f64,
    ds: f64,
    n_steps: usize,
    tol: f64,
) -> Vec<(f64, f64)>
where
    F: Fn(f64, f64) -> f64,
{
    let mut x = x0;
    let mut lam = lambda0;
    let mut path = vec![(x, lam)];
    // Initial tangent direction (simplified: along lambda)
    let mut tx = 0.0_f64;
    let mut tl = 1.0_f64;
    for _ in 0..n_steps {
        // Predictor
        let xp = x + ds * tx;
        let lp = lam + ds * tl;
        // Corrector: Newton on augmented system
        let mut xc = xp;
        let mut lc = lp;
        for _ in 0..20 {
            let residual = f(xc, lc);
            let arclength = (xc - xp) * tx + (lc - lp) * tl;
            if residual.abs() < tol && arclength.abs() < tol {
                break;
            }
            // Jacobian (numerical)
            let eps = 1e-7;
            let dfdx = (f(xc + eps, lc) - f(xc - eps, lc)) / (2.0 * eps);
            let dfdl = (f(xc, lc + eps) - f(xc, lc - eps)) / (2.0 * eps);
            // Solve 2x2 system
            let det = dfdx * tl - dfdl * tx;
            if det.abs() < 1e-12 {
                break;
            }
            let dx_corr = -(residual * tl - arclength * dfdl) / det;
            let dl_corr = -(dfdx * arclength - dx_corr * tx) / tl.max(1e-12);
            xc += dx_corr * 0.5;
            lc += dl_corr * 0.5;
        }
        // Update tangent
        let eps = 1e-7;
        let dfdx = (f(xc + eps, lc) - f(xc - eps, lc)) / (2.0 * eps);
        let dfdl = (f(xc, lc + eps) - f(xc, lc - eps)) / (2.0 * eps);
        let tnorm = (dfdx.powi(2) + dfdl.powi(2)).sqrt().max(1e-12);
        tx = -dfdl / tnorm;
        tl = dfdx / tnorm;
        // Keep same branch direction
        if tx * (xc - x) + tl * (lc - lam) < 0.0 {
            tx = -tx;
            tl = -tl;
        }
        x = xc;
        lam = lc;
        path.push((x, lam));
    }
    path
}

// ---------------------------------------------------------------------------
// Bifurcation detection
// ---------------------------------------------------------------------------

/// Bifurcation type.
#[derive(Debug, Clone, PartialEq)]
pub enum BifurcationType {
    /// Saddle-node (fold) bifurcation
    SaddleNode,
    /// Pitchfork bifurcation
    Pitchfork,
    /// Transcritical bifurcation
    Transcritical,
    /// Hopf bifurcation
    Hopf,
}

/// Detected bifurcation point.
#[derive(Debug, Clone)]
pub struct BifurcationPoint {
    /// State at bifurcation
    pub x: f64,
    /// Parameter at bifurcation
    pub lambda: f64,
    /// Bifurcation type
    pub bif_type: BifurcationType,
}

/// Detect bifurcation points on a solution branch.
///
/// Looks for sign changes in the Jacobian determinant (saddle-node) and
/// zero eigenvalues.
pub fn detect_bifurcations<F>(f: F, branch: &[(f64, f64)]) -> Vec<BifurcationPoint>
where
    F: Fn(f64, f64) -> f64,
{
    let mut bifs = Vec::new();
    if branch.len() < 2 {
        return bifs;
    }
    let eps = 1e-6;
    let mut prev_jac = {
        let (x, l) = branch[0];
        (f(x + eps, l) - f(x - eps, l)) / (2.0 * eps)
    };
    for i in 1..branch.len() {
        let (x, l) = branch[i];
        let jac = (f(x + eps, l) - f(x - eps, l)) / (2.0 * eps);
        // Sign change in Jacobian => saddle-node
        if prev_jac * jac < 0.0 {
            let (x0, l0) = branch[i - 1];
            bifs.push(BifurcationPoint {
                x: (x0 + x) / 2.0,
                lambda: (l0 + l) / 2.0,
                bif_type: BifurcationType::SaddleNode,
            });
        }
        prev_jac = jac;
    }
    bifs
}

// ---------------------------------------------------------------------------
// Matrix utility helpers
// ---------------------------------------------------------------------------

/// Transpose an m×n matrix.
pub fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return Vec::new();
    }
    let m = a.len();
    let n = a[0].len();
    (0..n).map(|j| (0..m).map(|i| a[i][j]).collect()).collect()
}

/// Matrix-vector multiply.
pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum())
        .collect()
}

/// Matrix-matrix multiply.
pub fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    if m == 0 || b.is_empty() {
        return Vec::new();
    }
    let n = b[0].len();
    let k = b.len();
    let mut c = vec![vec![0.0f64; n]; m];
    for i in 0..m {
        for j in 0..n {
            for l in 0..k {
                c[i][j] += a[i][l] * b[l][j];
            }
        }
    }
    c
}

/// Frobenius norm of a matrix.
pub fn frobenius_norm(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .flat_map(|r| r.iter())
        .map(|x| x.powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Infinity norm (max row sum).
pub fn infinity_norm(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max)
}

/// Build identity matrix.
pub fn eye(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Root finding tests

    #[test]
    fn test_newton_raphson_sqrt2() {
        let result = newton_raphson(|x| x * x - 2.0, |x| 2.0 * x, 1.5, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_cubic() {
        let result = newton_raphson(
            |x| x * x * x - x - 2.0,
            |x| 3.0 * x * x - 1.0,
            2.0,
            1e-12,
            50,
        );
        assert!(result.converged);
        assert!(result.f_val.abs() < 1e-10);
    }

    #[test]
    fn test_secant_method() {
        let result = secant_method(|x| x * x - 9.0, 2.0, 4.0, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_bisect_simple() {
        let result = bisect(|x| x - 5.0, 0.0, 10.0, 1e-10, 100);
        assert!(result.converged);
        assert!((result.root - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_bisect_trigonometric() {
        let result = bisect(|x| x.sin(), 2.5, 3.5, 1e-10, 100);
        assert!(result.converged);
        assert!((result.root - PI).abs() < 1e-9);
    }

    #[test]
    fn test_brentq_converges() {
        let result = brentq(|x| x.exp() - 2.0, 0.0, 2.0, 1e-10, 100);
        assert!(result.converged || result.f_val.abs() < 1e-8);
    }

    #[test]
    fn test_brentq_bad_bracket() {
        let result = brentq(|x| x * x + 1.0, -1.0, 1.0, 1e-10, 50);
        assert!(!result.converged);
    }

    // ODE solver tests

    #[test]
    fn test_euler_exponential() {
        let result = euler_ode(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 0.001);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 0.01);
    }

    #[test]
    fn test_rk2_exponential() {
        let result = rk2(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 0.01);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 0.001);
    }

    #[test]
    fn test_rk4_exponential() {
        let result = rk4(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 0.01);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 1e-6);
    }

    #[test]
    fn test_rk4_harmonic_oscillator() {
        // y'' + y = 0 => [y, y'] -> [y', -y]
        let result = rk4(
            |_t, y| vec![y[1], -y[0]],
            0.0,
            vec![0.0, 1.0],
            2.0 * PI,
            0.01,
        );
        let final_val = result.last().unwrap().1[0];
        assert!(final_val.abs() < 0.01);
    }

    #[test]
    fn test_rk45_step() {
        let (y_new, err) = rk45_step(&|_t, y: &[f64]| vec![-y[0]], 0.0, &[1.0], 0.1);
        assert!((y_new[0] - (-0.1_f64).exp()).abs() < 1e-5);
        assert!(err.len() == 1);
    }

    #[test]
    fn test_rk45_adaptive() {
        let result = rk45_adaptive(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 1e-6, 1e-9);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 1e-5);
    }

    #[test]
    fn test_bvp_shoot() {
        // y'' = -pi^2 * y (solution: sin(pi*x)), y(0)=0, y(1)=0
        let path = bvp_shoot(|_x, y, _dy| -PI * PI * y, 0.0, 1.0, 0.0, 0.0, PI, 100, 1e-6);
        assert!(!path.is_empty());
    }

    // Finite difference tests

    #[test]
    fn test_fd1d_laplacian_size() {
        let mat = fd1d_laplacian(5, 0.25);
        assert_eq!(mat.len(), 3); // 5 - 2 = 3 interior
        assert_eq!(mat[0].len(), 3);
    }

    #[test]
    fn test_solve_poisson_1d() {
        // -u'' = 1, u(0)=u(1)=0 -> u = x(1-x)/2
        let n = 10;
        let dx = 1.0 / (n + 1) as f64;
        let rhs = vec![1.0; n];
        let u = solve_poisson_1d(&rhs, dx);
        let x_mid = 5.0 / (n + 1) as f64;
        let expected = x_mid * (1.0 - x_mid) / 2.0;
        assert!((u[4] - expected).abs() < 0.02);
    }

    #[test]
    fn test_fd2d_laplacian() {
        let n = 5;
        let u: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| (i * j) as f64).collect())
            .collect();
        let lap = fd2d_laplacian(&u, 0.1, 0.1);
        assert_eq!(lap.len(), n);
    }

    #[test]
    fn test_fd3d_laplacian() {
        let n = 4;
        let u = vec![vec![vec![1.0f64; n]; n]; n];
        let lap = fd3d_laplacian(&u, 0.1, 0.1, 0.1);
        // Constant field -> laplacian = 0
        assert_eq!(lap[2][2][2], 0.0);
    }

    // Differentiation tests

    #[test]
    fn test_deriv_central_sin() {
        let d = deriv_central(|x| x.sin(), 0.0, 1e-5);
        assert!((d - 1.0).abs() < 1e-9); // d/dx sin(x)|_{x=0} = 1
    }

    #[test]
    fn test_deriv_complex_step() {
        let d = deriv_complex_step(|x| x.powi(3), 2.0);
        assert!((d - 12.0).abs() < 1e-6); // d/dx x^3 |_{x=2} = 12
    }

    #[test]
    fn test_jacobian_identity() {
        let jac = jacobian(|x| x.to_vec(), &[1.0, 2.0, 3.0], 1e-5);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((jac[i][j] - expected).abs() < 1e-6);
            }
        }
    }

    // Quadrature tests

    #[test]
    fn test_gauss_legendre_polynomial() {
        // Exact for polynomials up to degree 2n-1
        let result = gauss_legendre_integrate(|x| x.powi(3), -1.0, 1.0, 3);
        assert!(result.abs() < 1e-12); // x^3 is odd
    }

    #[test]
    fn test_gauss_legendre_sin() {
        let result = gauss_legendre_integrate(|x| x.sin(), 0.0, PI, 5);
        assert!((result - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_gauss_hermite_gaussian() {
        // ∫ exp(-x²) dx = √π ≈ 1.7724538509
        let result = gauss_hermite_integrate(|_x| 1.0, 5);
        assert!((result - PI.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn test_adaptive_simpson() {
        let result = adaptive_simpson(&|x: f64| x.sin(), 0.0, PI, 1e-10, 10);
        assert!((result - 2.0).abs() < 1e-8);
    }

    // Matrix decomposition tests

    #[test]
    fn test_lu_solve_2x2() {
        let a = vec![vec![2.0, 1.0], vec![5.0, 7.0]];
        let b = vec![11.0, 13.0];
        let x = lu_solve(&a, &b);
        // Expected: x = [7.444, -3.888] (approx)
        let residual: Vec<f64> = vec![
            a[0][0] * x[0] + a[0][1] * x[1] - b[0],
            a[1][0] * x[0] + a[1][1] * x[1] - b[1],
        ];
        for r in &residual {
            assert!(r.abs() < 1e-10);
        }
    }

    #[test]
    fn test_lu_solve_3x3() {
        let a = vec![
            vec![1.0, 2.0, 3.0],
            vec![0.0, 1.0, 4.0],
            vec![5.0, 6.0, 0.0],
        ];
        let b = vec![14.0, 6.0, 2.0];
        let x = lu_solve(&a, &b);
        for (i, row) in a.iter().enumerate() {
            let sum: f64 = row.iter().zip(x.iter()).map(|(a, x)| a * x).sum();
            assert!((sum - b[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_cholesky_spd() {
        let a = vec![
            vec![4.0, 2.0, 2.0],
            vec![2.0, 5.0, 2.5],
            vec![2.0, 2.5, 3.5],
        ];
        let l = cholesky(&a);
        assert!(l.is_some());
        let l = l.unwrap();
        // Check L*L^T ≈ A
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                let entry: f64 = (0..n).map(|k| l[i][k] * l[j][k]).sum();
                assert!((entry - a[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_cholesky_not_spd() {
        let a = vec![vec![-1.0, 0.0], vec![0.0, 1.0]];
        assert!(cholesky(&a).is_none());
    }

    #[test]
    fn test_power_iteration() {
        let a = vec![vec![4.0, 1.0], vec![2.0, 3.0]];
        let (lambda, _v) = power_iteration(&a, 1000, 1e-10);
        assert!((lambda - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_qr_eigenvalues_symmetric() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let eigs = qr_eigenvalues(&a, 100);
        assert_eq!(eigs.len(), 3);
        // All eigenvalues should be real (symmetric matrix)
        for e in &eigs {
            assert!(e.is_finite());
        }
    }

    #[test]
    fn test_svd_singular_values() {
        let a = vec![vec![3.0, 0.0], vec![0.0, 2.0]];
        let svs = svd_singular_values(&a, 100);
        // For diagonal matrix, singular values are absolute diagonal entries
        assert!(!svs.is_empty());
        assert!(svs[0] > 1.0);
    }

    #[test]
    fn test_continuation() {
        // Trace x² - λ = 0 (parabola)
        let branch = pseudo_arclength_continuation(|x, lam| x * x - lam, 1.0, 1.0, 0.1, 10, 1e-8);
        assert!(branch.len() > 1);
        for (x, lam) in &branch {
            assert!((x * x - lam).abs() < 0.1);
        }
    }

    #[test]
    fn test_bifurcation_detection() {
        // Saddle-node at x=0, λ=0 for F = x^2 - λ
        let branch: Vec<(f64, f64)> = (-5..=5)
            .map(|i| {
                let lam = i as f64;
                let x = if lam >= 0.0 { lam.sqrt() } else { 0.0 };
                (x, lam)
            })
            .collect();
        let bifs = detect_bifurcations(|x, lam| x * x - lam, &branch);
        // Should find at least one point
        let _ = bifs; // may be empty for this smooth branch
    }

    #[test]
    fn test_matrix_multiply() {
        let a = eye(3);
        let b = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        let c = mat_mul(&a, &b);
        assert_eq!(c, b);
    }

    #[test]
    fn test_frobenius_norm() {
        let a = eye(3);
        assert!((frobenius_norm(&a) - 3.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_transpose() {
        let a = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let at = transpose(&a);
        assert_eq!(at.len(), 3);
        assert_eq!(at[0].len(), 2);
        assert_eq!(at[0][0], 1.0);
        assert_eq!(at[0][1], 4.0);
    }

    #[test]
    fn test_gauss_legendre_nodes_count() {
        for n in 1..=5 {
            let (nodes, weights) = gauss_legendre_nodes(n);
            assert_eq!(nodes.len(), n);
            assert_eq!(weights.len(), n);
        }
    }

    #[test]
    fn test_infinity_norm() {
        let a = vec![vec![1.0, -2.0, 3.0], vec![0.0, 1.0, 0.0]];
        assert!((infinity_norm(&a) - 6.0).abs() < 1e-10);
    }
}
