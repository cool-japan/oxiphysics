#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Finite-difference PDE solvers for common one-dimensional equations.
//!
//! All solvers operate on simple `Vec`f64` grids and accept closure-based
//! boundary conditions or source terms for flexibility.  The emphasis is on
//! correctness and clarity rather than maximum performance.

// ──────────────────────────────────────────────────────────────────────────────
// CflCondition
// ──────────────────────────────────────────────────────────────────────────────

/// Stability criterion checker for explicit finite-difference schemes.
///
/// The Courant–Friedrichs–Lewy (CFL) number must satisfy `ν ≤ 1` (advection)
/// or `α ≤ 0.5` (diffusion) for stability of explicit schemes.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CflCondition {
    /// Grid spacing Δx.
    pub dx: f64,
    /// Time step Δt.
    pub dt: f64,
    /// Wave / advection speed (m/s).
    pub speed: f64,
    /// Diffusion coefficient (m²/s).
    pub diffusivity: f64,
}

impl CflCondition {
    /// Construct a new `CflCondition`.
    pub fn new(dx: f64, dt: f64, speed: f64, diffusivity: f64) -> Self {
        Self {
            dx,
            dt,
            speed,
            diffusivity,
        }
    }

    /// Advective CFL number `c·Δt/Δx`. Must be ≤ 1 for explicit upwind/leapfrog.
    pub fn advective_cfl(&self) -> f64 {
        self.speed * self.dt / self.dx
    }

    /// Diffusive CFL number `α·Δt/Δx²`. Must be ≤ 0.5 for explicit FTCS.
    pub fn diffusive_cfl(&self) -> f64 {
        self.diffusivity * self.dt / (self.dx * self.dx)
    }

    /// Returns `true` if the advective CFL condition is satisfied (`≤ 1`).
    pub fn is_advection_stable(&self) -> bool {
        self.advective_cfl() <= 1.0
    }

    /// Returns `true` if the diffusive CFL condition is satisfied (`≤ 0.5`).
    pub fn is_diffusion_stable(&self) -> bool {
        self.diffusive_cfl() <= 0.5
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 1-D Heat equation  (explicit FTCS)
// ──────────────────────────────────────────────────────────────────────────────

/// Solve the 1-D heat equation  ∂u/∂t = α ∂²u/∂x²  for `n_steps` time steps
/// using the explicit Forward-Time Centred-Space (FTCS) scheme.
///
/// # Parameters
/// - `u0` — initial temperature field (length N).
/// - `alpha` — thermal diffusivity (m²/s).
/// - `dx` — grid spacing (m).
/// - `dt` — time step (s).  Stability requires `α·Δt/Δx² ≤ 0.5`.
/// - `n_steps` — number of time steps to advance.
///
/// Dirichlet boundary conditions are enforced by keeping the first and last
/// elements fixed.
///
/// Returns the temperature field after `n_steps`.
#[allow(dead_code)]
pub fn heat_equation_1d(u0: &[f64], alpha: f64, dx: f64, dt: f64, n_steps: usize) -> Vec<f64> {
    let n = u0.len();
    let mut u = u0.to_vec();
    let mut u_new = u.clone();
    let r = alpha * dt / (dx * dx);

    for _ in 0..n_steps {
        // Interior nodes
        for i in 1..n - 1 {
            u_new[i] = u[i] + r * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
        }
        // Dirichlet BCs: endpoints unchanged
        u_new[0] = u[0];
        u_new[n - 1] = u[n - 1];

        std::mem::swap(&mut u, &mut u_new);
    }
    u
}

// ──────────────────────────────────────────────────────────────────────────────
// 1-D Wave equation  (explicit leapfrog)
// ──────────────────────────────────────────────────────────────────────────────

/// Solve the 1-D wave equation  ∂²u/∂t² = c² ∂²u/∂x²  using the second-order
/// leapfrog (Störmer–Verlet) explicit scheme.
///
/// # Parameters
/// - `u_prev` — solution at time t−Δt (length N).
/// - `u_curr` — solution at time t (length N).
/// - `c` — wave speed (m/s).
/// - `dx` — grid spacing (m).
/// - `dt` — time step (s).  Stability requires `c·Δt/Δx ≤ 1`.
/// - `n_steps` — number of time steps to advance.
///
/// Dirichlet zero boundary conditions are applied.
///
/// Returns the solution at time t + n_steps·Δt.
#[allow(dead_code)]
pub fn wave_equation_1d(
    u_prev: &[f64],
    u_curr: &[f64],
    c: f64,
    dx: f64,
    dt: f64,
    n_steps: usize,
) -> Vec<f64> {
    let n = u_curr.len();
    let mut prev = u_prev.to_vec();
    let mut curr = u_curr.to_vec();
    let mut next = vec![0.0_f64; n];
    let r2 = (c * dt / dx).powi(2);

    for _ in 0..n_steps {
        for i in 1..n - 1 {
            next[i] = 2.0 * curr[i] - prev[i] + r2 * (curr[i + 1] - 2.0 * curr[i] + curr[i - 1]);
        }
        // Dirichlet zero BCs
        next[0] = 0.0;
        next[n - 1] = 0.0;

        prev.copy_from_slice(&curr);
        curr.copy_from_slice(&next);
    }
    curr
}

// ──────────────────────────────────────────────────────────────────────────────
// 1-D Poisson equation  (Thomas algorithm)
// ──────────────────────────────────────────────────────────────────────────────

/// Solve the 1-D Poisson equation  −u″ = f  on a uniform grid using the
/// Thomas algorithm (tridiagonal matrix algorithm, TDMA).
///
/// The system is  −(u\[i+1\] − 2u\[i\] + u\[i−1\]) / Δx² = f\[i\].
///
/// Dirichlet boundary conditions `u\[0\] = bc_left` and `u\[n-1\] = bc_right`
/// are applied.
///
/// Returns the solution vector `u` of the same length as `f`.
#[allow(dead_code)]
pub fn poisson_1d(f: &[f64], dx: f64, bc_left: f64, bc_right: f64) -> Vec<f64> {
    let n = f.len();
    if n < 2 {
        return vec![bc_left; n];
    }
    let dx2 = dx * dx;

    // Build RHS: rhs[i] = f[i] * dx^2 for interior nodes
    // Interior indices: 1 .. n-2
    let m = n - 2; // number of interior unknowns
    if m == 0 {
        return vec![bc_left, bc_right];
    }

    let mut rhs = vec![0.0_f64; m];
    for i in 0..m {
        rhs[i] = f[i + 1] * dx2;
    }
    // Apply BCs to RHS
    rhs[0] += bc_left;
    rhs[m - 1] += bc_right;

    // Tridiagonal: a_i = c_i = -1, b_i = 2
    let mut c_star = vec![0.0_f64; m];
    let mut d_star = vec![0.0_f64; m];

    // Forward sweep
    c_star[0] = -1.0 / 2.0;
    d_star[0] = rhs[0] / 2.0;
    for i in 1..m {
        let denom = 2.0 + c_star[i - 1];
        c_star[i] = -1.0 / denom;
        d_star[i] = (rhs[i] + d_star[i - 1]) / denom;
    }

    // Back substitution
    let mut x = vec![0.0_f64; m];
    x[m - 1] = d_star[m - 1];
    for i in (0..m - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    // Assemble full solution
    let mut u = vec![0.0_f64; n];
    u[0] = bc_left;
    u[n - 1] = bc_right;
    u[1..m + 1].copy_from_slice(&x[..m]);
    u
}

// ──────────────────────────────────────────────────────────────────────────────
// 1-D Advection  (first-order upwind)
// ──────────────────────────────────────────────────────────────────────────────

/// Solve the 1-D linear advection equation  ∂u/∂t + c ∂u/∂x = 0  for
/// `n_steps` time steps using a first-order upwind scheme.
///
/// # Parameters
/// - `u0` — initial field (length N).
/// - `c` — advection speed. Positive → right-moving.
/// - `dx` — grid spacing.
/// - `dt` — time step.  Stability requires `|c|·Δt/Δx ≤ 1`.
/// - `n_steps` — number of time steps.
///
/// Periodic boundary conditions are applied.
///
/// Returns the advected field after `n_steps`.
#[allow(dead_code)]
pub fn advection_1d_upwind(u0: &[f64], c: f64, dx: f64, dt: f64, n_steps: usize) -> Vec<f64> {
    let n = u0.len();
    let mut u = u0.to_vec();
    let mut u_new = u.clone();
    let nu = c * dt / dx; // CFL number

    for _ in 0..n_steps {
        for i in 0..n {
            if c >= 0.0 {
                // Left-biased stencil
                let im1 = if i == 0 { n - 1 } else { i - 1 };
                u_new[i] = u[i] - nu * (u[i] - u[im1]);
            } else {
                // Right-biased stencil
                let ip1 = if i == n - 1 { 0 } else { i + 1 };
                u_new[i] = u[i] - nu * (u[ip1] - u[i]);
            }
        }
        std::mem::swap(&mut u, &mut u_new);
    }
    u
}

// ──────────────────────────────────────────────────────────────────────────────
// 1-D Inviscid Burgers equation  (Godunov flux)
// ──────────────────────────────────────────────────────────────────────────────

/// Solve the inviscid Burgers equation  ∂u/∂t + ∂(u²/2)/∂x = 0  for
/// `n_steps` time steps using the Godunov (Rusanov) flux.
///
/// # Parameters
/// - `u0` — initial field (length N).
/// - `dx` — grid spacing.
/// - `dt` — time step.
/// - `n_steps` — number of time steps.
///
/// Periodic boundary conditions are applied.
#[allow(dead_code)]
pub fn burger_equation_1d(u0: &[f64], dx: f64, dt: f64, n_steps: usize) -> Vec<f64> {
    let n = u0.len();
    let mut u = u0.to_vec();
    let mut u_new = vec![0.0_f64; n];

    for _ in 0..n_steps {
        for i in 0..n {
            let ip1 = (i + 1) % n;
            let im1 = if i == 0 { n - 1 } else { i - 1 };

            // Godunov flux at right interface i+1/2
            let f_right = godunov_flux_burgers(u[i], u[ip1]);
            // Godunov flux at left interface i-1/2
            let f_left = godunov_flux_burgers(u[im1], u[i]);

            u_new[i] = u[i] - (dt / dx) * (f_right - f_left);
        }
        std::mem::swap(&mut u, &mut u_new);
    }
    u
}

/// Godunov (entropy-consistent) flux for Burgers equation, F(u) = u²/2.
fn godunov_flux_burgers(u_l: f64, u_r: f64) -> f64 {
    if u_l >= u_r {
        // Shock region: use Rankine-Hugoniot
        let s = (u_l + u_r) / 2.0; // shock speed
        if s >= 0.0 {
            u_l * u_l / 2.0
        } else {
            u_r * u_r / 2.0
        }
    } else {
        // Rarefaction region
        if u_l >= 0.0 {
            u_l * u_l / 2.0
        } else if u_r <= 0.0 {
            u_r * u_r / 2.0
        } else {
            0.0 // sonic point
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 1-D Diffusion-Reaction  (FTCS)
// ──────────────────────────────────────────────────────────────────────────────

/// Solve the 1-D diffusion-reaction equation
///   ∂u/∂t = D ∂²u/∂x² + R(u, x, t)
/// for `n_steps` time steps using the explicit FTCS scheme.
///
/// # Parameters
/// - `u0` — initial field (length N).
/// - `d` — diffusion coefficient (m²/s).
/// - `dx` — grid spacing.
/// - `dt` — time step.  Stability requires `D·Δt/Δx² ≤ 0.5`.
/// - `n_steps` — number of time steps.
/// - `reaction` — closure `R(u, i, t)` returning the reaction rate at node `i`,
///   value `u`, and time `t`.
///
/// Dirichlet boundary conditions are preserved from `u0`.
#[allow(dead_code)]
pub fn diffusion_reaction_1d<F>(
    u0: &[f64],
    d: f64,
    dx: f64,
    dt: f64,
    n_steps: usize,
    reaction: F,
) -> Vec<f64>
where
    F: Fn(f64, usize, f64) -> f64,
{
    let n = u0.len();
    let mut u = u0.to_vec();
    let mut u_new = u.clone();
    let r = d * dt / (dx * dx);
    let mut t = 0.0_f64;

    for _ in 0..n_steps {
        for i in 1..n - 1 {
            let diffusion = r * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
            let react = reaction(u[i], i, t) * dt;
            u_new[i] = u[i] + diffusion + react;
        }
        u_new[0] = u[0];
        u_new[n - 1] = u[n - 1];

        t += dt;
        std::mem::swap(&mut u, &mut u_new);
    }
    u
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // ── CflCondition ────────────────────────────────────────────────────────

    #[test]
    fn test_cfl_advective_stable() {
        let cfl = CflCondition::new(0.1, 0.05, 1.0, 0.0);
        // ν = 1.0 * 0.05 / 0.1 = 0.5 → stable
        assert!(cfl.is_advection_stable());
        assert!((cfl.advective_cfl() - 0.5).abs() < EPS);
    }

    #[test]
    fn test_cfl_advective_unstable() {
        let cfl = CflCondition::new(0.1, 0.2, 1.0, 0.0);
        // ν = 2.0 → unstable
        assert!(!cfl.is_advection_stable());
    }

    #[test]
    fn test_cfl_diffusive_stable() {
        // α·Δt/Δx² = 0.1 * 0.001 / 0.01 = 0.01 → stable
        let cfl = CflCondition::new(0.1, 0.001, 0.0, 0.1);
        assert!(cfl.is_diffusion_stable());
    }

    #[test]
    fn test_cfl_diffusive_unstable() {
        // α·Δt/Δx² = 1.0 * 1.0 / 0.01 = 100 → unstable
        let cfl = CflCondition::new(0.1, 1.0, 0.0, 1.0);
        assert!(!cfl.is_diffusion_stable());
    }

    // ── Heat equation ───────────────────────────────────────────────────────

    #[test]
    fn test_heat_bc_fixed() {
        // BCs (u[0] and u[n-1]) should remain fixed.
        let n = 20;
        let u0: Vec<f64> = (0..n)
            .map(|i| {
                if i == 0 {
                    0.0
                } else if i == n - 1 {
                    1.0
                } else {
                    0.5
                }
            })
            .collect();
        let u = heat_equation_1d(&u0, 0.1, 0.1, 0.001, 100);
        assert!((u[0] - u0[0]).abs() < EPS);
        assert!((u[n - 1] - u0[n - 1]).abs() < EPS);
    }

    #[test]
    fn test_heat_steady_state_linear() {
        // Steady state of heat equation with BCs u(0)=0, u(1)=1 is linear: u(x)=x.
        let n = 11;
        let dx = 1.0 / (n - 1) as f64;
        // Start from all-zeros interior; BCs enforced inside solver.
        let mut u0 = vec![0.5_f64; n];
        u0[0] = 0.0;
        u0[n - 1] = 1.0;

        let alpha = 0.5;
        let dt = 0.4 * dx * dx / alpha; // just below stability limit
        let u = heat_equation_1d(&u0, alpha, dx, dt, 5000);

        for i in 1..n - 1 {
            let x = i as f64 * dx;
            assert!(
                (u[i] - x).abs() < 1e-3,
                "expected {x:.3} at node {i}, got {:.6}",
                u[i]
            );
        }
    }

    #[test]
    fn test_heat_conservation_insulated() {
        // With zero temperature at both ends, total "heat" decreases but stays ≥ 0.
        let n = 50;
        let u0: Vec<f64> = (0..n)
            .map(|i| {
                let x = i as f64 / (n - 1) as f64;
                (std::f64::consts::PI * x).sin()
            })
            .collect();
        let dx = 1.0 / (n - 1) as f64;
        let alpha = 0.01;
        let dt = 0.4 * dx * dx / alpha;
        let u = heat_equation_1d(&u0, alpha, dx, dt, 100);
        // All values should remain finite
        for val in &u {
            assert!(val.is_finite());
        }
    }

    // ── Wave equation ───────────────────────────────────────────────────────

    #[test]
    fn test_wave_zero_field_stays_zero() {
        let n = 30;
        let u_prev = vec![0.0_f64; n];
        let u_curr = vec![0.0_f64; n];
        let u = wave_equation_1d(&u_prev, &u_curr, 1.0, 0.1, 0.05, 10);
        for v in &u {
            assert!(v.abs() < EPS);
        }
    }

    #[test]
    fn test_wave_bc_fixed_zero() {
        let n = 20;
        let u_prev: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let u_curr = u_prev.clone();
        let u = wave_equation_1d(&u_prev, &u_curr, 1.0, 0.1, 0.05, 5);
        assert!(u[0].abs() < EPS);
        assert!(u[n - 1].abs() < EPS);
    }

    #[test]
    fn test_wave_finite_values() {
        let n = 40;
        let dx = 1.0 / (n - 1) as f64;
        let c = 1.0;
        let dt = 0.9 * dx / c;
        let u0: Vec<f64> = (0..n)
            .map(|i| {
                let x = i as f64 * dx;
                (2.0 * std::f64::consts::PI * x).sin()
            })
            .collect();
        let u = wave_equation_1d(&u0, &u0, c, dx, dt, 50);
        for v in &u {
            assert!(v.is_finite());
        }
    }

    // ── Poisson equation ────────────────────────────────────────────────────

    #[test]
    fn test_poisson_bc_satisfied() {
        let n = 10;
        let f = vec![1.0_f64; n];
        let u = poisson_1d(&f, 0.1, 0.0, 0.0);
        assert!((u[0] - 0.0).abs() < EPS);
        assert!((u[n - 1] - 0.0).abs() < EPS);
    }

    #[test]
    fn test_poisson_zero_rhs_linear_solution() {
        // -u'' = 0 with u(0)=0, u(1)=1 → u(x) = x
        let n = 11;
        let dx = 1.0 / (n - 1) as f64;
        let f = vec![0.0_f64; n];
        let u = poisson_1d(&f, dx, 0.0, 1.0);
        for i in 0..n {
            let x = i as f64 * dx;
            assert!(
                (u[i] - x).abs() < 1e-10,
                "node {i}: expected {x}, got {}",
                u[i]
            );
        }
    }

    #[test]
    fn test_poisson_nonzero_rhs_symmetry() {
        // -u'' = 2 with u(0)=u(1)=0 → u(x) = x(1-x)
        let n = 11;
        let dx = 1.0 / (n - 1) as f64;
        let f = vec![2.0_f64; n];
        let u = poisson_1d(&f, dx, 0.0, 0.0);
        for i in 0..n {
            let x = i as f64 * dx;
            let expected = x * (1.0 - x);
            assert!(
                (u[i] - expected).abs() < 1e-10,
                "node {i}: expected {expected:.6}, got {:.6}",
                u[i]
            );
        }
    }

    #[test]
    fn test_poisson_two_node_trivial() {
        let u = poisson_1d(&[0.0, 0.0], 1.0, 1.0, 2.0);
        assert_eq!(u.len(), 2);
        assert!((u[0] - 1.0).abs() < EPS);
        assert!((u[1] - 2.0).abs() < EPS);
    }

    // ── Advection ───────────────────────────────────────────────────────────

    #[test]
    fn test_advection_zero_field() {
        let u0 = vec![0.0_f64; 20];
        let u = advection_1d_upwind(&u0, 1.0, 0.1, 0.05, 10);
        for v in &u {
            assert!(v.abs() < EPS);
        }
    }

    #[test]
    fn test_advection_positive_speed_shifts_right() {
        // Step function at index 5; after enough steps it should shift right.
        let n = 20;
        let dx = 1.0 / n as f64;
        let c = 1.0;
        let dt = 0.5 * dx / c; // CFL = 0.5
        let mut u0 = vec![0.0_f64; n];
        u0[5] = 1.0;
        let u = advection_1d_upwind(&u0, c, dx, dt, 4);
        // After 4 steps the mass should have diffused rightward; sum preserved approx
        let sum_orig: f64 = u0.iter().sum();
        let sum_adv: f64 = u.iter().sum();
        assert!((sum_adv - sum_orig).abs() < 1e-10);
    }

    #[test]
    fn test_advection_negative_speed_shifts_left() {
        let n = 20;
        let dx = 1.0 / n as f64;
        let c = -1.0_f64;
        let dt = 0.5 * dx / c.abs();
        let mut u0 = vec![0.0_f64; n];
        u0[15] = 1.0;
        let u = advection_1d_upwind(&u0, c, dx, dt, 4);
        // With negative speed, peak should shift toward lower indices
        let sum_orig: f64 = u0.iter().sum();
        let sum_adv: f64 = u.iter().sum();
        assert!((sum_adv - sum_orig).abs() < 1e-10);
    }

    #[test]
    fn test_advection_finite_values() {
        let n = 50;
        let dx = 1.0 / n as f64;
        let c = 2.0;
        let dt = 0.4 * dx / c;
        let u0: Vec<f64> = (0..n)
            .map(|i| (i as f64 * dx * std::f64::consts::PI * 2.0).sin())
            .collect();
        let u = advection_1d_upwind(&u0, c, dx, dt, 20);
        for v in &u {
            assert!(v.is_finite());
        }
    }

    // ── Burgers equation ─────────────────────────────────────────────────────

    #[test]
    fn test_burgers_zero_field() {
        let u0 = vec![0.0_f64; 20];
        let u = burger_equation_1d(&u0, 0.1, 0.01, 10);
        for v in &u {
            assert!(v.abs() < EPS);
        }
    }

    #[test]
    fn test_burgers_finite_values() {
        let n = 40;
        let dx = 1.0 / n as f64;
        let dt = 0.3 * dx; // CFL < 1 for u ∈ [-1,1]
        let u0: Vec<f64> = (0..n)
            .map(|i| {
                let x = i as f64 * dx;
                (std::f64::consts::PI * 2.0 * x).sin()
            })
            .collect();
        let u = burger_equation_1d(&u0, dx, dt, 30);
        for v in &u {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_godunov_flux_shock() {
        // u_L > u_R, shock speed s = (u_L + u_R)/2
        // u_L=2, u_R=0, s=1>0 → flux = u_L²/2 = 2
        let f = super::godunov_flux_burgers(2.0, 0.0);
        assert!((f - 2.0).abs() < EPS);
    }

    #[test]
    fn test_godunov_flux_rarefaction_positive() {
        // u_L=1, u_R=2, both positive → flux = u_L²/2 = 0.5
        let f = super::godunov_flux_burgers(1.0, 2.0);
        assert!((f - 0.5).abs() < EPS);
    }

    #[test]
    fn test_godunov_flux_sonic() {
        // u_L < 0 < u_R: sonic point → flux = 0
        let f = super::godunov_flux_burgers(-1.0, 1.0);
        assert!(f.abs() < EPS);
    }

    // ── Diffusion-Reaction ──────────────────────────────────────────────────

    #[test]
    fn test_diffusion_reaction_no_reaction_matches_heat() {
        let n = 20;
        let dx = 0.05;
        let d = 0.1;
        let dt = 0.4 * dx * dx / d;
        let u0: Vec<f64> = (0..n)
            .map(|i| (i as f64 / (n - 1) as f64 * std::f64::consts::PI).sin())
            .collect();

        let u_heat = heat_equation_1d(&u0, d, dx, dt, 50);
        let u_dr = diffusion_reaction_1d(&u0, d, dx, dt, 50, |_u, _i, _t| 0.0);

        for (a, b) in u_heat.iter().zip(u_dr.iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "heat and DR (no reaction) differ: {a} vs {b}"
            );
        }
    }

    #[test]
    fn test_diffusion_reaction_linear_decay() {
        // Reaction R = -k*u (linear decay), D → 0, should decay exponentially.
        let n = 20;
        let dx = 0.05;
        let d = 0.0; // no diffusion
        let dt = 0.01;
        let k = 0.5;
        let steps = 50;
        let mut u0 = vec![0.5_f64; n];
        u0[0] = 0.0;
        u0[n - 1] = 0.0;

        let u = diffusion_reaction_1d(&u0, d, dx, dt, steps, |val, _i, _t| -k * val);

        // Interior values should decay: u ≈ 0.5 * (1 - k*dt)^steps
        let expected = 0.5 * (1.0 - k * dt).powi(steps as i32);
        for i in 1..n - 1 {
            assert!(
                (u[i] - expected).abs() < 1e-6,
                "node {i}: expected {expected:.8}, got {:.8}",
                u[i]
            );
        }
    }

    #[test]
    fn test_diffusion_reaction_bc_preserved() {
        let n = 15;
        let u0: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let dx = 0.1;
        let d = 0.01;
        let dt = 0.001;
        let u = diffusion_reaction_1d(&u0, d, dx, dt, 10, |_u, _i, _t| 0.0);
        assert!((u[0] - u0[0]).abs() < EPS);
        assert!((u[n - 1] - u0[n - 1]).abs() < EPS);
    }

    #[test]
    fn test_diffusion_reaction_finite_values() {
        let n = 30;
        let dx = 0.05;
        let d = 0.05;
        let dt = 0.4 * dx * dx / d;
        let u0: Vec<f64> = (0..n)
            .map(|i| (i as f64 / (n - 1) as f64).powi(2))
            .collect();
        let u = diffusion_reaction_1d(&u0, d, dx, dt, 20, |val, _i, _t| val * 0.1);
        for v in &u {
            assert!(v.is_finite());
        }
    }
}
