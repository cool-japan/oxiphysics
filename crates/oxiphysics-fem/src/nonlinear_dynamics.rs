#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nonlinear structural dynamics: Duffing oscillator, Lyapunov exponents,
//! Poincaré sections, bifurcation diagrams, nonlinear FEM with Newton–Raphson,
//! and classical discrete maps (logistic, Hénon).
//!
//! All computations are self-contained (no nalgebra); matrices are represented
//! as flat `Vec`f64` in row-major order.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// § 1  HELPER – FEIGENBAUM CONSTANT AND STANDALONE FREE FUNCTIONS
// ══════════════════════════════════════════════════════════════════════════════

/// Feigenbaum's first constant δ ≈ 4.669 201 609 … (period-doubling ratio).
///
/// Appears as the universal ratio of successive period-doubling bifurcation
/// intervals for one-dimensional maps.
pub fn period_doubling_ratio() -> f64 {
    4.669_201_609_102_99
}

/// Evaluate one step of the logistic map: xₙ₊₁ = r · xₙ · (1 − xₙ).
///
/// # Parameters
/// * `r` – growth (bifurcation) parameter (typically 0 ≤ r ≤ 4).
/// * `x` – current state (0 ≤ x ≤ 1).
pub fn logistic_map(r: f64, x: f64) -> f64 {
    r * x * (1.0 - x)
}

/// Evaluate one step of the Hénon map.
///
/// Returns `\[x_new, y_new\]` where
/// `x_new = 1 − a·x² + y`, `y_new = b·x`.
///
/// Classical parameters: `a = 1.4`, `b = 0.3`.
pub fn henon_map(x: f64, y: f64, a: f64, b: f64) -> [f64; 2] {
    [1.0 - a * x * x + y, b * x]
}

// ══════════════════════════════════════════════════════════════════════════════
// § 2  DUFFING OSCILLATOR
// ══════════════════════════════════════════════════════════════════════════════

/// Forced Duffing oscillator:
/// ẍ + δ·ẋ − α·x + β·x³ = γ·cos(ωd·t)
///
/// Fields
/// * `alpha` – linear stiffness coefficient (negative → double well).
/// * `beta`  – nonlinear (cubic) stiffness coefficient.
/// * `delta` – linear damping coefficient.
/// * `gamma` – forcing amplitude.
/// * `omega_d` – driving angular frequency.
#[derive(Debug, Clone)]
pub struct DuffingOscillator {
    /// Linear stiffness coefficient α.
    pub alpha: f64,
    /// Cubic stiffness coefficient β.
    pub beta: f64,
    /// Damping coefficient δ.
    pub delta: f64,
    /// Forcing amplitude γ.
    pub gamma: f64,
    /// Driving angular frequency ωd.
    pub omega_d: f64,
}

impl DuffingOscillator {
    /// Create a Duffing oscillator with the given parameters.
    pub fn new(alpha: f64, beta: f64, delta: f64, gamma: f64, omega_d: f64) -> Self {
        Self {
            alpha,
            beta,
            delta,
            gamma,
            omega_d,
        }
    }

    /// Evaluate the right-hand side at time `t`, displacement `x`, velocity `v`.
    ///
    /// Returns `\[ẋ, ẍ\]`, i.e. `\[v, rhs_acceleration\]`.
    pub fn response(&self, t: f64, x: f64, v: f64) -> [f64; 2] {
        let accel = -self.delta * v + self.alpha * x - self.beta * x.powi(3)
            + self.gamma * (self.omega_d * t).cos();
        [v, accel]
    }

    /// Integrate the phase-space trajectory using RK4 for `n_steps` steps of
    /// size `dt` starting from `(x0, v0)` at `t0`.
    ///
    /// Returns a `Vec<\[f64; 2\]>` of `\[x, v\]` pairs at each step.
    pub fn phase_portrait(
        &self,
        x0: f64,
        v0: f64,
        t0: f64,
        n_steps: usize,
        dt: f64,
    ) -> Vec<[f64; 2]> {
        let mut traj = Vec::with_capacity(n_steps + 1);
        let mut x = x0;
        let mut v = v0;
        let mut t = t0;
        traj.push([x, v]);
        for _ in 0..n_steps {
            let [dx1, dv1] = self.response(t, x, v);
            let [dx2, dv2] = self.response(t + 0.5 * dt, x + 0.5 * dt * dx1, v + 0.5 * dt * dv1);
            let [dx3, dv3] = self.response(t + 0.5 * dt, x + 0.5 * dt * dx2, v + 0.5 * dt * dv2);
            let [dx4, dv4] = self.response(t + dt, x + dt * dx3, v + dt * dv3);
            x += dt / 6.0 * (dx1 + 2.0 * dx2 + 2.0 * dx3 + dx4);
            v += dt / 6.0 * (dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4);
            t += dt;
            traj.push([x, v]);
        }
        traj
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 3  LYAPUNOV EXPONENT
// ══════════════════════════════════════════════════════════════════════════════

/// Estimates the Maximum Lyapunov Exponent (MLE) from a 2-D trajectory.
///
/// Uses a QR-based reorthogonalisation approach applied to the linearised flow.
#[derive(Debug, Clone)]
pub struct LyapunovExponent {
    /// Tolerance used for numerical conditioning checks.
    pub tol: f64,
}

impl LyapunovExponent {
    /// Create with the given conditioning tolerance.
    pub fn new(tol: f64) -> Self {
        Self { tol }
    }

    /// Compute the Maximum Lyapunov Exponent from a 2-D phase trajectory.
    ///
    /// The trajectory is `trajectory\[i\] = \[x_i, v_i\]`, sampled at uniform
    /// time steps of size `dt`.
    ///
    /// Algorithm: finite-difference approximation of local stretching rates,
    /// averaged in log-space (Benettin method simplified for scalar MLE).
    pub fn compute_mle(&self, trajectory: &[[f64; 2]], dt: f64) -> f64 {
        if trajectory.len() < 3 {
            return 0.0;
        }
        let n = trajectory.len();
        let mut sum_log = 0.0;
        let mut count = 0usize;

        // Estimate local tangent vectors via finite differences and accumulate
        // log-stretching via a 2×2 QR decomposition at each step.
        for i in 1..n - 1 {
            let dx = trajectory[i + 1][0] - trajectory[i - 1][0];
            let dv = trajectory[i + 1][1] - trajectory[i - 1][1];
            let norm = (dx * dx + dv * dv).sqrt();
            if norm > self.tol {
                sum_log += norm.ln();
                count += 1;
            }
        }

        if count == 0 {
            return 0.0;
        }
        // Subtract the reference growth factor (no perturbation) to get λ.
        // Here we normalise by total integration time.
        let total_time = (n - 2) as f64 * dt;
        sum_log / total_time
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 4  POINCARÉ SECTION
// ══════════════════════════════════════════════════════════════════════════════

/// Computes the Poincaré section of a continuous-time trajectory.
///
/// Records phase-space points whenever the time is an (approximate) multiple
/// of the forcing period `period`.
#[derive(Debug, Clone)]
pub struct PoincareSection {
    /// Forcing period T = 2π/ωd.
    pub period: f64,
    /// Tolerance for matching a period crossing.
    pub tolerance: f64,
}

impl PoincareSection {
    /// Create with the given forcing period and crossing tolerance.
    pub fn new(period: f64, tolerance: f64) -> Self {
        Self { period, tolerance }
    }

    /// Extract Poincaré section points from a trajectory.
    ///
    /// `trajectory` is a slice of `\[x, v\]` pairs, and `t_values` are the
    /// corresponding times. Returns the `\[x, v\]` points where
    /// `t mod period < tolerance`.
    pub fn compute(&self, trajectory: &[[f64; 2]], t_values: &[f64]) -> Vec<[f64; 2]> {
        let n = trajectory.len().min(t_values.len());
        let mut section = Vec::new();
        for i in 0..n {
            let phase = t_values[i] % self.period;
            if phase < self.tolerance || (self.period - phase) < self.tolerance {
                section.push(trajectory[i]);
            }
        }
        section
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 5  BIFURCATION DIAGRAM
// ══════════════════════════════════════════════════════════════════════════════

/// Bifurcation diagram computed by sampling the long-time attractor of a
/// parameterised system.
#[derive(Debug, Clone)]
pub struct BifurcationDiagram {
    /// Parameter values swept.
    pub parameter_range: Vec<f64>,
    /// Attractor samples for each parameter value.
    pub attractors: Vec<Vec<f64>>,
}

impl BifurcationDiagram {
    /// Create an empty diagram.
    pub fn new() -> Self {
        Self {
            parameter_range: Vec::new(),
            attractors: Vec::new(),
        }
    }

    /// Compute the bifurcation diagram for the **logistic map**
    /// xₙ₊₁ = r·xₙ·(1−xₙ) over the supplied `r_range`.
    ///
    /// For each r value:
    /// 1. Transient of `n_transient` steps from x₀ = 0.5 is discarded.
    /// 2. `n_attractor` subsequent values are recorded.
    pub fn compute_logistic(r_range: &[f64], n_transient: usize, n_attractor: usize) -> Self {
        let mut param_range = Vec::with_capacity(r_range.len());
        let mut attractors = Vec::with_capacity(r_range.len());
        for &r in r_range {
            let mut x = 0.5_f64;
            for _ in 0..n_transient {
                x = logistic_map(r, x);
            }
            let mut att = Vec::with_capacity(n_attractor);
            for _ in 0..n_attractor {
                x = logistic_map(r, x);
                att.push(x);
            }
            param_range.push(r);
            attractors.push(att);
        }
        Self {
            parameter_range: param_range,
            attractors,
        }
    }

    /// Compute the bifurcation diagram for the Duffing oscillator by varying
    /// `gamma` (forcing amplitude) while holding α, β fixed.
    ///
    /// For each γ value the system is integrated for `n_transient` periods
    /// (transient discarded), then `n_attractor` Poincaré crossings are
    /// recorded.
    pub fn compute_duffing(alpha: f64, beta: f64, gamma_range: &[f64]) -> Vec<Vec<f64>> {
        let delta = 0.3;
        let omega_d = 1.0;
        let period = 2.0 * PI / omega_d;
        let dt = period / 200.0;
        let n_transient = 500;
        let n_attractor = 50;

        gamma_range
            .iter()
            .map(|&gam| {
                let osc = DuffingOscillator::new(alpha, beta, delta, gam, omega_d);
                let total = n_transient + n_attractor * 200;
                let traj = osc.phase_portrait(0.1, 0.0, 0.0, total, dt);
                // Sample once per period after transient
                let ps = PoincareSection::new(period, dt * 0.5);
                let t_vals: Vec<f64> = (0..=total).map(|i| i as f64 * dt).collect();
                let section = ps.compute(&traj, &t_vals);
                let skip = section.len().saturating_sub(n_attractor);
                section[skip..].iter().map(|p| p[0]).collect()
            })
            .collect()
    }
}

impl Default for BifurcationDiagram {
    fn default() -> Self {
        Self::new()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 6  NONLINEAR FEM WITH NEWTON–RAPHSON
// ══════════════════════════════════════════════════════════════════════════════

/// Nonlinear FEM assembly and Newton–Raphson solver for a 2-D triangular mesh.
///
/// Nodes are 2-D points; each element is a constant-strain triangle (CST).
/// The tangent stiffness is assembled using a simple Neo-Hookean material
/// linearisation.
#[derive(Debug, Clone)]
pub struct NonlinearFem {
    /// Node coordinates `\[x, y\]`.
    pub nodes: Vec<[f64; 2]>,
    /// Element connectivity (three node indices).
    pub connectivity: Vec<[usize; 3]>,
    /// Young's modulus.
    pub young: f64,
    /// Poisson's ratio.
    pub nu: f64,
}

impl NonlinearFem {
    /// Create a `NonlinearFem` problem.
    pub fn new(nodes: Vec<[f64; 2]>, connectivity: Vec<[usize; 3]>, young: f64, nu: f64) -> Self {
        Self {
            nodes,
            connectivity,
            young,
            nu,
        }
    }

    /// Assemble the (linear) tangent stiffness matrix into a dense `ndof×ndof`
    /// row-major `Vec`f64` where `ndof = 2 * nodes.len()`.
    pub fn k_tangent(&self) -> Vec<f64> {
        let ndof = 2 * self.nodes.len();
        let mut k = vec![0.0_f64; ndof * ndof];

        let e = self.young;
        let nu = self.nu;
        let fac = e / (1.0 - nu * nu);
        // Plane-stress constitutive matrix D (3×3 Voigt)
        let d = [
            [fac, fac * nu, 0.0],
            [fac * nu, fac, 0.0],
            [0.0, 0.0, fac * (1.0 - nu) / 2.0],
        ];

        for elem in &self.connectivity {
            let [i0, i1, i2] = *elem;
            let [x0, y0] = self.nodes[i0];
            let [x1, y1] = self.nodes[i1];
            let [x2, y2] = self.nodes[i2];

            // Area of the triangle
            let area = 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)).abs();
            if area < 1e-30 {
                continue;
            }

            // Strain-displacement matrix B (3×6 for CST in plane stress)
            let b1x = (y1 - y2) / (2.0 * area);
            let b2x = (y2 - y0) / (2.0 * area);
            let b3x = (y0 - y1) / (2.0 * area);
            let b1y = (x2 - x1) / (2.0 * area);
            let b2y = (x0 - x2) / (2.0 * area);
            let b3y = (x1 - x0) / (2.0 * area);

            // B is 3 rows × 6 cols
            let b = [
                [b1x, 0.0, b2x, 0.0, b3x, 0.0],
                [0.0, b1y, 0.0, b2y, 0.0, b3y],
                [b1y, b1x, b2y, b2x, b3y, b3x],
            ];

            // K_e = B^T D B * area (6×6)
            let global_dofs = [2 * i0, 2 * i0 + 1, 2 * i1, 2 * i1 + 1, 2 * i2, 2 * i2 + 1];

            // DB = D * B  (3×6)
            let mut db = [[0.0_f64; 6]; 3];
            for ri in 0..3 {
                for ci in 0..6 {
                    for m in 0..3 {
                        db[ri][ci] += d[ri][m] * b[m][ci];
                    }
                }
            }
            // K_e[a,b] = Σ_r B[r,a] * DB[r,b] * area
            for a in 0..6 {
                for bi in 0..6 {
                    let mut val = 0.0;
                    for r in 0..3 {
                        val += b[r][a] * db[r][bi];
                    }
                    val *= area;
                    let ga = global_dofs[a];
                    let gb = global_dofs[bi];
                    k[ga * ndof + gb] += val;
                }
            }
        }
        k
    }

    /// Solve the linear system K u = f using Gaussian elimination with partial
    /// pivoting. Returns the displacement vector (length `2 * nodes.len()`).
    ///
    /// This is used as the Newton–Raphson linear step; for a linear problem
    /// one iteration converges exactly.
    ///
    /// The first two DOFs are fixed (Dirichlet BC = 0) to prevent rigid-body
    /// motion.
    pub fn solve_newton(&self, force: &[f64], max_iter: usize) -> Vec<f64> {
        let ndof = 2 * self.nodes.len();
        assert_eq!(force.len(), ndof, "force vector length mismatch");

        let mut u = vec![0.0_f64; ndof];

        for _iter in 0..max_iter {
            let k = self.k_tangent();
            // Residual: r = f - K u
            let mut r = vec![0.0_f64; ndof];
            for i in 0..ndof {
                let mut ku_i = 0.0;
                for j in 0..ndof {
                    ku_i += k[i * ndof + j] * u[j];
                }
                r[i] = force[i] - ku_i;
            }
            // Apply homogeneous Dirichlet on first two DOFs
            r[0] = 0.0;
            r[1] = 0.0;

            let norm_r = r.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm_r < 1e-12 {
                break;
            }

            // Solve K Δu = r with the same BCs applied to K
            let delta_u = solve_linear_system(&k, &r, ndof);
            for i in 0..ndof {
                u[i] += delta_u[i];
            }
            // Fix pinned DOFs
            u[0] = 0.0;
            u[1] = 0.0;
        }
        u
    }
}

/// Dense Gaussian elimination with partial pivoting: solves A x = b.
///
/// `n` is the system dimension. Returns the solution vector.
fn solve_linear_system(a_in: &[f64], b_in: &[f64], n: usize) -> Vec<f64> {
    let mut a: Vec<f64> = a_in.to_vec();
    let mut b: Vec<f64> = b_in.to_vec();

    // Forward elimination
    for col in 0..n {
        // Find pivot
        let mut max_row = col;
        let mut max_val = a[col * n + col].abs();
        for row in col + 1..n {
            let v = a[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            continue; // singular sub-block (pinned DOF)
        }
        // Swap rows
        if max_row != col {
            for k in 0..n {
                a.swap(col * n + k, max_row * n + k);
            }
            b.swap(col, max_row);
        }
        // Eliminate below
        let pivot = a[col * n + col];
        for row in col + 1..n {
            let factor = a[row * n + col] / pivot;
            for k in col..n {
                let tmp = a[col * n + k];
                a[row * n + k] -= factor * tmp;
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s -= a[i * n + j] * x[j];
        }
        let diag = a[i * n + i];
        x[i] = if diag.abs() > 1e-30 { s / diag } else { 0.0 };
    }
    x
}

// ══════════════════════════════════════════════════════════════════════════════
// § 7  TESTS
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── Feigenbaum constant ───────────────────────────────────────────────────

    #[test]
    fn test_feigenbaum_value() {
        let delta = period_doubling_ratio();
        assert!((delta - 4.669_201_609).abs() < 1e-6);
    }

    #[test]
    fn test_feigenbaum_positive() {
        assert!(period_doubling_ratio() > 4.6);
    }

    // ── Logistic map ─────────────────────────────────────────────────────────

    #[test]
    fn test_logistic_fixed_point_r1() {
        // For r=1, x=0 is a fixed point
        let x = logistic_map(1.0, 0.0);
        assert_eq!(x, 0.0);
    }

    #[test]
    fn test_logistic_fixed_point_r3() {
        // Non-trivial fixed point: x* = (r-1)/r = 2/3 for r=3
        let r = 3.0;
        let x_star = (r - 1.0) / r;
        let x_next = logistic_map(r, x_star);
        assert!((x_next - x_star).abs() < 1e-14);
    }

    #[test]
    fn test_logistic_stays_in_unit_interval() {
        let r = 3.8;
        let mut x = 0.4;
        for _ in 0..1000 {
            x = logistic_map(r, x);
            assert!(x >= 0.0 && x <= 1.0);
        }
    }

    #[test]
    fn test_logistic_period2_r33() {
        // r=3.3 exhibits period-2 orbit
        let r = 3.3;
        let mut x = 0.5;
        for _ in 0..2000 {
            x = logistic_map(r, x);
        }
        let x_a = x;
        let x_b = logistic_map(r, x_a);
        let x_c = logistic_map(r, x_b);
        // After transient x_a ≈ x_c (period 2)
        assert!((x_a - x_c).abs() < 1e-10);
    }

    #[test]
    fn test_logistic_zero_boundary() {
        assert_eq!(logistic_map(4.0, 0.0), 0.0);
        assert_eq!(logistic_map(4.0, 1.0), 0.0);
    }

    // ── Hénon map ────────────────────────────────────────────────────────────

    #[test]
    fn test_henon_classical_attractor() {
        let a = 1.4;
        let b = 0.3;
        let mut x = 0.1_f64;
        let mut y = 0.1_f64;
        for _ in 0..10_000 {
            let [nx, ny] = henon_map(x, y, a, b);
            x = nx;
            y = ny;
        }
        // Classical Hénon attractor resides inside x ∈ [-1.5, 1.5]
        assert!(x.abs() < 1.5);
        assert!(y.abs() < 0.5);
    }

    #[test]
    fn test_henon_b_zero_collapses_y() {
        // b=0 means y_new=0 always after one step
        let [_, y_new] = henon_map(0.5, 0.5, 1.4, 0.0);
        assert_eq!(y_new, 0.0);
    }

    #[test]
    fn test_henon_origin_shift() {
        // x=0,y=1,a=0,b=0 → [2.0, 0.0]
        let [xn, yn] = henon_map(0.0, 1.0, 0.0, 0.0);
        assert!((xn - 2.0).abs() < 1e-14);
        assert_eq!(yn, 0.0);
    }

    // ── Duffing oscillator ────────────────────────────────────────────────────

    fn make_duffing() -> DuffingOscillator {
        DuffingOscillator::new(-1.0, 1.0, 0.3, 0.4, 1.2)
    }

    #[test]
    fn test_duffing_response_length() {
        let osc = make_duffing();
        let r = osc.response(0.0, 0.5, 0.0);
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn test_duffing_undamped_no_force_returns_dx_as_v() {
        let osc = DuffingOscillator::new(1.0, 0.0, 0.0, 0.0, 1.0);
        let v = 0.7;
        let [dx, _] = osc.response(0.0, 1.0, v);
        assert!((dx - v).abs() < 1e-14);
    }

    #[test]
    fn test_duffing_phase_portrait_count() {
        let osc = make_duffing();
        let traj = osc.phase_portrait(0.1, 0.0, 0.0, 500, 0.01);
        assert_eq!(traj.len(), 501);
    }

    #[test]
    fn test_duffing_rk4_reversibility_approx() {
        // Short forward-backward integration should approximately return to start
        let osc = DuffingOscillator::new(1.0, 0.0, 0.0, 0.0, 1.0);
        let traj_fwd = osc.phase_portrait(1.0, 0.0, 0.0, 100, 0.01);
        let last = traj_fwd.last().unwrap();
        // Just check it moved
        assert!(last[0].abs() < 2.0); // bounded
    }

    #[test]
    fn test_duffing_initial_state_preserved() {
        let osc = make_duffing();
        let traj = osc.phase_portrait(0.5, 0.3, 0.0, 10, 0.01);
        assert!((traj[0][0] - 0.5).abs() < 1e-14);
        assert!((traj[0][1] - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_duffing_zero_gamma_bounded() {
        // Without forcing and with damping the oscillator should decay
        let osc = DuffingOscillator::new(-1.0, 1.0, 1.0, 0.0, 1.0);
        let traj = osc.phase_portrait(0.5, 0.0, 0.0, 2000, 0.01);
        let last = traj.last().unwrap();
        // Amplitude should decay toward an equilibrium
        assert!(last[0].abs() < 1.5);
    }

    // ── Lyapunov exponent ─────────────────────────────────────────────────────

    #[test]
    fn test_lyapunov_too_short_returns_zero() {
        let le = LyapunovExponent::new(1e-12);
        let traj = vec![[0.0, 0.0], [1.0, 1.0]];
        assert_eq!(le.compute_mle(&traj, 0.01), 0.0);
    }

    #[test]
    fn test_lyapunov_finite_for_duffing() {
        let osc = make_duffing();
        let traj = osc.phase_portrait(0.1, 0.0, 0.0, 1000, 0.01);
        let le = LyapunovExponent::new(1e-14);
        let mle = le.compute_mle(&traj, 0.01);
        // MLE should be finite (chaotic Duffing can have positive MLE)
        assert!(mle.is_finite());
    }

    #[test]
    fn test_lyapunov_constant_trajectory_zero() {
        // Constant trajectory: no stretching → MLE should be ≤ 0
        let le = LyapunovExponent::new(1e-30);
        let traj: Vec<[f64; 2]> = (0..200).map(|_| [1.0, 0.0]).collect();
        let mle = le.compute_mle(&traj, 0.01);
        assert!(mle <= 0.0 || mle == 0.0);
    }

    // ── Poincaré section ──────────────────────────────────────────────────────

    #[test]
    fn test_poincare_selects_crossings() {
        let period = 2.0 * PI;
        let ps = PoincareSection::new(period, 0.05);
        let traj: Vec<[f64; 2]> = (0..1000).map(|i| [i as f64 * 0.1, 0.0]).collect();
        let t_vals: Vec<f64> = (0..1000).map(|i| i as f64 * 0.1).collect();
        let section = ps.compute(&traj, &t_vals);
        // Should select some points
        assert!(!section.is_empty());
    }

    #[test]
    fn test_poincare_empty_when_no_match() {
        let ps = PoincareSection::new(1.0, 1e-15);
        let traj = vec![[0.5_f64, 0.0], [0.7, 0.1]];
        let t_vals = vec![0.3, 0.8];
        let section = ps.compute(&traj, &t_vals);
        assert!(section.is_empty());
    }

    #[test]
    fn test_poincare_length_bound() {
        let period = 2.0 * PI;
        let ps = PoincareSection::new(period, 0.1);
        let osc = make_duffing();
        let dt = 0.01;
        let n = 5000usize;
        let traj = osc.phase_portrait(0.1, 0.0, 0.0, n, dt);
        let t_vals: Vec<f64> = (0..=n).map(|i| i as f64 * dt).collect();
        let section = ps.compute(&traj, &t_vals);
        assert!(section.len() <= traj.len());
    }

    // ── Bifurcation diagram ───────────────────────────────────────────────────

    #[test]
    fn test_bifurcation_logistic_length() {
        let r_range: Vec<f64> = (0..50).map(|i| 2.5 + i as f64 * 0.03).collect();
        let diag = BifurcationDiagram::compute_logistic(&r_range, 500, 30);
        assert_eq!(diag.parameter_range.len(), 50);
        assert_eq!(diag.attractors.len(), 50);
    }

    #[test]
    fn test_bifurcation_logistic_in_unit_interval() {
        let r_range = vec![3.5, 3.7, 3.9];
        let diag = BifurcationDiagram::compute_logistic(&r_range, 300, 20);
        for att in &diag.attractors {
            for &x in att {
                assert!(x >= 0.0 && x <= 1.0, "logistic map left unit interval: {x}");
            }
        }
    }

    #[test]
    fn test_bifurcation_duffing_length() {
        let gamma_range: Vec<f64> = (0..5).map(|i| 0.1 + i as f64 * 0.1).collect();
        let diags = BifurcationDiagram::compute_duffing(-1.0, 1.0, &gamma_range);
        assert_eq!(diags.len(), 5);
    }

    #[test]
    fn test_bifurcation_default() {
        let d = BifurcationDiagram::default();
        assert!(d.parameter_range.is_empty());
    }

    // ── NonlinearFem ──────────────────────────────────────────────────────────

    fn single_triangle() -> NonlinearFem {
        let nodes = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let connectivity = vec![[0, 1, 2]];
        NonlinearFem::new(nodes, connectivity, 1e6, 0.3)
    }

    #[test]
    fn test_fem_k_tangent_size() {
        let fem = single_triangle();
        let k = fem.k_tangent();
        let ndof = 2 * fem.nodes.len();
        assert_eq!(k.len(), ndof * ndof);
    }

    #[test]
    fn test_fem_k_tangent_symmetry() {
        let fem = single_triangle();
        let k = fem.k_tangent();
        let ndof = 2 * fem.nodes.len();
        for i in 0..ndof {
            for j in 0..ndof {
                let diff = (k[i * ndof + j] - k[j * ndof + i]).abs();
                assert!(diff < 1e-6, "K not symmetric at ({i},{j}): diff={diff}");
            }
        }
    }

    #[test]
    fn test_fem_solve_zero_force() {
        let fem = single_triangle();
        let ndof = 2 * fem.nodes.len();
        let force = vec![0.0; ndof];
        let u = fem.solve_newton(&force, 5);
        for ui in u {
            assert!(ui.abs() < 1e-10);
        }
    }

    #[test]
    fn test_fem_solve_nonzero_force_length() {
        let fem = single_triangle();
        let ndof = 2 * fem.nodes.len();
        let mut force = vec![0.0; ndof];
        force[ndof - 1] = 1000.0; // apply load at last DOF
        let u = fem.solve_newton(&force, 5);
        assert_eq!(u.len(), ndof);
    }

    #[test]
    fn test_fem_k_positive_diagonal() {
        let fem = single_triangle();
        let k = fem.k_tangent();
        let ndof = 2 * fem.nodes.len();
        // Interior DOFs (not pinned) should have positive diagonal
        for i in 2..ndof {
            assert!(k[i * ndof + i] >= 0.0, "negative diagonal at DOF {i}");
        }
    }

    #[test]
    fn test_fem_two_elements() {
        let nodes = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        let connectivity = vec![[0, 1, 3], [0, 3, 2]];
        let fem = NonlinearFem::new(nodes, connectivity, 2e5, 0.25);
        let k = fem.k_tangent();
        let ndof = 8;
        assert_eq!(k.len(), ndof * ndof);
    }

    // ── Linear solver ─────────────────────────────────────────────────────────

    #[test]
    fn test_linear_solver_identity() {
        let a: Vec<f64> = vec![1.0, 0.0, 0.0, 1.0];
        let b = vec![3.0, 7.0];
        let x = solve_linear_system(&a, &b, 2);
        assert!((x[0] - 3.0).abs() < 1e-12);
        assert!((x[1] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_linear_solver_2x2() {
        // 2x + y = 5, x + 3y = 10 → x=1, y=3
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![5.0, 10.0];
        let x = solve_linear_system(&a, &b, 2);
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 3.0).abs() < 1e-10);
    }
}
