//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{OdeSystem, rk45_step, solve_linear_system};

/// Dormand-Prince RK45 integrator (struct form).
///
/// Computes an adaptive step and returns a [`StepResult`] containing the new
/// state, error estimate, and acceptance flag.
pub struct DormandPrince45 {
    /// Absolute tolerance for step acceptance.
    pub atol: f64,
    /// Relative tolerance for step acceptance.
    pub rtol: f64,
}
impl DormandPrince45 {
    /// Create a new `DormandPrince45` integrator with given tolerances.
    pub fn new(atol: f64, rtol: f64) -> Self {
        Self { atol, rtol }
    }
    /// Perform a single step of length `dt` starting from state `y`.
    ///
    /// Returns a [`StepResult`] indicating whether the step should be accepted.
    pub fn step(&self, y: &[f64], f: impl Fn(&[f64]) -> Vec<f64>, dt: f64) -> StepResult {
        let (y_new, err) = rk45_step(y, f, dt, self.atol);
        let tol_norm = (self.atol + self.rtol * y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max))
            .max(self.atol);
        let err_norm = err / tol_norm;
        let accepted = err_norm <= 1.0 || err == 0.0;
        StepResult {
            y: y_new,
            error_estimate: err,
            accepted,
        }
    }
}
/// Basic position-Verlet integrator (Störmer's method).
///
/// Uses two consecutive positions `x` and `x_prev` to step forward:
/// `x_new = 2*x - x_prev + a * dt^2`
pub struct Verlet;
impl Verlet {
    /// Create a new `Verlet` integrator.
    pub fn new() -> Self {
        Self
    }
    /// Advance the position Verlet by one step.
    ///
    /// Updates `x` to the new position and moves `x_prev` to the old `x`.
    ///
    /// Formula: `x_new = 2*x - x_prev + a * dt²`
    pub fn step(x: &mut Vec<f64>, x_prev: &mut Vec<f64>, a: &[f64], dt: f64) {
        let x_new: Vec<f64> = x
            .iter()
            .zip(x_prev.iter())
            .zip(a.iter())
            .map(|((xi, xp), ai)| 2.0 * xi - xp + ai * dt * dt)
            .collect();
        *x_prev = x.clone();
        *x = x_new;
    }
}
/// Adaptive Dormand-Prince RK45 integrator that returns the full trajectory.
///
/// Distinct from [`DormandPrince45`] (single-step struct): this struct drives a
/// complete integration from `t0` to `tf` with automatic step-size control and
/// returns every accepted `(t, y)` pair.
pub struct DormandPrince {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Maximum allowed step size.
    pub h_max: f64,
}
impl DormandPrince {
    /// Create a new `DormandPrince` integrator.
    pub fn new(rtol: f64, atol: f64, h_max: f64) -> Self {
        Self { rtol, atol, h_max }
    }
    /// Integrate `f(t, y) -> dy/dt` from `t0` to `tf` starting at `y0`.
    ///
    /// Returns a vector of `(t, y)` pairs at every accepted step.
    pub fn integrate<F>(&self, f: F, t0: f64, tf: f64, y0: Vec<f64>) -> Vec<(f64, Vec<f64>)>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let _n = y0.len();
        let mut t = t0;
        let mut y = y0;
        let mut h = ((tf - t0) / 100.0).min(self.h_max).max(1e-12);
        let mut result = vec![(t, y.clone())];
        let safety = 0.9_f64;
        let min_scale = 0.2_f64;
        let max_scale = 10.0_f64;
        while t < tf {
            if t + h > tf {
                h = tf - t;
            }
            let f_at_t = |yv: &[f64]| f(t, yv);
            let (y_new, err) = rk45_step(&y, f_at_t, h, self.atol);
            let tol_norm = (self.atol
                + self.rtol * y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max))
            .max(self.atol);
            if err <= tol_norm || err == 0.0 {
                t += h;
                y = y_new;
                result.push((t, y.clone()));
                let factor = if err == 0.0 {
                    max_scale
                } else {
                    let en = err / tol_norm;
                    (safety * en.powf(-0.2)).clamp(min_scale, max_scale)
                };
                h = (h * factor).min(self.h_max);
            } else {
                let en = err / tol_norm;
                let factor = (safety * en.powf(-0.25)).clamp(min_scale, 1.0);
                h *= factor;
            }
            if h < 1e-15 * (tf - t0).abs().max(1.0) {
                break;
            }
        }
        result
    }
}
impl DormandPrince {
    /// Compute the DOPRI5 embedded error norm for a single step.
    ///
    /// Given the state `y` and the seven Dormand-Prince stage derivatives
    /// `k1..k7`, returns the RMS error norm scaling each component by
    /// `atol + rtol * |y_i|`.
    ///
    /// The DOPRI5 error coefficients are the difference between the 5th- and
    /// 4th-order solutions in the standard Butcher tableau.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_error_norm(
        &self,
        y: &[f64],
        k1: &[f64],
        k3: &[f64],
        k4: &[f64],
        k5: &[f64],
        k6: &[f64],
        k7: &[f64],
        h: f64,
    ) -> f64 {
        let n = y.len();
        if n == 0 {
            return 0.0;
        }
        let e1 = 71.0 / 57600.0;
        let e3 = -71.0 / 16695.0;
        let e4 = 71.0 / 1920.0;
        let e5 = -17253.0 / 339200.0;
        let e6 = 22.0 / 525.0;
        let e7 = -1.0 / 40.0;
        let sum_sq: f64 = (0..n)
            .map(|i| {
                let err_i = h
                    * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
                let sc = self.atol + self.rtol * y[i].abs();
                let r = err_i / sc;
                r * r
            })
            .sum();
        (sum_sq / n as f64).sqrt()
    }
}
/// Bulirsch-Stoer extrapolation integrator.
///
/// Uses the modified midpoint rule at increasing substep counts and polynomial
/// (Richardson) extrapolation to achieve high accuracy with few function
/// evaluations.
pub struct BulirschStoer {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Maximum extrapolation order (= number of different step sequences tried).
    pub max_order: usize,
}
impl BulirschStoer {
    /// Create a new `BulirschStoer` integrator.
    pub fn new(rtol: f64, atol: f64, max_order: usize) -> Self {
        Self {
            rtol,
            atol,
            max_order,
        }
    }
    /// Modified midpoint rule: integrate from `t` to `t + h` using `n` substeps.
    ///
    /// Returns the state at `t + h`.
    pub fn midpoint_method(
        f: &impl Fn(f64, &[f64]) -> Vec<f64>,
        t: f64,
        y: &[f64],
        h: f64,
        n: usize,
    ) -> Vec<f64> {
        let dim = y.len();
        let dt = h / n as f64;
        let mut z0 = y.to_vec();
        let k1 = f(t, &z0);
        let mut z1: Vec<f64> = (0..dim).map(|i| z0[i] + dt * k1[i]).collect();
        for m in 0..(n - 1) {
            let tm = t + (m + 1) as f64 * dt;
            let k = f(tm, &z1);
            let z2: Vec<f64> = (0..dim).map(|i| z0[i] + 2.0 * dt * k[i]).collect();
            z0 = z1;
            z1 = z2;
        }
        let k_last = f(t + h, &z1);
        (0..dim)
            .map(|i| 0.5 * (z0[i] + z1[i] + dt * k_last[i]))
            .collect()
    }
    /// Richardson extrapolation of a triangular array `t_table`.
    ///
    /// `t_table[k]` is the midpoint result using `n_k` substeps.
    /// Builds the Neville table and returns the highest-order estimate.
    pub fn extrapolate(t_table: &[Vec<f64>], n: usize) -> Vec<f64> {
        let dim = t_table[0].len();
        let mut table: Vec<Vec<f64>> = t_table[..n].to_vec();
        for j in 1..n {
            for i in (j..n).rev() {
                let ni = 2 * (i + 1);
                let nj = 2 * (i - j + 1);
                let ratio = (ni as f64 / nj as f64).powi(2);
                let prev = table[i - 1].clone();
                for d in 0..dim {
                    table[i][d] = table[i][d] + (table[i][d] - prev[d]) / (ratio - 1.0);
                }
            }
        }
        table[n - 1].clone()
    }
    /// Integrate `f` from `t0` to `tf` starting at `y0`.
    ///
    /// Returns `(t_final, y_final)` pairs at each accepted macro-step.
    pub fn integrate<F>(&self, f: &F, t0: f64, tf: f64, y0: Vec<f64>) -> Vec<(f64, Vec<f64>)>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut t = t0;
        let mut y = y0;
        let mut h = (tf - t0) / 10.0;
        let max_ord = self.max_order.max(2);
        let mut result = vec![(t, y.clone())];
        while t < tf {
            if t + h > tf {
                h = tf - t;
            }
            let mut t_table: Vec<Vec<f64>> = Vec::with_capacity(max_ord);
            for k in 1..=max_ord {
                let n_steps = 2 * k;
                let yk = Self::midpoint_method(f, t, &y, h, n_steps);
                t_table.push(yk);
            }
            let y_extrap = Self::extrapolate(&t_table, max_ord);
            let y_prev = Self::extrapolate(&t_table, max_ord.saturating_sub(1).max(1));
            let err: f64 = y_extrap
                .iter()
                .zip(y_prev.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            let tol = self.atol + self.rtol * y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            if err <= tol || max_ord <= 2 {
                t += h;
                y = y_extrap;
                result.push((t, y.clone()));
                let factor = if err == 0.0 {
                    2.0
                } else {
                    (0.9 * (tol / err).powf(1.0 / (2 * max_ord) as f64)).clamp(0.2, 5.0)
                };
                h *= factor;
            } else {
                h *= 0.5;
            }
            if h < 1e-15 * (tf - t0).abs().max(1.0) {
                break;
            }
        }
        result
    }
}
impl BulirschStoer {
    /// Richardson extrapolation from a sequence of midpoint estimates.
    ///
    /// `estimates[k]` is the midpoint result using `n_steps[k]` substeps.
    /// Applies Aitken-Neville polynomial extrapolation and returns the
    /// highest-order estimate together with an L2 error norm (difference
    /// between the last two extrapolation columns).
    ///
    /// Returns `(best_estimate, error_norm)`.
    pub fn richardson_extrapolate(estimates: &[Vec<f64>], n_steps: &[usize]) -> (Vec<f64>, f64) {
        let m = estimates.len();
        if m == 0 {
            return (vec![], 0.0);
        }
        if m == 1 {
            return (estimates[0].clone(), 0.0);
        }
        let dim = estimates[0].len();
        let mut t: Vec<Vec<f64>> = estimates.to_vec();
        for j in 1..m {
            for i in j..m {
                let ni_sq = (n_steps[i] as f64).powi(2);
                let nj_sq = (n_steps[i - j] as f64).powi(2);
                if (ni_sq - nj_sq).abs() < f64::EPSILON {
                    continue;
                }
                let ratio = ni_sq / nj_sq;
                let prev = t[i - 1].clone();
                for d in 0..dim {
                    t[i][d] = (ratio * t[i][d] - prev[d]) / (ratio - 1.0);
                }
            }
        }
        let best = t[m - 1].clone();
        let second = t[m - 2].clone();
        let err: f64 = best
            .iter()
            .zip(second.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        (best, err)
    }
}
/// PI step-size controller for adaptive ODE integration.
///
/// Uses a proportional-integral (PI) controller to select the next step size
/// based on the ratio of the current error estimate to the tolerance.
///
/// Reference: Hairer, Nørsett & Wanner, "Solving ODEs I", Section II.4.
pub struct PiStepController {
    /// Absolute tolerance.
    pub atol: f64,
    /// Relative tolerance.
    pub rtol: f64,
    /// Minimum allowed step size.
    pub dt_min: f64,
    /// Maximum allowed step size.
    pub dt_max: f64,
    /// Safety factor (< 1, typically 0.8–0.9).
    pub safety: f64,
    /// Previous error ratio (for PI control).
    pub(super) prev_err_ratio: f64,
    /// Proportional exponent (β₁ = 0.7/p for method of order p).
    pub(super) beta1: f64,
    /// Integral exponent (β₂ = -0.4/p).
    pub(super) beta2: f64,
}
impl PiStepController {
    /// Create a new PI controller.
    ///
    /// # Arguments
    /// * `atol`   - Absolute tolerance.
    /// * `rtol`   - Relative tolerance.
    /// * `dt_min` - Minimum step size.
    /// * `dt_max` - Maximum step size.
    /// * `safety` - Safety factor (e.g., 0.9).
    pub fn new(atol: f64, rtol: f64, dt_min: f64, dt_max: f64, safety: f64) -> Self {
        Self {
            atol,
            rtol,
            dt_min,
            dt_max,
            safety,
            prev_err_ratio: 1.0,
            beta1: 0.7 / 5.0,
            beta2: -0.4 / 5.0,
        }
    }
    /// Compute the next step size and decide acceptance.
    ///
    /// # Arguments
    /// * `dt`  - Current step size used.
    /// * `err` - Normalized error estimate (mixed absolute/relative norm).
    ///
    /// Returns `(new_dt, accepted)`.
    pub fn control(&self, dt: f64, err: f64) -> (f64, bool) {
        let err_ratio = if err < 1e-300 { 1e-300 } else { err };
        let factor =
            self.safety * err_ratio.powf(-self.beta1) * self.prev_err_ratio.powf(-self.beta2);
        let factor = factor.clamp(0.1, 5.0);
        let new_dt = (dt * factor).clamp(self.dt_min, self.dt_max);
        let accepted = err <= 1.0;
        (new_dt, accepted)
    }
}
/// Runge-Kutta-Fehlberg 4(5) adaptive integrator.
///
/// Uses the classic Fehlberg coefficients to produce an embedded 4th/5th-order
/// pair, adapting the step size so the local truncation error stays within
/// `atol + rtol * |y|`.
pub struct RungeKuttaFehlberg {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Minimum permitted step size.
    pub h_min: f64,
    /// Maximum permitted step size.
    pub h_max: f64,
}
impl RungeKuttaFehlberg {
    /// Create a new RKF45 integrator with the given tolerances and step bounds.
    pub fn new(rtol: f64, atol: f64, h_min: f64, h_max: f64) -> Self {
        Self {
            rtol,
            atol,
            h_min,
            h_max,
        }
    }
    /// Perform one RKF45 step from state `y` at time `t` with step size `h`.
    ///
    /// Returns `(y4, y5, h_new)` where `y4` is the 4th-order solution,
    /// `y5` the 5th-order solution, and `h_new` a suggested next step size.
    ///
    /// The RKF45 Butcher tableau uses the classic Fehlberg coefficients.
    pub fn step<F>(&self, f: &F, t: f64, y: &[f64], h: f64) -> (Vec<f64>, Vec<f64>, f64)
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = y.len();
        let k1 = f(t, y);
        let y2: Vec<f64> = (0..n).map(|i| y[i] + h * (1.0 / 4.0) * k1[i]).collect();
        let k2 = f(t + h / 4.0, &y2);
        let y3: Vec<f64> = (0..n)
            .map(|i| y[i] + h * (3.0 / 32.0 * k1[i] + 9.0 / 32.0 * k2[i]))
            .collect();
        let k3 = f(t + 3.0 * h / 8.0, &y3);
        let y4s: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (1932.0 / 2197.0 * k1[i] - 7200.0 / 2197.0 * k2[i] + 7296.0 / 2197.0 * k3[i])
            })
            .collect();
        let k4 = f(t + 12.0 * h / 13.0, &y4s);
        let y5s: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (439.0 / 216.0 * k1[i] - 8.0 * k2[i] + 3680.0 / 513.0 * k3[i]
                        - 845.0 / 4104.0 * k4[i])
            })
            .collect();
        let k5 = f(t + h, &y5s);
        let y6s: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (-8.0 / 27.0 * k1[i] + 2.0 * k2[i] - 3544.0 / 2565.0 * k3[i]
                        + 1859.0 / 4104.0 * k4[i]
                        - 11.0 / 40.0 * k5[i])
            })
            .collect();
        let k6 = f(t + h / 2.0, &y6s);
        let y4: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (25.0 / 216.0 * k1[i] + 1408.0 / 2565.0 * k3[i] + 2197.0 / 4104.0 * k4[i]
                        - 1.0 / 5.0 * k5[i])
            })
            .collect();
        let y5: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (16.0 / 135.0 * k1[i] + 6656.0 / 12825.0 * k3[i] + 28561.0 / 56430.0 * k4[i]
                        - 9.0 / 50.0 * k5[i]
                        + 2.0 / 55.0 * k6[i])
            })
            .collect();
        let err_rms = if n == 0 {
            0.0
        } else {
            let sum_sq: f64 = (0..n)
                .map(|i| {
                    let sc = self.atol + self.rtol * y[i].abs().max(y4[i].abs());
                    let e = (y5[i] - y4[i]) / sc;
                    e * e
                })
                .sum();
            (sum_sq / n as f64).sqrt()
        };
        let h_new = if err_rms == 0.0 {
            (h * 5.0).min(self.h_max)
        } else {
            let factor = 0.9 * err_rms.powf(-0.2);
            (h * factor.clamp(0.1, 5.0)).clamp(self.h_min, self.h_max)
        };
        (y4, y5, h_new)
    }
    /// Integrate `f(t, y)` from `t0` to `tf` starting from `y0`.
    ///
    /// Returns a vector of `(t, y)` pairs at every accepted step.
    pub fn integrate<F>(&self, f: F, t0: f64, tf: f64, y0: Vec<f64>) -> Vec<(f64, Vec<f64>)>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut t = t0;
        let mut y = y0;
        let mut h = ((tf - t0) / 100.0).clamp(self.h_min, self.h_max);
        let mut result = vec![(t, y.clone())];
        while t < tf {
            if t + h > tf {
                h = tf - t;
            }
            let (y4, y5, h_new) = self.step(&f, t, &y, h);
            let err: f64 = if y4.is_empty() {
                0.0
            } else {
                let s: f64 = y4
                    .iter()
                    .zip(y5.iter())
                    .map(|(a, b)| {
                        let sc = self.atol + self.rtol * a.abs().max(b.abs());
                        let e = (b - a) / sc;
                        e * e
                    })
                    .sum();
                (s / y4.len() as f64).sqrt()
            };
            if err <= 1.0 || h <= self.h_min {
                t += h;
                y = y4;
                result.push((t, y.clone()));
            }
            h = h_new.min((tf - t).max(0.0));
            if h < self.h_min && t < tf {
                h = self.h_min;
            }
        }
        result
    }
}
/// 4th-order Adams-Bashforth explicit multi-step integrator.
///
/// Requires 4 previous derivative evaluations (bootstrapped via RK4).
/// Call [`Adams4::push`] to record derivative values and [`Adams4::step`] to
/// advance once the history buffer is full.
pub struct Adams4 {
    /// History buffer: `[(t, f(t, y))]` oldest-first, capacity 4.
    pub(super) history: Vec<Vec<f64>>,
}
impl Adams4 {
    /// Create a new `Adams4` integrator with an empty history buffer.
    pub fn new() -> Self {
        Self {
            history: Vec::with_capacity(4),
        }
    }
    /// Push a derivative evaluation `dydt` into the history buffer (oldest entry
    /// is dropped when the buffer exceeds 4 entries).
    pub fn push(&mut self, dydt: Vec<f64>) {
        if self.history.len() >= 4 {
            self.history.remove(0);
        }
        self.history.push(dydt);
    }
    /// Returns `true` if the history buffer contains 4 evaluations (ready to step).
    pub fn ready(&self) -> bool {
        self.history.len() == 4
    }
    /// Perform one Adams-Bashforth 4th-order step.
    ///
    /// Panics if `ready()` returns `false`.
    ///
    /// # Formula
    /// `y_new = y + (dt/24) * (55*f3 - 59*f2 + 37*f1 - 9*f0)`
    pub fn step(&self, y: &[f64], dt: f64) -> Vec<f64> {
        assert!(
            self.ready(),
            "Adams4: need 4 history entries before stepping"
        );
        let n = y.len();
        let f0 = &self.history[0];
        let f1 = &self.history[1];
        let f2 = &self.history[2];
        let f3 = &self.history[3];
        (0..n)
            .map(|i| {
                y[i] + (dt / 24.0) * (55.0 * f3[i] - 59.0 * f2[i] + 37.0 * f1[i] - 9.0 * f0[i])
            })
            .collect()
    }
}
/// Adams-Bashforth-Moulton 4th-order predictor-corrector integrator.
///
/// Requires 4 history points (bootstrap with RK4). Uses Adams-Bashforth 4 as
/// predictor and Adams-Moulton 4 as corrector (PECE scheme).
pub struct AdamsBashforthMoulton4 {
    /// Circular history buffer of (t, y, f) tuples, oldest first.
    pub(super) history: std::collections::VecDeque<(f64, Vec<f64>, Vec<f64>)>,
}
impl AdamsBashforthMoulton4 {
    /// Create a new Adams-Bashforth-Moulton4 solver.
    pub fn new() -> Self {
        Self {
            history: std::collections::VecDeque::with_capacity(4),
        }
    }
    /// Push a history entry (t, y, f(t,y)).
    pub fn push_history(&mut self, t: f64, y: Vec<f64>, fy: Vec<f64>) {
        if self.history.len() == 4 {
            self.history.pop_front();
        }
        self.history.push_back((t, y, fy));
    }
    /// Returns true when enough history is available for a full step.
    pub fn ready(&self) -> bool {
        self.history.len() == 4
    }
    /// Perform one Adams-Bashforth-Moulton PECE step.
    ///
    /// # Arguments
    /// * `f`  - RHS function `f(t, y) -> Vec`f64`.
    /// * `t`  - Current time (= history.back().t + dt conceptually).
    /// * `y`  - Current state (last history entry's y should match).
    /// * `dt` - Step size.
    pub fn step<F>(&self, f: &F, t: f64, y: &[f64], dt: f64) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = y.len();
        let h = self.history.iter().collect::<Vec<_>>();
        if h.len() < 4 {
            let fy = f(t, y);
            return (0..n).map(|i| y[i] + dt * fy[i]).collect();
        }
        let f3 = &h[3].2;
        let f2 = &h[2].2;
        let f1 = &h[1].2;
        let f0 = &h[0].2;
        let y_pred: Vec<f64> = (0..n)
            .map(|i| y[i] + dt / 24.0 * (55.0 * f3[i] - 59.0 * f2[i] + 37.0 * f1[i] - 9.0 * f0[i]))
            .collect();
        let t_new = t + dt;
        let f_pred = f(t_new, &y_pred);
        (0..n)
            .map(|i| y[i] + dt / 24.0 * (9.0 * f_pred[i] + 19.0 * f3[i] - 5.0 * f2[i] + f1[i]))
            .collect()
    }
}
/// Adaptive ODE integrator built on top of [`DormandPrince45`].
///
/// Wraps any [`OdeSystem`] and integrates it from `t0` to `t_end` with
/// automatic step-size control via error estimation and the PI controller
/// safety factors from DormandPrince45.
pub struct AdaptiveIntegrator {
    /// Absolute tolerance.
    pub atol: f64,
    /// Relative tolerance.
    pub rtol: f64,
    /// Minimum step size (integration stops if `dt < dt_min`).
    pub dt_min: f64,
    /// Maximum step size.
    pub dt_max: f64,
}
impl AdaptiveIntegrator {
    /// Create a new `AdaptiveIntegrator` with default tolerances.
    pub fn new(atol: f64, rtol: f64) -> Self {
        Self {
            atol,
            rtol,
            dt_min: 1e-15,
            dt_max: f64::INFINITY,
        }
    }
    /// Integrate `system` from time `t0` to `t_end` starting at state `y0`
    /// with initial step `dt_init`.
    ///
    /// Returns a `Vec<(t, y)>` of accepted steps including the start point.
    pub fn integrate(
        &self,
        system: &impl OdeSystem,
        y0: Vec<f64>,
        t0: f64,
        t_end: f64,
        dt_init: f64,
    ) -> Vec<(f64, Vec<f64>)> {
        let stepper = DormandPrince45::new(self.atol, self.rtol);
        let n = y0.len();
        let mut t = t0;
        let mut y = y0;
        let mut dt = dt_init.min(self.dt_max);
        let mut result = vec![(t, y.clone())];
        let safety = 0.9_f64;
        let min_scale = 0.2_f64;
        let max_scale = 10.0_f64;
        while t < t_end {
            if t + dt > t_end {
                dt = t_end - t;
            }
            let f_closure = |yv: &[f64]| {
                let mut dydt = vec![0.0_f64; n];
                system.ode_rhs(t, yv, &mut dydt);
                dydt
            };
            let res = stepper.step(&y, f_closure, dt);
            let tol_norm = (self.atol
                + self.rtol * y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max))
            .max(self.atol);
            let err_norm = res.error_estimate / tol_norm;
            if res.accepted {
                t += dt;
                y = res.y;
                result.push((t, y.clone()));
                let scale = if res.error_estimate == 0.0 {
                    max_scale
                } else {
                    (safety * err_norm.powf(-0.2)).clamp(min_scale, max_scale)
                };
                dt = (dt * scale).min(self.dt_max);
            } else {
                let scale = (safety * err_norm.powf(-0.25)).clamp(min_scale, 1.0);
                dt *= scale;
            }
            if dt < self.dt_min {
                break;
            }
        }
        result
    }
}
/// Event detection utilities for ODE integration.
///
/// Detects zero crossings of a scalar function using coarse scanning
/// followed by bisection refinement.
pub struct EventDetector;
impl EventDetector {
    /// Find all zero crossings of `g(t)` on `\[t_start, t_end\]`.
    ///
    /// First scans with step `scan_dt`, then refines each sign change
    /// with bisection to tolerance `bisect_tol`.
    pub fn find_zero_crossings<G>(
        g: G,
        t_start: f64,
        t_end: f64,
        scan_dt: f64,
        bisect_tol: f64,
    ) -> Vec<f64>
    where
        G: Fn(f64) -> f64,
    {
        let mut crossings = Vec::new();
        let mut t = t_start;
        let mut g_prev = g(t);
        while t + scan_dt <= t_end {
            let t_next = (t + scan_dt).min(t_end);
            let g_next = g(t_next);
            if g_prev * g_next < 0.0 {
                let root = Self::bisect(&g, t, t_next, bisect_tol);
                crossings.push(root);
            }
            g_prev = g_next;
            t = t_next;
        }
        crossings
    }
    /// Bisection root-finding for `g(t) = 0` on `\[a, b\]` (assumes sign change).
    pub fn bisect<G>(g: &G, mut a: f64, mut b: f64, tol: f64) -> f64
    where
        G: Fn(f64) -> f64,
    {
        let max_iter = 100usize;
        for _ in 0..max_iter {
            let mid = 0.5 * (a + b);
            if (b - a).abs() < tol {
                return mid;
            }
            let gm = g(mid);
            let ga = g(a);
            if ga * gm <= 0.0 {
                b = mid;
            } else {
                a = mid;
            }
        }
        0.5 * (a + b)
    }
}
/// Kick-drift LeapFrog (symplectic) integrator.
///
/// Separate kick (velocity half-update) and drift (position full-update) steps
/// allow operator-splitting composition for symplectic integration.
pub struct LeapFrog;
impl LeapFrog {
    /// Create a new `LeapFrog` integrator.
    pub fn new() -> Self {
        Self
    }
    /// Velocity kick: `v += a * dt`
    pub fn kick(v: &mut Vec<f64>, a: &[f64], dt: f64) {
        for (vi, ai) in v.iter_mut().zip(a.iter()) {
            *vi += ai * dt;
        }
    }
    /// Position drift: `x += v * dt`
    pub fn drift(x: &mut Vec<f64>, v: &[f64], dt: f64) {
        for (xi, vi) in x.iter_mut().zip(v.iter()) {
            *xi += vi * dt;
        }
    }
}
/// Implicit trapezoidal method (Crank–Nicolson) for stiff ODE systems.
///
/// Solves `y_{n+1} = y_n + (dt/2) * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))`
/// using Newton iteration with a finite-difference Jacobian.
pub struct StiffOdeSolver {
    /// Implicit weight γ in `\[0.5, 1\]`; 0.5 = trapezoidal, 1 = implicit Euler.
    pub gamma: f64,
}
impl StiffOdeSolver {
    /// Create a new `StiffOdeSolver`.
    pub fn new(gamma: f64) -> Self {
        Self { gamma }
    }
    /// Perform one implicit step of size `dt` from state `y` at time `t`.
    ///
    /// Uses Newton iteration with a finite-difference Jacobian approximation.
    pub fn step_fixed<F>(&self, f: F, t: f64, y: &[f64], dt: f64) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = y.len();
        let t_new = t + dt;
        let f0 = f(t, y);
        let gamma = self.gamma;
        let mut z: Vec<f64> = (0..n).map(|i| y[i] + dt * f0[i]).collect();
        let fd_eps = 1e-7;
        let max_iter = 50usize;
        let tol = 1e-10_f64;
        for _ in 0..max_iter {
            let fz = f(t_new, &z);
            let g: Vec<f64> = (0..n)
                .map(|i| z[i] - y[i] - dt * ((1.0 - gamma) * f0[i] + gamma * fz[i]))
                .collect();
            let g_norm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if g_norm < tol {
                break;
            }
            let mut jac = vec![vec![0.0_f64; n]; n];
            for j in 0..n {
                let mut z_pert = z.clone();
                z_pert[j] += fd_eps;
                let f_pert = f(t_new, &z_pert);
                for i in 0..n {
                    let dfdz = (f_pert[i] - fz[i]) / fd_eps;
                    jac[i][j] = if i == j { 1.0 } else { 0.0 } - dt * gamma * dfdz;
                }
            }
            let neg_g: Vec<f64> = g.iter().map(|v| -v).collect();
            if let Some(dz) = solve_linear_system(&jac, &neg_g) {
                for i in 0..n {
                    z[i] += dz[i];
                }
            } else {
                break;
            }
        }
        z
    }
}
/// Backward Differentiation Formula of order 2 (BDF2) implicit integrator.
///
/// Suitable for stiff ODEs. Uses Newton iteration with a finite-difference
/// Jacobian approximation. Requires one previous step to be stored.
pub struct Bdf2 {
    /// Previous state (needed for the BDF2 formula).
    pub y_prev: Option<Vec<f64>>,
    /// Maximum Newton iterations.
    pub max_iter: usize,
    /// Newton convergence tolerance.
    pub newton_tol: f64,
}
impl Bdf2 {
    /// Create a new `Bdf2` integrator.
    pub fn new() -> Self {
        Self {
            y_prev: None,
            max_iter: 50,
            newton_tol: 1e-10,
        }
    }
    /// Prime the integrator with the initial state `y0` (stores `y_prev`).
    pub fn prime(&mut self, y0: Vec<f64>) {
        self.y_prev = Some(y0);
    }
    /// Perform one BDF2 step: solve
    /// `(3/2)*y_new - 2*y_cur + (1/2)*y_prev = dt * f(t_new, y_new)`
    ///
    /// Uses Newton iteration with a finite-difference Jacobian.
    /// Returns the new state, or falls back to BDF1 (implicit Euler) if
    /// `y_prev` is not set.
    pub fn step(
        &mut self,
        y: &[f64],
        f: impl Fn(f64, &[f64]) -> Vec<f64>,
        t_new: f64,
        dt: f64,
    ) -> Vec<f64> {
        let n = y.len();
        let h = dt;
        let (alpha, beta, rhs_const) = match &self.y_prev {
            Some(yp) => {
                let rc: Vec<f64> = (0..n).map(|i| 2.0 * y[i] - 0.5 * yp[i]).collect();
                (1.5_f64, h, rc)
            }
            None => {
                let rc = y.to_vec();
                (1.0_f64, h, rc)
            }
        };
        let mut z = y.to_vec();
        let fd_eps = 1e-7;
        for _ in 0..self.max_iter {
            let fz = f(t_new, &z);
            let g: Vec<f64> = (0..n)
                .map(|i| alpha * z[i] - rhs_const[i] - beta * fz[i])
                .collect();
            let g_norm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if g_norm < self.newton_tol {
                break;
            }
            let mut jac = vec![vec![0.0_f64; n]; n];
            for j in 0..n {
                let mut z_pert = z.clone();
                z_pert[j] += fd_eps;
                let f_pert = f(t_new, &z_pert);
                for i in 0..n {
                    let dfdz_ij = (f_pert[i] - fz[i]) / fd_eps;
                    jac[i][j] = if i == j { alpha } else { 0.0 } - beta * dfdz_ij;
                }
            }
            let neg_g: Vec<f64> = g.iter().map(|v| -v).collect();
            if let Some(dz) = solve_linear_system(&jac, &neg_g) {
                for i in 0..n {
                    z[i] += dz[i];
                }
            } else {
                break;
            }
        }
        self.y_prev = Some(y.to_vec());
        z
    }
}
/// Implicit Euler integrator using Newton iteration for the nonlinear solve.
///
/// Solves `y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1})` via Newton's method
/// with a finite-difference Jacobian.
pub struct ImplicitEulerNewton {
    /// Maximum Newton iterations.
    pub max_iter: usize,
    /// Convergence tolerance on the residual norm.
    pub tol: f64,
}
impl ImplicitEulerNewton {
    /// Create a new `ImplicitEulerNewton` solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }
    /// Perform one implicit Euler step of size `dt`.
    pub fn step<F>(&self, f: F, t: f64, y: &[f64], dt: f64) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = y.len();
        let t_new = t + dt;
        let fd_eps = 1e-7_f64;
        let f0 = f(t, y);
        let mut z: Vec<f64> = (0..n).map(|i| y[i] + dt * f0[i]).collect();
        for _ in 0..self.max_iter {
            let fz = f(t_new, &z);
            let g: Vec<f64> = (0..n).map(|i| z[i] - y[i] - dt * fz[i]).collect();
            let g_norm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if g_norm < self.tol {
                break;
            }
            let mut jac = vec![vec![0.0_f64; n]; n];
            for j in 0..n {
                let mut z_pert = z.clone();
                z_pert[j] += fd_eps;
                let f_pert = f(t_new, &z_pert);
                for i in 0..n {
                    let dfdz = (f_pert[i] - fz[i]) / fd_eps;
                    jac[i][j] = if i == j { 1.0 } else { 0.0 } - dt * dfdz;
                }
            }
            let neg_g: Vec<f64> = g.iter().map(|v| -v).collect();
            if let Some(dz) = solve_linear_system(&jac, &neg_g) {
                for i in 0..n {
                    z[i] += dz[i];
                }
            } else {
                break;
            }
        }
        z
    }
}
/// Backward Differentiation Formula of order 2 (BDF2).
///
/// BDF2 step: `(3/2)*y_{n+1} - 2*y_n + (1/2)*y_{n-1} = dt * f(t_{n+1}, y_{n+1})`
///
/// Uses Newton iteration with finite-difference Jacobian.
pub struct BdfOrder2 {
    /// Maximum Newton iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}
impl BdfOrder2 {
    /// Create a new BDF2 solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }
    /// Perform the first step using implicit Euler (bootstrapping BDF2).
    pub fn first_step<F>(&self, f: F, t: f64, y: &[f64], dt: f64) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        ImplicitEulerNewton::new(self.max_iter, self.tol).step(f, t, y, dt)
    }
    /// Perform one BDF2 step given the two previous states.
    ///
    /// # Arguments
    /// * `f`     - RHS function `f(t, y)`.
    /// * `t`     - Current time `t_n` (step will be to `t_{n+1} = t_n + dt`).
    /// * `y_nm1` - State at `t_{n-1}`.
    /// * `y_n`   - State at `t_n`.
    /// * `dt`    - Step size.
    pub fn step<F>(&self, f: F, t: f64, y_nm1: &[f64], y_n: &[f64], dt: f64) -> Vec<f64>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let n = y_n.len();
        let t_new = t + dt;
        let fd_eps = 1e-7_f64;
        let mut z: Vec<f64> = (0..n).map(|i| 2.0 * y_n[i] - y_nm1[i]).collect();
        for _ in 0..self.max_iter {
            let fz = f(t_new, &z);
            let g: Vec<f64> = (0..n)
                .map(|i| 1.5 * z[i] - 2.0 * y_n[i] + 0.5 * y_nm1[i] - dt * fz[i])
                .collect();
            let g_norm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if g_norm < self.tol {
                break;
            }
            let mut jac = vec![vec![0.0_f64; n]; n];
            for j in 0..n {
                let mut z_pert = z.clone();
                z_pert[j] += fd_eps;
                let f_pert = f(t_new, &z_pert);
                for i in 0..n {
                    let dfdz = (f_pert[i] - fz[i]) / fd_eps;
                    jac[i][j] = if i == j { 1.5 } else { 0.0 } - dt * dfdz;
                }
            }
            let neg_g: Vec<f64> = g.iter().map(|v| -v).collect();
            if let Some(dz) = solve_linear_system(&jac, &neg_g) {
                for i in 0..n {
                    z[i] += dz[i];
                }
            } else {
                break;
            }
        }
        z
    }
}
/// Runge-Kutta-Fehlberg 4(5) integrator.
///
/// Uses the classic RKF45 Butcher tableau. Returns `(y_5th, error_estimate)`
/// where the error is the RMS difference between the 4th- and 5th-order
/// solutions.
pub struct Fehlberg45;
impl Fehlberg45 {
    /// Create a new `Fehlberg45` integrator.
    pub fn new() -> Self {
        Self
    }
    /// Perform a single RKF45 step from state `y` with step size `dt`.
    pub fn step(&self, y: &[f64], f: impl Fn(&[f64]) -> Vec<f64>, dt: f64) -> (Vec<f64>, f64) {
        let n = y.len();
        let a21 = 1.0 / 4.0;
        let a31 = 3.0 / 32.0;
        let a32 = 9.0 / 32.0;
        let a41 = 1932.0 / 2197.0;
        let a42 = -7200.0 / 2197.0;
        let a43 = 7296.0 / 2197.0;
        let a51 = 439.0 / 216.0;
        let a52 = -8.0;
        let a53 = 3680.0 / 513.0;
        let a54 = -845.0 / 4104.0;
        let a61 = -8.0 / 27.0;
        let a62 = 2.0;
        let a63 = -3544.0 / 2565.0;
        let a64 = 1859.0 / 4104.0;
        let a65 = -11.0 / 40.0;
        let b1 = 16.0 / 135.0;
        let b3 = 6656.0 / 12825.0;
        let b4 = 28561.0 / 56430.0;
        let b5 = -9.0 / 50.0;
        let b6 = 2.0 / 55.0;
        let c1 = 25.0 / 216.0;
        let c3 = 1408.0 / 2565.0;
        let c4 = 2197.0 / 4104.0;
        let c5 = -1.0 / 5.0;
        let k1 = f(y);
        let y2: Vec<f64> = (0..n).map(|i| y[i] + dt * a21 * k1[i]).collect();
        let k2 = f(&y2);
        let y3: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (a31 * k1[i] + a32 * k2[i]))
            .collect();
        let k3 = f(&y3);
        let y4: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]))
            .collect();
        let k4 = f(&y4);
        let y5: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]))
            .collect();
        let k5 = f(&y5);
        let y6: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + dt * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i])
            })
            .collect();
        let k6 = f(&y6);
        let y5th: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]))
            .collect();
        let y4th: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i]))
            .collect();
        let err = {
            let sum_sq: f64 = y5th
                .iter()
                .zip(y4th.iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum();
            (sum_sq / n as f64).sqrt()
        };
        (y5th, err)
    }
}
/// Result of a single adaptive ODE step.
#[derive(Debug, Clone)]
pub struct StepResult {
    /// The new state vector after the step.
    pub y: Vec<f64>,
    /// RMS error estimate (difference between high- and low-order solutions).
    pub error_estimate: f64,
    /// Whether the step was accepted (error within tolerance).
    pub accepted: bool,
}
/// Dense output segment for continuous interpolation between ODE steps.
///
/// Uses a 3rd-order Hermite interpolant based on the endpoint states and
/// derivatives (k1, k6 from Dormand-Prince).
pub struct DenseOutputSegment {
    /// Start time of this segment.
    pub t0: f64,
    /// End time of this segment.
    pub t1: f64,
    /// State at `t0`.
    pub y0: Vec<f64>,
    /// State at `t1`.
    pub y1: Vec<f64>,
    /// Derivative at `t0` (first stage slope k1).
    pub k1: Vec<f64>,
    /// Derivative at `t1` (last stage slope k6).
    pub k6: Vec<f64>,
}
impl DenseOutputSegment {
    /// Create a new dense output segment.
    pub fn new(t0: f64, t1: f64, y0: Vec<f64>, y1: Vec<f64>, k1: Vec<f64>, k6: Vec<f64>) -> Self {
        Self {
            t0,
            t1,
            y0,
            y1,
            k1,
            k6,
        }
    }
    /// Evaluate the interpolant at time `t` in `\[t0, t1\]`.
    ///
    /// Uses cubic Hermite interpolation:
    /// `y(t0 + θ*h) = (1-θ)*y0 + θ*y1 + θ*(θ-1)*\[(1-2θ)*(y1-y0) + (θ-1)*h*k1 + θ*h*k6\]`
    pub fn evaluate(&self, t: f64) -> Vec<f64> {
        let h = self.t1 - self.t0;
        let theta = if h.abs() < 1e-15 {
            0.0
        } else {
            (t - self.t0) / h
        };
        let n = self.y0.len();
        (0..n)
            .map(|i| {
                let dy = self.y1[i] - self.y0[i];
                (1.0 - theta) * self.y0[i]
                    + theta * self.y1[i]
                    + theta
                        * (theta - 1.0)
                        * ((1.0 - 2.0 * theta) * dy
                            + (theta - 1.0) * h * self.k1[i]
                            + theta * h * self.k6[i])
            })
            .collect()
    }
}
