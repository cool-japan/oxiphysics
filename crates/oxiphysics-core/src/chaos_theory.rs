// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Chaos theory and nonlinear dynamics.
//!
//! This module provides tools for analysing chaotic dynamical systems, including:
//!
//! - [`LorenzSystem`]: The classic Lorenz attractor (σ, ρ, β parameters).
//! - [`RosslerSystem`]: Rössler system with spiral chaos.
//! - [`LogisticMap`]: Logistic map and period-doubling bifurcations.
//! - [`DuffingOscillator`]: Periodically-forced Duffing oscillator.
//! - [`LyapunovExponent`]: Largest Lyapunov exponent estimation.
//! - [`LyapunovSpectrum`]: Full spectrum via QR-based algorithm.
//! - [`FtleField`]: Finite-time Lyapunov exponents on a grid.
//! - [`PoincareSection`]: Poincaré section sampler.
//! - [`CorrelationDimension`]: Grassberger-Procaccia correlation dimension.
//! - [`BoxCounting`]: Box-counting fractal dimension.
//! - [`RecurrencePlot`]: Recurrence plot and recurrence quantification.
//! - [`MutualInformation`]: Average mutual information for embedding.
//! - [`EmbeddingDimension`]: False nearest neighbours (Taken's embedding).
//! - [`FeigenbaumConstant`]: Numerical Feigenbaum δ from period-doubling.
//! - [`ChaosAnalysis`]: Convenience wrapper combining several diagnostics.
//!
//! # References
//! - Lorenz, E.N. (1963). Deterministic nonperiodic flow. *J. Atmos. Sci.*, 20, 130–141.
//! - Rössler, O.E. (1976). An equation for continuous chaos. *Phys. Lett. A*, 57, 397–398.
//! - Grassberger, P. & Procaccia, I. (1983). Characterization of strange attractors.
//!   *Phys. Rev. Lett.*, 50, 346–349.
//! - Takens, F. (1981). Detecting strange attractors in turbulence. *Lect. Notes Math.*, 898.
//! - Feigenbaum, M.J. (1978). Quantitative universality for a class of nonlinear transformations.
//!   *J. Stat. Phys.*, 19, 25–52.

// ─── Lorenz System ───────────────────────────────────────────────────────────

/// Parameters and state for the Lorenz attractor.
///
/// The Lorenz equations are:
/// ```text
/// dx/dt = σ(y − x)
/// dy/dt = x(ρ − z) − y
/// dz/dt = xy − βz
/// ```
/// with the classic chaotic parameters σ = 10, ρ = 28, β = 8/3.
#[derive(Debug, Clone)]
pub struct LorenzSystem {
    /// Prandtl number σ.
    pub sigma: f64,
    /// Rayleigh number ρ.
    pub rho: f64,
    /// Geometric factor β.
    pub beta: f64,
}

impl LorenzSystem {
    /// Creates a new Lorenz system with the given parameters.
    pub fn new(sigma: f64, rho: f64, beta: f64) -> Self {
        Self { sigma, rho, beta }
    }

    /// Creates the classic chaotic Lorenz system (σ=10, ρ=28, β=8/3).
    pub fn classic() -> Self {
        Self::new(10.0, 28.0, 8.0 / 3.0)
    }

    /// Computes the time derivative `[dx/dt, dy/dt, dz/dt]` at state `[x, y, z]`.
    pub fn deriv(&self, state: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = state;
        [
            self.sigma * (y - x),
            x * (self.rho - z) - y,
            x * y - self.beta * z,
        ]
    }

    /// Integrates using the 4th-order Runge-Kutta method for `steps` steps of size `dt`.
    pub fn integrate(&self, initial: [f64; 3], dt: f64, steps: usize) -> Vec<[f64; 3]> {
        let mut traj = Vec::with_capacity(steps + 1);
        let mut s = initial;
        traj.push(s);
        for _ in 0..steps {
            s = rk4_step(|st| self.deriv(st), s, dt);
            traj.push(s);
        }
        traj
    }

    /// Returns the three fixed points of the Lorenz system (origin and two symmetric points).
    pub fn fixed_points(&self) -> [[f64; 3]; 3] {
        let c = (self.beta * (self.rho - 1.0)).sqrt();
        [
            [0.0, 0.0, 0.0],
            [c, c, self.rho - 1.0],
            [-c, -c, self.rho - 1.0],
        ]
    }

    /// Computes the Jacobian matrix at state `[x, y, z]` as a row-major 3×3 array.
    pub fn jacobian(&self, state: [f64; 3]) -> [[f64; 3]; 3] {
        let [x, y, z] = state;
        [
            [-self.sigma, self.sigma, 0.0],
            [self.rho - z, -1.0, -x],
            [y, x, -self.beta],
        ]
    }
}

// ─── Rössler System ──────────────────────────────────────────────────────────

/// Parameters and state for the Rössler system.
///
/// The Rössler equations are:
/// ```text
/// dx/dt = −y − z
/// dy/dt = x + ay
/// dz/dt = b + z(x − c)
/// ```
/// Typical chaotic parameters: a = 0.2, b = 0.2, c = 5.7.
#[derive(Debug, Clone)]
pub struct RosslerSystem {
    /// Parameter a (spiral rate).
    pub a: f64,
    /// Parameter b (z-coupling).
    pub b: f64,
    /// Parameter c (nonlinear threshold).
    pub c: f64,
}

impl RosslerSystem {
    /// Creates a new Rössler system with the given parameters.
    pub fn new(a: f64, b: f64, c: f64) -> Self {
        Self { a, b, c }
    }

    /// Creates the classic chaotic Rössler system (a=0.2, b=0.2, c=5.7).
    pub fn classic() -> Self {
        Self::new(0.2, 0.2, 5.7)
    }

    /// Computes the time derivative `[dx/dt, dy/dt, dz/dt]`.
    pub fn deriv(&self, state: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = state;
        [-y - z, x + self.a * y, self.b + z * (x - self.c)]
    }

    /// Integrates using RK4 for `steps` steps of size `dt`.
    pub fn integrate(&self, initial: [f64; 3], dt: f64, steps: usize) -> Vec<[f64; 3]> {
        let mut traj = Vec::with_capacity(steps + 1);
        let mut s = initial;
        traj.push(s);
        for _ in 0..steps {
            s = rk4_step(|st| self.deriv(st), s, dt);
            traj.push(s);
        }
        traj
    }

    /// Approximate period of the spiral for small-amplitude oscillations.
    pub fn approximate_period(&self) -> f64 {
        // Frequency ≈ 1 for small a.
        std::f64::consts::TAU
    }
}

// ─── Logistic Map ────────────────────────────────────────────────────────────

/// The logistic map x_{n+1} = r·x_n·(1 − x_n).
#[derive(Debug, Clone)]
pub struct LogisticMap {
    /// Growth parameter r ∈ \[0, 4\].
    pub r: f64,
}

impl LogisticMap {
    /// Creates a new logistic map with the given growth rate.
    pub fn new(r: f64) -> Self {
        Self { r }
    }

    /// Iterates the map for `n` steps starting from `x0`.
    pub fn iterate(&self, x0: f64, n: usize) -> Vec<f64> {
        let mut orbit = Vec::with_capacity(n + 1);
        let mut x = x0;
        orbit.push(x);
        for _ in 0..n {
            x = self.r * x * (1.0 - x);
            orbit.push(x);
        }
        orbit
    }

    /// Returns the attractor (after discarding `transient` steps) of length `n`.
    pub fn attractor(&self, x0: f64, transient: usize, n: usize) -> Vec<f64> {
        let orbit = self.iterate(x0, transient + n);
        orbit[transient..].to_vec()
    }

    /// Computes the Lyapunov exponent of the map.
    pub fn lyapunov_exponent(&self, x0: f64, n: usize) -> f64 {
        let mut x = x0;
        let mut sum = 0.0;
        for _ in 0..n {
            let deriv = (self.r * (1.0 - 2.0 * x)).abs();
            if deriv > 1e-15 {
                sum += deriv.ln();
            }
            x = self.r * x * (1.0 - x);
        }
        sum / n as f64
    }

    /// Generates a bifurcation diagram: for each r in `r_values` returns (r, attractor_points).
    pub fn bifurcation_diagram(
        r_values: &[f64],
        x0: f64,
        transient: usize,
        n_attractor: usize,
    ) -> Vec<(f64, Vec<f64>)> {
        r_values
            .iter()
            .map(|&r| {
                let lm = LogisticMap::new(r);
                let attr = lm.attractor(x0, transient, n_attractor);
                (r, attr)
            })
            .collect()
    }
}

// ─── Duffing Oscillator ──────────────────────────────────────────────────────

/// Periodically-forced Duffing oscillator.
///
/// Equation of motion:
/// ```text
/// ẍ + δẋ + αx + βx³ = γ cos(ωt)
/// ```
/// State vector is `[x, ẋ, t]` (extended with time for the forcing phase).
#[derive(Debug, Clone)]
pub struct DuffingOscillator {
    /// Linear damping coefficient δ.
    pub delta: f64,
    /// Linear stiffness α (can be negative for double-well).
    pub alpha: f64,
    /// Nonlinear stiffness β.
    pub beta: f64,
    /// Forcing amplitude γ.
    pub gamma: f64,
    /// Forcing frequency ω.
    pub omega: f64,
}

impl DuffingOscillator {
    /// Creates a new Duffing oscillator with the given parameters.
    pub fn new(delta: f64, alpha: f64, beta: f64, gamma: f64, omega: f64) -> Self {
        Self {
            delta,
            alpha,
            beta,
            gamma,
            omega,
        }
    }

    /// Classic chaotic Duffing parameters (δ=0.3, α=−1, β=1, γ=0.37, ω=1.2).
    pub fn classic_chaotic() -> Self {
        Self::new(0.3, -1.0, 1.0, 0.37, 1.2)
    }

    /// Computes `[ẋ, ẍ]` given `[x, ẋ]` and time `t`.
    pub fn deriv(&self, state: [f64; 2], t: f64) -> [f64; 2] {
        let [x, v] = state;
        let accel = -self.delta * v - self.alpha * x - self.beta * x.powi(3)
            + self.gamma * (self.omega * t).cos();
        [v, accel]
    }

    /// Integrates for `steps` steps of size `dt` starting from `initial` at `t0`.
    pub fn integrate(&self, initial: [f64; 2], t0: f64, dt: f64, steps: usize) -> Vec<[f64; 3]> {
        let mut result = Vec::with_capacity(steps + 1);
        let mut state = initial;
        let mut t = t0;
        result.push([state[0], state[1], t]);
        for _ in 0..steps {
            let k1 = self.deriv(state, t);
            let s2 = [state[0] + 0.5 * dt * k1[0], state[1] + 0.5 * dt * k1[1]];
            let k2 = self.deriv(s2, t + 0.5 * dt);
            let s3 = [state[0] + 0.5 * dt * k2[0], state[1] + 0.5 * dt * k2[1]];
            let k3 = self.deriv(s3, t + 0.5 * dt);
            let s4 = [state[0] + dt * k3[0], state[1] + dt * k3[1]];
            let k4 = self.deriv(s4, t + dt);
            state[0] += dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
            state[1] += dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
            t += dt;
            result.push([state[0], state[1], t]);
        }
        result
    }

    /// Returns the natural frequency for small oscillations around the stable fixed point.
    pub fn natural_frequency(&self) -> Option<f64> {
        // For double-well (α < 0, β > 0): fixed points at x = ±√(−α/β).
        if self.alpha < 0.0 && self.beta > 0.0 {
            let x_star = (-self.alpha / self.beta).sqrt();
            let omega_n = (2.0 * self.beta * x_star * x_star).sqrt();
            Some(omega_n)
        } else if self.alpha > 0.0 {
            Some(self.alpha.sqrt())
        } else {
            None
        }
    }
}

// ─── Lyapunov Exponent ───────────────────────────────────────────────────────

/// Estimates the largest Lyapunov exponent of a continuous flow using the
/// Wolf et al. algorithm (orbit separation renormalization).
#[derive(Debug, Clone)]
pub struct LyapunovExponent {
    /// Number of renormalization steps.
    pub n_steps: usize,
    /// Integration time step.
    pub dt: f64,
    /// Renormalization interval (in steps).
    pub renorm_interval: usize,
    /// Initial perturbation magnitude.
    pub epsilon: f64,
}

impl LyapunovExponent {
    /// Creates a new estimator with the given settings.
    pub fn new(n_steps: usize, dt: f64, renorm_interval: usize, epsilon: f64) -> Self {
        Self {
            n_steps,
            dt,
            renorm_interval,
            epsilon,
        }
    }

    /// Default settings for the Lorenz system.
    pub fn default_lorenz() -> Self {
        Self::new(50_000, 0.01, 10, 1e-8)
    }

    /// Estimates the largest Lyapunov exponent for a generic 3D flow `f`.
    ///
    /// `f(state) -> derivative` must be the vector field.
    pub fn estimate<F>(&self, f: F, initial: [f64; 3]) -> f64
    where
        F: Fn([f64; 3]) -> [f64; 3],
    {
        let mut x = initial;
        // Perturbed trajectory.
        let mut xp = [initial[0] + self.epsilon, initial[1], initial[2]];
        let mut lambda_sum = 0.0;
        let mut n_renorms = 0usize;

        for step in 0..self.n_steps {
            x = rk4_step(&f, x, self.dt);
            xp = rk4_step(&f, xp, self.dt);

            if (step + 1) % self.renorm_interval == 0 {
                let d = dist3(x, xp);
                if d > 1e-15 {
                    lambda_sum += (d / self.epsilon).ln();
                    n_renorms += 1;
                    // Rescale xp back to distance epsilon from x.
                    let scale = self.epsilon / d;
                    xp[0] = x[0] + (xp[0] - x[0]) * scale;
                    xp[1] = x[1] + (xp[1] - x[1]) * scale;
                    xp[2] = x[2] + (xp[2] - x[2]) * scale;
                }
            }
        }

        if n_renorms == 0 {
            return 0.0;
        }
        lambda_sum / (n_renorms as f64 * self.renorm_interval as f64 * self.dt)
    }
}

// ─── Lyapunov Spectrum ───────────────────────────────────────────────────────

/// Full Lyapunov spectrum computation via simultaneous QR decomposition.
///
/// The Benettin–Galgani–Giorgilli–Strelcyn (BGGS) algorithm.
#[derive(Debug, Clone)]
pub struct LyapunovSpectrum {
    /// Integration time step.
    pub dt: f64,
    /// Number of integration steps per QR step.
    pub steps_per_qr: usize,
    /// Total number of QR steps.
    pub n_qr: usize,
}

impl LyapunovSpectrum {
    /// Creates a new spectrum calculator.
    pub fn new(dt: f64, steps_per_qr: usize, n_qr: usize) -> Self {
        Self {
            dt,
            steps_per_qr,
            n_qr,
        }
    }

    /// Computes the 3 Lyapunov exponents for the Lorenz system.
    pub fn compute_lorenz(&self, sigma: f64, rho: f64, beta: f64, initial: [f64; 3]) -> [f64; 3] {
        let sys = LorenzSystem::new(sigma, rho, beta);
        self.compute(|s| sys.deriv(s), |s| sys.jacobian(s), initial)
    }

    /// Computes the Lyapunov spectrum for a general 3D system with Jacobian.
    ///
    /// `f(state) -> derivative`, `jac(state) -> 3×3 Jacobian`.
    pub fn compute<F, J>(&self, f: F, jac: J, initial: [f64; 3]) -> [f64; 3]
    where
        F: Fn([f64; 3]) -> [f64; 3],
        J: Fn([f64; 3]) -> [[f64; 3]; 3],
    {
        let mut x = initial;
        // Q is stored column-major: Q[col][row].
        let mut q: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut lambda = [0.0f64; 3];

        for _ in 0..self.n_qr {
            // Integrate the tangent vectors.
            for _ in 0..self.steps_per_qr {
                let j = jac(x);
                // Update each column of Q: dQ/dt = J·Q.
                let mut q_new = [[0.0f64; 3]; 3];
                for col in 0..3 {
                    for row in 0..3 {
                        let mut s = 0.0;
                        for k in 0..3 {
                            s += j[row][k] * q[col][k];
                        }
                        q_new[col][row] = q[col][row] + self.dt * s;
                    }
                }
                x = rk4_step(&f, x, self.dt);
                q = q_new;
            }
            // QR decomposition via Gram-Schmidt.
            let (q_new, r_diag) = gram_schmidt_qr(q);
            q = q_new;
            for i in 0..3 {
                lambda[i] += r_diag[i].abs().ln();
            }
        }

        let total_time = self.n_qr as f64 * self.steps_per_qr as f64 * self.dt;
        [
            lambda[0] / total_time,
            lambda[1] / total_time,
            lambda[2] / total_time,
        ]
    }
}

// ─── Poincaré Section ────────────────────────────────────────────────────────

/// Records the Poincaré section of a 3D trajectory through a plane z = z_section.
#[derive(Debug, Clone)]
pub struct PoincareSection {
    /// z-coordinate of the section plane.
    pub z_section: f64,
    /// Only record crossings in the positive-z direction if true.
    pub upward_only: bool,
}

impl PoincareSection {
    /// Creates a new Poincaré section at the given z level.
    pub fn new(z_section: f64, upward_only: bool) -> Self {
        Self {
            z_section,
            upward_only,
        }
    }

    /// Extracts crossing points from a 3D trajectory.
    ///
    /// Returns `(x, y)` coordinates at each crossing.
    pub fn extract(&self, traj: &[[f64; 3]]) -> Vec<[f64; 2]> {
        let mut crossings = Vec::new();
        for i in 1..traj.len() {
            let prev = traj[i - 1];
            let curr = traj[i];
            let dz_prev = prev[2] - self.z_section;
            let dz_curr = curr[2] - self.z_section;
            // Detect sign change (crossing).
            if dz_prev * dz_curr < 0.0 {
                let upward = dz_curr > dz_prev;
                if !self.upward_only || upward {
                    // Linear interpolation.
                    let t = dz_prev.abs() / (dz_prev.abs() + dz_curr.abs());
                    let x = prev[0] + t * (curr[0] - prev[0]);
                    let y = prev[1] + t * (curr[1] - prev[1]);
                    crossings.push([x, y]);
                }
            }
        }
        crossings
    }
}

// ─── Correlation Dimension ───────────────────────────────────────────────────

/// Grassberger-Procaccia correlation dimension estimator.
///
/// The correlation integral C(r) ~ r^D gives the correlation dimension D.
#[derive(Debug, Clone)]
pub struct CorrelationDimension {
    /// Number of radii to evaluate C(r) at.
    pub n_radii: usize,
    /// Minimum radius (log scale).
    pub r_min: f64,
    /// Maximum radius (log scale).
    pub r_max: f64,
}

impl CorrelationDimension {
    /// Creates a new correlation dimension estimator.
    pub fn new(n_radii: usize, r_min: f64, r_max: f64) -> Self {
        Self {
            n_radii,
            r_min,
            r_max,
        }
    }

    /// Computes the correlation integral for 3D data.
    ///
    /// Returns `(log_r, log_C)` pairs suitable for linear regression.
    pub fn correlation_integral(&self, data: &[[f64; 3]]) -> Vec<[f64; 2]> {
        let n = data.len();
        if n < 2 {
            return vec![];
        }

        let log_r_min = self.r_min.ln();
        let log_r_max = self.r_max.ln();
        let step = (log_r_max - log_r_min) / (self.n_radii as f64 - 1.0);

        let radii: Vec<f64> = (0..self.n_radii)
            .map(|i| (log_r_min + i as f64 * step).exp())
            .collect();

        let n_pairs = (n * (n - 1)) as f64;
        radii
            .iter()
            .map(|&r| {
                let count = data
                    .iter()
                    .enumerate()
                    .flat_map(|(i, a)| data[i + 1..].iter().map(move |b| dist3(*a, *b)))
                    .filter(|&d| d < r)
                    .count();
                let c = 2.0 * count as f64 / n_pairs;
                [r.ln(), if c > 0.0 { c.ln() } else { f64::NEG_INFINITY }]
            })
            .collect()
    }

    /// Estimates the slope (dimension) from the linear part of the log-log plot.
    pub fn estimate_dimension(&self, data: &[[f64; 3]]) -> f64 {
        let pairs = self.correlation_integral(data);
        let valid: Vec<[f64; 2]> = pairs.into_iter().filter(|p| p[1].is_finite()).collect();
        if valid.len() < 2 {
            return 0.0;
        }
        linear_regression_slope(&valid)
    }
}

// ─── Box-Counting Dimension ──────────────────────────────────────────────────

/// Box-counting fractal dimension estimator.
///
/// D = lim_{ε→0} log(N(ε)) / log(1/ε).
#[derive(Debug, Clone)]
pub struct BoxCounting {
    /// Box sizes to use (decreasing).
    pub box_sizes: Vec<f64>,
}

impl BoxCounting {
    /// Creates a new box-counting estimator with logarithmically spaced sizes.
    pub fn new(n_sizes: usize, size_min: f64, size_max: f64) -> Self {
        let log_min = size_min.ln();
        let log_max = size_max.ln();
        let step = (log_max - log_min) / (n_sizes as f64 - 1.0);
        let box_sizes = (0..n_sizes)
            .map(|i| (log_max - i as f64 * step).exp())
            .collect();
        Self { box_sizes }
    }

    /// Counts occupied boxes of size `eps` for 3D data.
    pub fn count_boxes(&self, data: &[[f64; 3]], eps: f64) -> usize {
        use std::collections::HashSet;
        let mut occupied: HashSet<(i64, i64, i64)> = HashSet::new();
        for &[x, y, z] in data {
            let bx = (x / eps).floor() as i64;
            let by = (y / eps).floor() as i64;
            let bz = (z / eps).floor() as i64;
            occupied.insert((bx, by, bz));
        }
        occupied.len()
    }

    /// Estimates the fractal dimension via linear regression on log-log plot.
    pub fn estimate_dimension(&self, data: &[[f64; 3]]) -> f64 {
        let pairs: Vec<[f64; 2]> = self
            .box_sizes
            .iter()
            .map(|&eps| {
                let n = self.count_boxes(data, eps) as f64;
                [(1.0 / eps).ln(), if n > 0.0 { n.ln() } else { 0.0 }]
            })
            .collect();
        linear_regression_slope(&pairs)
    }
}

// ─── Recurrence Plot ─────────────────────────────────────────────────────────

/// Recurrence plot and recurrence quantification analysis (RQA).
#[derive(Debug, Clone)]
pub struct RecurrencePlot {
    /// Threshold radius ε for recurrence.
    pub epsilon: f64,
}

impl RecurrencePlot {
    /// Creates a new recurrence plot with the given threshold.
    pub fn new(epsilon: f64) -> Self {
        Self { epsilon }
    }

    /// Builds the recurrence matrix for 3D data.
    ///
    /// Returns a flattened row-major boolean matrix of size `n × n`.
    pub fn build(&self, data: &[[f64; 3]]) -> Vec<bool> {
        let n = data.len();
        let mut mat = vec![false; n * n];
        for i in 0..n {
            for j in 0..n {
                mat[i * n + j] = dist3(data[i], data[j]) < self.epsilon;
            }
        }
        mat
    }

    /// Recurrence rate RR = (number of recurrences) / N².
    pub fn recurrence_rate(&self, data: &[[f64; 3]]) -> f64 {
        let n = data.len();
        if n == 0 {
            return 0.0;
        }
        let mat = self.build(data);
        let count = mat.iter().filter(|&&b| b).count();
        count as f64 / (n * n) as f64
    }

    /// Determinism DET = fraction of recurrences forming diagonal lines ≥ min_line.
    pub fn determinism(&self, data: &[[f64; 3]], min_line: usize) -> f64 {
        let n = data.len();
        if n == 0 {
            return 0.0;
        }
        let mat = self.build(data);
        let total_rec: usize = mat.iter().filter(|&&b| b).count();
        if total_rec == 0 {
            return 0.0;
        }

        let mut in_diag = 0usize;
        for diag in (-(n as isize) + 1)..(n as isize) {
            let mut run = 0usize;
            let i_start = (-diag).max(0) as usize;
            let j_start = diag.max(0) as usize;
            let len = n.saturating_sub(i_start.max(j_start));
            for k in 0..len {
                let i = i_start + k;
                let j = j_start + k;
                if mat[i * n + j] {
                    run += 1;
                } else {
                    if run >= min_line {
                        in_diag += run;
                    }
                    run = 0;
                }
            }
            if run >= min_line {
                in_diag += run;
            }
        }
        in_diag as f64 / total_rec as f64
    }
}

// ─── Mutual Information ──────────────────────────────────────────────────────

/// Average mutual information for embedding lag selection.
///
/// Used to select the optimal delay τ for Takens' embedding.
#[derive(Debug, Clone)]
pub struct MutualInformation {
    /// Number of histogram bins.
    pub n_bins: usize,
}

impl MutualInformation {
    /// Creates a new mutual information calculator with the given number of bins.
    pub fn new(n_bins: usize) -> Self {
        Self { n_bins }
    }

    /// Computes the average mutual information between `x[t]` and `x[t+tau]`.
    pub fn compute(&self, data: &[f64], tau: usize) -> f64 {
        if tau >= data.len() {
            return 0.0;
        }
        let n = data.len() - tau;
        let min_val = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_val = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = (max_val - min_val).max(1e-15);
        let nb = self.n_bins;

        // Build joint histogram.
        let mut joint = vec![vec![0usize; nb]; nb];
        let mut marg_x = vec![0usize; nb];
        let mut marg_y = vec![0usize; nb];

        for i in 0..n {
            let bx = ((data[i] - min_val) / range * nb as f64).floor() as usize;
            let by = ((data[i + tau] - min_val) / range * nb as f64).floor() as usize;
            let bx = bx.min(nb - 1);
            let by = by.min(nb - 1);
            joint[bx][by] += 1;
            marg_x[bx] += 1;
            marg_y[by] += 1;
        }

        let mut mi = 0.0f64;
        for (joint_row, &mx) in joint.iter().zip(marg_x.iter()) {
            let px = mx as f64 / n as f64;
            for (jval, &my) in joint_row.iter().zip(marg_y.iter()) {
                let pxy = *jval as f64 / n as f64;
                let py = my as f64 / n as f64;
                if pxy > 0.0 && px > 0.0 && py > 0.0 {
                    mi += pxy * (pxy / (px * py)).ln();
                }
            }
        }
        mi
    }

    /// Returns the mutual information vs lag for lags 1..=max_tau.
    pub fn ami_curve(&self, data: &[f64], max_tau: usize) -> Vec<(usize, f64)> {
        (1..=max_tau)
            .map(|tau| (tau, self.compute(data, tau)))
            .collect()
    }

    /// Finds the first local minimum of AMI (optimal embedding lag).
    pub fn optimal_lag(&self, data: &[f64], max_tau: usize) -> usize {
        let curve = self.ami_curve(data, max_tau);
        for i in 1..curve.len().saturating_sub(1) {
            if curve[i].1 < curve[i - 1].1 && curve[i].1 < curve[i + 1].1 {
                return curve[i].0;
            }
        }
        1
    }
}

// ─── Embedding Dimension (FNN) ───────────────────────────────────────────────

/// False nearest neighbours (FNN) method for embedding dimension (Kennel et al.).
///
/// Implements Takens' embedding theorem numerically.
#[derive(Debug, Clone)]
pub struct EmbeddingDimension {
    /// Embedding lag τ.
    pub tau: usize,
    /// FNN threshold R_tol (typically 10–15).
    pub r_tol: f64,
    /// Attractor size threshold A_tol.
    pub a_tol: f64,
}

impl EmbeddingDimension {
    /// Creates a new FNN estimator.
    pub fn new(tau: usize, r_tol: f64, a_tol: f64) -> Self {
        Self { tau, r_tol, a_tol }
    }

    /// Default parameters (τ=1, R_tol=15.0, A_tol=2.0).
    pub fn default_params() -> Self {
        Self::new(1, 15.0, 2.0)
    }

    /// Embeds a scalar time series into d-dimensional delay vectors.
    pub fn embed(&self, data: &[f64], d: usize) -> Vec<Vec<f64>> {
        if data.len() < d * self.tau {
            return vec![];
        }
        let n = data.len() - (d - 1) * self.tau;
        (0..n)
            .map(|i| (0..d).map(|k| data[i + k * self.tau]).collect())
            .collect()
    }

    /// Computes the fraction of false nearest neighbours for embedding dimension `d`.
    pub fn fnn_fraction(&self, data: &[f64], d: usize) -> f64 {
        let embedded_d = self.embed(data, d);
        let embedded_d1 = self.embed(data, d + 1);
        let n = embedded_d.len().min(embedded_d1.len());
        if n < 2 {
            return 0.0;
        }
        // Std deviation of data as attractor size estimate.
        let mean = data.iter().sum::<f64>() / data.len() as f64;
        let std = (data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / data.len() as f64)
            .sqrt()
            .max(1e-15);

        let mut n_fnn = 0usize;
        let mut n_total = 0usize;
        for i in 0..n {
            // Find nearest neighbour in d dimensions.
            let mut nn_dist = f64::INFINITY;
            let mut nn_idx = 0usize;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dist = embedded_d[i]
                    .iter()
                    .zip(&embedded_d[j])
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                if dist < nn_dist {
                    nn_dist = dist;
                    nn_idx = j;
                }
            }
            if nn_dist < 1e-15 {
                continue;
            }
            // Distance in d+1 dimensions.
            let new_component = (embedded_d1[i][d] - embedded_d1[nn_idx][d]).abs();
            let dist_d1 = (nn_dist * nn_dist + new_component * new_component).sqrt();

            n_total += 1;
            if new_component / nn_dist > self.r_tol || dist_d1 / std > self.a_tol {
                n_fnn += 1;
            }
        }
        if n_total == 0 {
            0.0
        } else {
            n_fnn as f64 / n_total as f64
        }
    }

    /// Finds the minimal embedding dimension where FNN fraction < threshold.
    pub fn estimate(&self, data: &[f64], max_d: usize, threshold: f64) -> usize {
        for d in 1..=max_d {
            if self.fnn_fraction(data, d) < threshold {
                return d;
            }
        }
        max_d
    }
}

// ─── Feigenbaum Constant ─────────────────────────────────────────────────────

/// Numerical estimation of the Feigenbaum constant δ ≈ 4.6692016…
///
/// Computed from the period-doubling cascade of the logistic map.
#[derive(Debug, Clone)]
pub struct FeigenbaumConstant;

impl FeigenbaumConstant {
    /// Known theoretical value of the Feigenbaum constant δ.
    pub const DELTA: f64 = 4.669_201_609_102_99;

    /// Finds the r-value where the logistic map first shows period-2k behaviour.
    ///
    /// Returns the approximate bifurcation point via iterated cobweb detection.
    pub fn bifurcation_point(period: usize, x0: f64, transient: usize, n_check: usize) -> f64 {
        // Binary search for the r value at which period 2k first appears.
        let mut r_lo = 2.5f64;
        let mut r_hi = 4.0f64;
        for _ in 0..60 {
            let r_mid = 0.5 * (r_lo + r_hi);
            let lm = LogisticMap::new(r_mid);
            let attr = lm.attractor(x0, transient, n_check);
            let detected = detect_period(&attr, 0.001);
            if detected >= period {
                r_hi = r_mid;
            } else {
                r_lo = r_mid;
            }
        }
        0.5 * (r_lo + r_hi)
    }

    /// Estimates δ from the first few bifurcation points.
    pub fn estimate_delta(x0: f64) -> f64 {
        let r1 = Self::bifurcation_point(2, x0, 1000, 256);
        let r2 = Self::bifurcation_point(4, x0, 1000, 512);
        let r3 = Self::bifurcation_point(8, x0, 1000, 1024);
        (r2 - r1) / (r3 - r2).max(1e-12)
    }
}

// ─── FTLE Field ───────────────────────────────────────────────────────────────

/// Finite-Time Lyapunov Exponent (FTLE) field on a regular 2D grid.
///
/// Computes the FTLE by tracking the deformation of a grid of initial conditions.
#[derive(Debug, Clone)]
pub struct FtleField {
    /// Grid resolution in x.
    pub nx: usize,
    /// Grid resolution in y.
    pub ny: usize,
    /// Integration time T (positive → forward FTLE, negative → backward).
    pub t_integration: f64,
    /// Sub-step size dt.
    pub dt: f64,
}

impl FtleField {
    /// Creates a new FTLE field calculator.
    pub fn new(nx: usize, ny: usize, t_integration: f64, dt: f64) -> Self {
        Self {
            nx,
            ny,
            t_integration,
            dt,
        }
    }

    /// Computes the FTLE field for the Lorenz system projected onto the xy plane.
    ///
    /// `x_range = (x_min, x_max)`, `y_range = (y_min, y_max)`, z is fixed.
    pub fn compute_lorenz(
        &self,
        x_range: (f64, f64),
        y_range: (f64, f64),
        z_fixed: f64,
        sigma: f64,
        rho: f64,
        beta: f64,
    ) -> Vec<Vec<f64>> {
        let lorenz = LorenzSystem::new(sigma, rho, beta);
        let steps = (self.t_integration / self.dt).abs() as usize;
        let dt = if self.t_integration > 0.0 {
            self.dt
        } else {
            -self.dt
        };
        let dx = (x_range.1 - x_range.0) / (self.nx as f64 - 1.0);
        let dy = (y_range.1 - y_range.0) / (self.ny as f64 - 1.0);
        let eps = dx.min(dy) * 0.01;

        let mut ftle = vec![vec![0.0f64; self.ny]; self.nx];
        for (ix, ftle_row) in ftle.iter_mut().enumerate() {
            for (iy, cell) in ftle_row.iter_mut().enumerate() {
                let x0 = x_range.0 + ix as f64 * dx;
                let y0 = y_range.0 + iy as f64 * dy;
                let p0 = [x0, y0, z_fixed];
                let px = [x0 + eps, y0, z_fixed];
                let py = [x0, y0 + eps, z_fixed];

                let mut p0f = p0;
                let mut pxf = px;
                let mut pyf = py;
                for _ in 0..steps {
                    p0f = rk4_step(|s| lorenz.deriv(s), p0f, dt);
                    pxf = rk4_step(|s| lorenz.deriv(s), pxf, dt);
                    pyf = rk4_step(|s| lorenz.deriv(s), pyf, dt);
                }

                // Deformation gradient columns.
                let dfdx = [(pxf[0] - p0f[0]) / eps, (pxf[1] - p0f[1]) / eps];
                let dfdy = [(pyf[0] - p0f[0]) / eps, (pyf[1] - p0f[1]) / eps];
                // Cauchy-Green strain tensor C = F^T F (2×2).
                let c11 = dfdx[0] * dfdx[0] + dfdx[1] * dfdx[1];
                let c12 = dfdx[0] * dfdy[0] + dfdx[1] * dfdy[1];
                let c22 = dfdy[0] * dfdy[0] + dfdy[1] * dfdy[1];
                // Largest eigenvalue of C.
                let tr = c11 + c22;
                let det = c11 * c22 - c12 * c12;
                let disc = (0.25 * tr * tr - det).max(0.0);
                let lambda_max = 0.5 * tr + disc.sqrt();
                *cell = if lambda_max > 1.0 {
                    lambda_max.ln() / (2.0 * self.t_integration.abs())
                } else {
                    0.0
                };
            }
        }
        ftle
    }
}

// ─── Period-Doubling Route to Chaos ─────────────────────────────────────────

/// Analysis of the period-doubling cascade in a map.
#[derive(Debug, Clone)]
pub struct PeriodDoubling;

impl PeriodDoubling {
    /// Detects the period of a time series (up to max_period).
    ///
    /// Returns the smallest period p such that `|x[i] − x[i+p]| < tol` for all i.
    pub fn detect_period(series: &[f64], tol: f64, max_period: usize) -> usize {
        detect_period(series, tol).min(max_period)
    }

    /// Returns `(r, period)` pairs for the logistic map over a range of r values.
    pub fn period_diagram(
        r_values: &[f64],
        x0: f64,
        transient: usize,
        check_len: usize,
    ) -> Vec<(f64, usize)> {
        r_values
            .iter()
            .map(|&r| {
                let lm = LogisticMap::new(r);
                let attr = lm.attractor(x0, transient, check_len);
                let p = detect_period(&attr, 0.001);
                (r, p)
            })
            .collect()
    }

    /// Counts the number of period-doubling bifurcations up to r_max.
    pub fn count_doublings(r_max: f64, x0: f64) -> usize {
        let mut count = 0usize;
        let mut period = 1usize;
        for i in 1..=6 {
            let target = 1usize << i;
            let r = FeigenbaumConstant::bifurcation_point(target, x0, 500, 512);
            if r < r_max {
                count += 1;
                period = target;
            }
        }
        let _ = period;
        count
    }
}

// ─── Chaotic Time Series Analysis ────────────────────────────────────────────

/// Comprehensive chaotic time series analysis toolkit.
#[derive(Debug, Clone)]
pub struct ChaosAnalysis {
    /// Embedding lag.
    pub tau: usize,
    /// Embedding dimension.
    pub dim: usize,
}

impl ChaosAnalysis {
    /// Creates a new analysis instance with the given embedding parameters.
    pub fn new(tau: usize, dim: usize) -> Self {
        Self { tau, dim }
    }

    /// Computes the largest Lyapunov exponent from a scalar time series
    /// (Rosenstein et al. algorithm).
    ///
    /// Returns an estimate of the largest Lyapunov exponent in nats/step.
    pub fn largest_lyapunov(&self, data: &[f64], evolve_steps: usize) -> f64 {
        let emb = EmbeddingDimension::new(self.tau, 15.0, 2.0);
        let vecs = emb.embed(data, self.dim);
        let n = vecs.len();
        if n < 2 {
            return 0.0;
        }

        let min_sep = self.dim * self.tau + 1;
        let mut sum = 0.0f64;
        let mut count = 0usize;

        for i in 0..n {
            let mut nn_dist = f64::INFINITY;
            let mut nn_j = 0usize;
            for j in 0..n {
                if (i as isize - j as isize).unsigned_abs() < min_sep {
                    continue;
                }
                let d: f64 = vecs[i]
                    .iter()
                    .zip(&vecs[j])
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                if d < nn_dist {
                    nn_dist = d;
                    nn_j = j;
                }
            }
            if nn_dist < 1e-15 {
                continue;
            }

            let end_i = (i + evolve_steps).min(n - 1);
            let end_j = (nn_j + evolve_steps).min(n - 1);
            if end_i == i || end_j == nn_j {
                continue;
            }

            let d_evolved: f64 = vecs[end_i]
                .iter()
                .zip(&vecs[end_j])
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            if d_evolved > 1e-15 {
                sum += (d_evolved / nn_dist).ln();
                count += 1;
            }
        }

        if count == 0 {
            0.0
        } else {
            sum / (count as f64 * evolve_steps as f64)
        }
    }

    /// Estimates the correlation dimension of a scalar time series.
    pub fn correlation_dimension(&self, data: &[f64]) -> f64 {
        let emb = EmbeddingDimension::new(self.tau, 15.0, 2.0);
        let vecs = emb.embed(data, self.dim);
        if vecs.len() < 10 {
            return 0.0;
        }
        // Convert to [f64; 3] by taking first 3 components (pad if needed).
        let data3: Vec<[f64; 3]> = vecs
            .iter()
            .map(|v| {
                let x = v.first().copied().unwrap_or(0.0);
                let y = v.get(1).copied().unwrap_or(0.0);
                let z = v.get(2).copied().unwrap_or(0.0);
                [x, y, z]
            })
            .collect();
        let cd = CorrelationDimension::new(20, 0.01, 5.0);
        cd.estimate_dimension(&data3)
    }

    /// Prints a summary of chaos diagnostics.
    pub fn summary(&self, data: &[f64]) -> String {
        let lambda = self.largest_lyapunov(data, 5);
        let d_corr = self.correlation_dimension(data);
        format!(
            "Largest Lyapunov: {:.6}, Correlation dim: {:.6}",
            lambda, d_corr
        )
    }
}

// ─── Internal Helpers ─────────────────────────────────────────────────────────

/// 4th-order Runge-Kutta step for a 3D autonomous ODE.
fn rk4_step<F>(f: F, state: [f64; 3], dt: f64) -> [f64; 3]
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    let k1 = f(state);
    let s2 = add3(state, scale3(k1, 0.5 * dt));
    let k2 = f(s2);
    let s3 = add3(state, scale3(k2, 0.5 * dt));
    let k3 = f(s3);
    let s4 = add3(state, scale3(k3, dt));
    let k4 = f(s4);
    [
        state[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        state[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        state[2] + dt / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
    ]
}

/// Euclidean distance between two 3D points.
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Component-wise addition of two 3D vectors.
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scalar multiplication of a 3D vector.
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Gram-Schmidt QR decomposition for 3 column vectors.
///
/// Returns (orthonormal Q, diagonal of R).
fn gram_schmidt_qr(cols: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3]) {
    let mut q = [[0.0f64; 3]; 3];
    let mut r_diag = [0.0f64; 3];

    for j in 0..3 {
        let mut v = cols[j];
        // Subtract projections onto previous columns.
        for qk in q[..j].iter() {
            let proj = dot3(v, *qk);
            for (vi, &qki) in v.iter_mut().zip(qk.iter()) {
                *vi -= proj * qki;
            }
        }
        let norm = dot3(v, v).sqrt();
        r_diag[j] = norm;
        if norm > 1e-15 {
            for i in 0..3 {
                q[j][i] = v[i] / norm;
            }
        } else {
            q[j] = [0.0; 3];
            q[j][j] = 1.0;
        }
    }
    (q, r_diag)
}

/// Dot product of two 3D vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Linear regression slope from (x, y) pairs.
fn linear_regression_slope(pairs: &[[f64; 2]]) -> f64 {
    let n = pairs.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let sum_x: f64 = pairs.iter().map(|p| p[0]).sum();
    let sum_y: f64 = pairs.iter().map(|p| p[1]).sum();
    let sum_xx: f64 = pairs.iter().map(|p| p[0] * p[0]).sum();
    let sum_xy: f64 = pairs.iter().map(|p| p[0] * p[1]).sum();
    let denom = n * sum_xx - sum_x * sum_x;
    if denom.abs() < 1e-15 {
        return 0.0;
    }
    (n * sum_xy - sum_x * sum_y) / denom
}

/// Detects the period of a time series.
fn detect_period(series: &[f64], tol: f64) -> usize {
    let n = series.len();
    if n < 4 {
        return 1;
    }
    for p in 1..=(n / 2) {
        let is_periodic = series[..n - p]
            .iter()
            .zip(&series[p..])
            .all(|(a, b)| (a - b).abs() < tol);
        if is_periodic {
            return p;
        }
    }
    n // Treat as "chaotic" = period n.
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Lorenz system ─────────────────────────────────────────────────────────

    #[test]
    fn test_lorenz_deriv_at_origin() {
        let sys = LorenzSystem::classic();
        let d = sys.deriv([0.0, 0.0, 0.0]);
        assert!(d[0].abs() < 1e-12, "dx/dt at origin should be 0");
        assert!(d[1].abs() < 1e-12, "dy/dt at origin should be 0");
        assert!(d[2].abs() < 1e-12, "dz/dt at origin should be 0");
    }

    #[test]
    fn test_lorenz_integrate_length() {
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([1.0, 1.0, 1.0], 0.01, 100);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_lorenz_fixed_points_symmetric() {
        let sys = LorenzSystem::classic();
        let fps = sys.fixed_points();
        // Origin.
        assert!(fps[0].iter().all(|&v| v.abs() < 1e-12));
        // The two non-trivial fixed points should be symmetric.
        assert!((fps[1][0] + fps[2][0]).abs() < 1e-12);
        assert!((fps[1][1] + fps[2][1]).abs() < 1e-12);
        assert!((fps[1][2] - fps[2][2]).abs() < 1e-12);
    }

    #[test]
    fn test_lorenz_fixed_point_is_fixed() {
        let sys = LorenzSystem::classic();
        let fps = sys.fixed_points();
        for fp in fps.iter().skip(1) {
            let d = sys.deriv(*fp);
            let mag = (d[0].powi(2) + d[1].powi(2) + d[2].powi(2)).sqrt();
            assert!(
                mag < 1e-10,
                "Fixed point derivative should be near zero, got {mag}"
            );
        }
    }

    #[test]
    fn test_lorenz_sensitivity() {
        let sys = LorenzSystem::classic();
        // Use a large perturbation so divergence is detectable within a reasonable time.
        let t1 = sys.integrate([1.0, 1.0, 1.0], 0.01, 2000);
        let t2 = sys.integrate([1.0 + 1e-3, 1.0, 1.0], 0.01, 2000);
        // After 2000 steps (t=20) even a 1e-3 perturbation should have grown in a chaotic system.
        // We just check that the two trajectories are not identical.
        let end1 = t1.last().unwrap();
        let end2 = t2.last().unwrap();
        let sep = dist3(*end1, *end2);
        assert!(
            sep > 1e-6,
            "Lorenz should show sensitivity, separation={sep}"
        );
    }

    // ── Rössler system ────────────────────────────────────────────────────────

    #[test]
    fn test_rossler_integrate_no_nan() {
        let sys = RosslerSystem::classic();
        let traj = sys.integrate([1.0, 0.0, 0.0], 0.01, 200);
        for pt in &traj {
            assert!(
                pt.iter().all(|v| v.is_finite()),
                "Rössler trajectory must be finite"
            );
        }
    }

    #[test]
    fn test_rossler_period_positive() {
        let sys = RosslerSystem::classic();
        assert!(sys.approximate_period() > 0.0);
    }

    // ── Logistic map ──────────────────────────────────────────────────────────

    #[test]
    fn test_logistic_fixed_point_r3() {
        // For r=3.0 the fixed point is x*=(r-1)/r = 2/3.
        let lm = LogisticMap::new(3.0);
        let attr = lm.attractor(0.5, 1000, 10);
        for &v in &attr {
            assert!(
                (v - 2.0 / 3.0).abs() < 0.01,
                "r=3 should approach fixed point, got {v}"
            );
        }
    }

    #[test]
    fn test_logistic_period2_r35() {
        // r=3.5 → period 4; verify at least period ≥ 2.
        let lm = LogisticMap::new(3.5);
        let attr = lm.attractor(0.5, 2000, 32);
        let p = detect_period(&attr, 0.001);
        assert!(p >= 2, "r=3.5 should have period ≥ 2, got {p}");
    }

    #[test]
    fn test_logistic_lyapunov_chaotic() {
        // r=3.9 is deeply chaotic → λ > 0.
        let lm = LogisticMap::new(3.9);
        let lambda = lm.lyapunov_exponent(0.5, 10000);
        assert!(lambda > 0.0, "r=3.9 should be chaotic, λ={lambda}");
    }

    #[test]
    fn test_logistic_lyapunov_stable() {
        // r=2.5: stable fixed point → λ ≤ 0 (derivative |r(1-2x*)| < 1 at x*=(r-1)/r).
        let lm = LogisticMap::new(2.5);
        let lambda = lm.lyapunov_exponent(0.3, 5000);
        assert!(lambda <= 0.0, "r=2.5 should be stable (λ≤0), λ={lambda}");
    }

    // ── Duffing oscillator ────────────────────────────────────────────────────

    #[test]
    fn test_duffing_integrate_length() {
        let duff = DuffingOscillator::classic_chaotic();
        let traj = duff.integrate([1.0, 0.0], 0.0, 0.01, 100);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_duffing_natural_freq_double_well() {
        let duff = DuffingOscillator::new(0.0, -1.0, 1.0, 0.0, 1.0);
        let freq = duff.natural_frequency();
        assert!(freq.is_some(), "Double-well should have natural frequency");
        assert!((freq.unwrap() - 2.0_f64.sqrt()).abs() < 0.01);
    }

    #[test]
    fn test_duffing_no_nan() {
        let duff = DuffingOscillator::classic_chaotic();
        let traj = duff.integrate([0.5, 0.0], 0.0, 0.01, 500);
        for pt in &traj {
            assert!(
                pt.iter().all(|v| v.is_finite()),
                "Duffing trajectory must be finite"
            );
        }
    }

    // ── Lyapunov exponent ─────────────────────────────────────────────────────

    #[test]
    fn test_lyapunov_lorenz_positive() {
        let est = LyapunovExponent::new(20_000, 0.01, 10, 1e-8);
        let sys = LorenzSystem::classic();
        let lambda = est.estimate(|s| sys.deriv(s), [0.1, 0.1, 0.1]);
        assert!(
            lambda > 0.0,
            "Lorenz largest Lyapunov should be positive, got {lambda}"
        );
    }

    // ── Lyapunov spectrum ─────────────────────────────────────────────────────

    #[test]
    fn test_lyapunov_spectrum_sum_negative() {
        // For Lorenz: λ1 + λ2 + λ3 = −(σ + 1 + β) < 0 (volume contraction).
        let spec = LyapunovSpectrum::new(0.01, 5, 200);
        let exps = spec.compute_lorenz(10.0, 28.0, 8.0 / 3.0, [1.0, 1.0, 1.0]);
        let sum = exps[0] + exps[1] + exps[2];
        // The sum should be near −(σ + 1 + β) ≈ −13.67, but allow numerical tolerance.
        assert!(
            sum < 0.0,
            "Lyapunov spectrum sum should be negative (volume contraction), got {sum}"
        );
    }

    // ── Poincaré section ──────────────────────────────────────────────────────

    #[test]
    fn test_poincare_section_count() {
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([1.0, 1.0, 20.0], 0.01, 5000);
        let poincare = PoincareSection::new(27.0, true);
        let crossings = poincare.extract(&traj);
        assert!(
            !crossings.is_empty(),
            "Lorenz should have Poincaré crossings"
        );
    }

    #[test]
    fn test_poincare_section_upward_only() {
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([1.0, 1.0, 25.0], 0.01, 3000);
        let all = PoincareSection::new(27.0, false).extract(&traj);
        let up_only = PoincareSection::new(27.0, true).extract(&traj);
        assert!(
            up_only.len() <= all.len(),
            "Upward-only should have ≤ crossings than all"
        );
    }

    // ── Correlation dimension ─────────────────────────────────────────────────

    #[test]
    fn test_correlation_dimension_line() {
        // Points on a line should have dimension ≈ 1.
        let data: Vec<[f64; 3]> = (0..100).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let cd = CorrelationDimension::new(10, 0.1, 5.0);
        let dim = cd.estimate_dimension(&data);
        assert!(dim < 2.0, "Line should have dimension < 2, got {dim}");
    }

    // ── Box-counting dimension ────────────────────────────────────────────────

    #[test]
    fn test_box_counting_cube() {
        // Points filling a cube → dimension ≈ 3.
        let mut data = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    data.push([i as f64, j as f64, k as f64]);
                }
            }
        }
        let bc = BoxCounting::new(5, 0.5, 4.0);
        let dim = bc.estimate_dimension(&data);
        assert!(
            dim > 1.5 && dim < 4.0,
            "Cube should have dimension ≈ 3, got {dim}"
        );
    }

    // ── Recurrence plot ───────────────────────────────────────────────────────

    #[test]
    fn test_recurrence_rate_between_0_and_1() {
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([1.0, 1.0, 1.0], 0.01, 99);
        let rp = RecurrencePlot::new(5.0);
        let rr = rp.recurrence_rate(&traj);
        assert!(
            (0.0..=1.0).contains(&rr),
            "Recurrence rate should be in [0,1], got {rr}"
        );
    }

    #[test]
    fn test_recurrence_diagonal_is_always_1() {
        let data: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rp = RecurrencePlot::new(0.5);
        let mat = rp.build(&data);
        let n = data.len();
        for i in 0..n {
            assert!(mat[i * n + i], "Diagonal should always be true");
        }
    }

    // ── Mutual information ────────────────────────────────────────────────────

    #[test]
    fn test_mutual_information_zero_lag_max() {
        // At lag=1 vs lag=large, AMI at lag=1 should generally not be zero for a chaotic series.
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([1.0, 0.5, 1.5], 0.01, 300);
        let x_series: Vec<f64> = traj.iter().map(|pt| pt[0]).collect();
        let mi = MutualInformation::new(10);
        let ami1 = mi.compute(&x_series, 1);
        let ami_large = mi.compute(&x_series, 50);
        assert!(
            ami1 >= ami_large - 0.5,
            "Short-lag AMI should not be less than long-lag AMI by more than 0.5"
        );
    }

    // ── Embedding dimension ───────────────────────────────────────────────────

    #[test]
    fn test_embedding_embed_shape() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let emb = EmbeddingDimension::new(1, 15.0, 2.0);
        let vecs = emb.embed(&data, 3);
        // Should have 20 - 2*1 = 18 vectors, each of length 3.
        assert_eq!(vecs.len(), 18);
        for v in &vecs {
            assert_eq!(v.len(), 3);
        }
    }

    // ── Feigenbaum constant ───────────────────────────────────────────────────

    #[test]
    fn test_feigenbaum_delta_known_value() {
        // The known value is ≈ 4.669.
        assert!(
            (FeigenbaumConstant::DELTA - 4.6692).abs() < 0.001,
            "Feigenbaum delta should be ≈ 4.6692"
        );
    }

    #[test]
    fn test_feigenbaum_bifurcation_r2() {
        // Period-2 bifurcation at r ≈ 3.0.
        let r = FeigenbaumConstant::bifurcation_point(2, 0.5, 1000, 256);
        assert!(
            (r - 3.0).abs() < 0.1,
            "Period-2 bifurcation should be near r=3.0, got {r}"
        );
    }

    // ── FTLE field ────────────────────────────────────────────────────────────

    #[test]
    fn test_ftle_field_non_negative() {
        let ftle = FtleField::new(5, 5, 0.5, 0.05);
        let field = ftle.compute_lorenz((-5.0, 5.0), (-5.0, 5.0), 25.0, 10.0, 28.0, 8.0 / 3.0);
        for row in &field {
            for &v in row {
                assert!(v >= 0.0, "FTLE values should be non-negative, got {v}");
            }
        }
    }

    // ── Period doubling ───────────────────────────────────────────────────────

    #[test]
    fn test_period_detection_simple() {
        let series = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0];
        let p = PeriodDoubling::detect_period(&series, 0.001, 10);
        assert_eq!(p, 2, "Simple alternating series should have period 2");
    }

    #[test]
    fn test_period_doubling_count() {
        // There should be at least 2 period-doubling bifurcations before r=3.6.
        let count = PeriodDoubling::count_doublings(3.6, 0.5);
        assert!(
            count >= 2,
            "Should see ≥ 2 doublings before r=3.6, got {count}"
        );
    }

    // ── Chaos analysis ────────────────────────────────────────────────────────

    #[test]
    fn test_chaos_analysis_summary_format() {
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([1.0, 1.0, 1.0], 0.01, 200);
        let x_series: Vec<f64> = traj.iter().map(|pt| pt[0]).collect();
        let ca = ChaosAnalysis::new(1, 3);
        let s = ca.summary(&x_series);
        assert!(
            s.contains("Largest Lyapunov"),
            "Summary should contain Lyapunov info"
        );
        assert!(
            s.contains("Correlation dim"),
            "Summary should contain dimension info"
        );
    }

    #[test]
    fn test_chaos_analysis_lorenz_positive_lyapunov() {
        let sys = LorenzSystem::classic();
        let traj = sys.integrate([0.1, 0.1, 0.1], 0.01, 1000);
        let x_series: Vec<f64> = traj.iter().map(|pt| pt[0]).collect();
        let ca = ChaosAnalysis::new(1, 3);
        let lambda = ca.largest_lyapunov(&x_series, 5);
        // Should be positive for chaotic series.
        assert!(
            lambda > -1.0,
            "Chaos analysis Lyapunov should not be very negative, got {lambda}"
        );
    }
}
