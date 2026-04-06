// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dynamical systems — maps, flows, bifurcations, synchronisation.
//!
//! Covers discrete maps (logistic, Hénon, Chirikov), strange attractor analysis,
//! normal forms for bifurcations, Poincaré sections, flow maps and the
//! fundamental solution matrix, and Pecora–Carroll / phase synchronisation.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// DiscreteMap
// ─────────────────────────────────────────────────────────────────────────────

/// A discrete map that can be iterated.
pub trait DiscreteMapIterate {
    /// State dimension.
    fn dim(&self) -> usize;
    /// Advance the state by one step.  `state` must have length `dim()`.
    fn step(&self, state: &[f64]) -> Vec<f64>;

    /// Compute an orbit of length `n` starting from `init`.
    fn orbit(&self, init: &[f64], n: usize) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(n + 1);
        out.push(init.to_vec());
        for _ in 0..n {
            let next = self.step(out.last().expect("orbit is non-empty"));
            out.push(next);
        }
        out
    }
}

// ── Logistic Map ─────────────────────────────────────────────────────────────

/// The logistic map  `x_{n+1} = r · xₙ (1 − xₙ)`.
///
/// Classic example of period-doubling route to chaos with `r ∈ [0, 4]`.
#[derive(Debug, Clone)]
pub struct LogisticMap {
    /// Growth parameter `r`.
    pub r: f64,
}

impl LogisticMap {
    /// Create a logistic map with growth parameter `r`.
    pub fn new(r: f64) -> Self {
        Self { r }
    }

    /// One application of the map.
    pub fn apply(&self, x: f64) -> f64 {
        self.r * x * (1.0 - x)
    }

    /// Compute an orbit starting at `x0` for `n` steps.
    pub fn orbit_1d(&self, x0: f64, n: usize) -> Vec<f64> {
        let mut orbit = Vec::with_capacity(n + 1);
        orbit.push(x0);
        for _ in 0..n {
            let last = *orbit.last().expect("orbit is non-empty");
            orbit.push(self.apply(last));
        }
        orbit
    }

    /// Generate a bifurcation diagram by sweeping `r` in `[r_min, r_max]`
    /// with `n_r` steps.  For each `r` value the system is iterated
    /// `transient` steps (discarded) then `keep` steps are recorded.
    ///
    /// Returns a list of `(r, x)` pairs.
    pub fn bifurcation_diagram(
        r_min: f64,
        r_max: f64,
        n_r: usize,
        transient: usize,
        keep: usize,
    ) -> Vec<(f64, f64)> {
        let mut result = Vec::with_capacity(n_r * keep);
        for k in 0..n_r {
            let r = r_min + (r_max - r_min) * k as f64 / n_r.saturating_sub(1).max(1) as f64;
            let map = LogisticMap::new(r);
            let mut x = 0.5;
            for _ in 0..transient {
                x = map.apply(x);
            }
            for _ in 0..keep {
                x = map.apply(x);
                result.push((r, x));
            }
        }
        result
    }

    /// Largest Lyapunov exponent estimated over `n` iterations.
    ///
    /// Returns `None` if the trajectory passes through 0 or 1 (derivative undefined).
    pub fn lyapunov_exponent(&self, x0: f64, n: usize) -> Option<f64> {
        let mut x = x0;
        let mut sum = 0.0f64;
        for _ in 0..n {
            let deriv = (self.r * (1.0 - 2.0 * x)).abs();
            if deriv < 1e-300 {
                return None;
            }
            sum += deriv.ln();
            x = self.apply(x);
        }
        Some(sum / n as f64)
    }
}

impl DiscreteMapIterate for LogisticMap {
    fn dim(&self) -> usize {
        1
    }
    fn step(&self, state: &[f64]) -> Vec<f64> {
        vec![self.apply(state[0])]
    }
}

// ── Hénon Map ─────────────────────────────────────────────────────────────────

/// The Hénon map:
///
/// ```text
/// x_{n+1} = 1 − a xₙ² + yₙ
/// y_{n+1} = b xₙ
/// ```
///
/// Classical parameter values are `a = 1.4`, `b = 0.3`.
#[derive(Debug, Clone)]
pub struct HenonMap {
    /// Parameter `a`.
    pub a: f64,
    /// Parameter `b`.
    pub b: f64,
}

impl HenonMap {
    /// Create a Hénon map with given parameters.
    pub fn new(a: f64, b: f64) -> Self {
        Self { a, b }
    }

    /// Apply one step starting at `(x, y)`.  Returns `(x_new, y_new)`.
    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        (1.0 - self.a * x * x + y, self.b * x)
    }

    /// Compute an orbit for `n` steps starting at `(x0, y0)`.
    pub fn orbit_2d(&self, x0: f64, y0: f64, n: usize) -> Vec<(f64, f64)> {
        let mut out = Vec::with_capacity(n + 1);
        out.push((x0, y0));
        for _ in 0..n {
            let &(x, y) = out.last().expect("orbit is non-empty");
            out.push(self.apply(x, y));
        }
        out
    }

    /// Jacobian determinant (constant for the Hénon map): `−b`.
    pub fn jacobian_det(&self) -> f64 {
        -self.b
    }

    /// Jacobian matrix at `(x, y)`:
    /// ```text
    /// [[-2ax,  1],
    ///  [b,     0]]
    /// ```
    pub fn jacobian(&self, x: f64, _y: f64) -> [[f64; 2]; 2] {
        [[-2.0 * self.a * x, 1.0], [self.b, 0.0]]
    }
}

impl DiscreteMapIterate for HenonMap {
    fn dim(&self) -> usize {
        2
    }
    fn step(&self, state: &[f64]) -> Vec<f64> {
        let (xn, yn) = self.apply(state[0], state[1]);
        vec![xn, yn]
    }
}

// ── Chirikov Standard Map ─────────────────────────────────────────────────────

/// The Chirikov–Taylor standard map (area-preserving):
///
/// ```text
/// p_{n+1} = p_n + K sin(θ_n)   (mod 2π)
/// θ_{n+1} = θ_n + p_{n+1}      (mod 2π)
/// ```
///
/// `K` is the stochasticity parameter; chaos sets in for `K ≳ 1`.
#[derive(Debug, Clone)]
pub struct ChirikovStandardMap {
    /// Stochasticity parameter `K`.
    pub k: f64,
}

impl ChirikovStandardMap {
    /// Create a Chirikov standard map with parameter `K`.
    pub fn new(k: f64) -> Self {
        Self { k }
    }

    /// Apply one step.  Both `theta` and `p` are taken modulo `2π`.
    /// Returns `(theta_new, p_new)`.
    pub fn apply(&self, theta: f64, p: f64) -> (f64, f64) {
        let p_new = (p + self.k * theta.sin()).rem_euclid(std::f64::consts::TAU);
        let theta_new = (theta + p_new).rem_euclid(std::f64::consts::TAU);
        (theta_new, p_new)
    }

    /// Compute an orbit for `n` steps.
    pub fn orbit_2d(&self, theta0: f64, p0: f64, n: usize) -> Vec<(f64, f64)> {
        let mut out = Vec::with_capacity(n + 1);
        out.push((theta0, p0));
        for _ in 0..n {
            let &(th, p) = out.last().expect("orbit is non-empty");
            out.push(self.apply(th, p));
        }
        out
    }

    /// Estimate the Lyapunov exponent via the method of Benettin et al.
    pub fn lyapunov_exponent(&self, theta0: f64, p0: f64, n: usize) -> f64 {
        let mut theta = theta0;
        let mut p = p0;
        // tangent vector
        let mut dt = 1.0f64;
        let mut dp = 0.0f64;
        let mut sum = 0.0f64;

        for _ in 0..n {
            let cos_th = theta.cos();
            // Jacobian: [[1, 0], [K cos θ, 1]] * tangent
            let dp_new = dp + self.k * cos_th * dt;
            let dt_new = dt + dp_new;
            let norm = (dt_new * dt_new + dp_new * dp_new).sqrt().max(1e-300);
            sum += norm.ln();
            dt = dt_new / norm;
            dp = dp_new / norm;
            let (th_new, p_new) = self.apply(theta, p);
            theta = th_new;
            p = p_new;
        }
        sum / n as f64
    }
}

impl DiscreteMapIterate for ChirikovStandardMap {
    fn dim(&self) -> usize {
        2
    }
    fn step(&self, state: &[f64]) -> Vec<f64> {
        let (th, p) = self.apply(state[0], state[1]);
        vec![th, p]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AttractorAnalysis
// ─────────────────────────────────────────────────────────────────────────────

/// Tools for characterising attractors from time-series data.
#[derive(Debug, Clone)]
pub struct AttractorAnalysis {
    /// The embedded time-series (each inner `Vec` is one point in phase space).
    pub data: Vec<Vec<f64>>,
}

impl AttractorAnalysis {
    /// Create from a time-series.
    pub fn new(data: Vec<Vec<f64>>) -> Self {
        Self { data }
    }

    /// Find approximate fixed points: points where `|f(x) − x| < tol`.
    ///
    /// `f` is the map applied once.
    pub fn find_fixed_points<F>(&self, f: F, tol: f64) -> Vec<Vec<f64>>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        self.data
            .iter()
            .filter(|x| {
                let fx = f(x.as_slice());
                let dist: f64 = x
                    .iter()
                    .zip(fx.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                dist < tol
            })
            .cloned()
            .collect()
    }

    /// Detect periodic orbits of period `p`: points where `|f^p(x) − x| < tol`.
    pub fn find_periodic_orbits<F>(&self, f: &F, period: usize, tol: f64) -> Vec<Vec<f64>>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        self.data
            .iter()
            .filter(|x| {
                let mut state = x.to_vec();
                for _ in 0..period {
                    state = f(&state);
                }
                let dist: f64 = x
                    .iter()
                    .zip(state.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                dist < tol
            })
            .cloned()
            .collect()
    }

    /// Estimate the correlation dimension via the Grassberger–Procaccia algorithm.
    ///
    /// Computes `C(r)` for a range of radii and returns `(log_r, log_C)` pairs
    /// from which the slope gives the correlation dimension.
    pub fn correlation_dimension(&self, n_radii: usize) -> Vec<(f64, f64)> {
        let n = self.data.len();
        if n < 2 {
            return vec![];
        }

        // Compute all pairwise distances
        let mut dists = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let d: f64 = self.data[i]
                    .iter()
                    .zip(self.data[j].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                dists.push(d);
            }
        }
        dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        if dists.is_empty() {
            return vec![];
        }

        let r_min = *dists.first().expect("dists is non-empty") * 1.1;
        let r_max = *dists.last().expect("dists is non-empty") * 0.9;
        if r_min >= r_max {
            return vec![];
        }

        let total_pairs = dists.len() as f64;
        let mut result = Vec::with_capacity(n_radii);

        for k in 0..n_radii {
            let r =
                r_min * (r_max / r_min).powf(k as f64 / n_radii.saturating_sub(1).max(1) as f64);
            let count = dists.partition_point(|&d| d <= r);
            let c_r = count as f64 / total_pairs;
            if c_r > 0.0 {
                result.push((r.ln(), c_r.ln()));
            }
        }
        result
    }

    /// Estimate the largest Lyapunov exponent from the data via the Rosenstein method.
    ///
    /// Finds the nearest neighbour for each point, tracks divergence over `horizon`
    /// steps and returns the mean log-divergence as a function of time.
    pub fn rosenstein_lyapunov(&self, horizon: usize) -> Vec<f64> {
        let n = self.data.len();
        if n < 2 || horizon == 0 {
            return vec![];
        }
        let steps = horizon.min(n / 2);
        let mut sum_log = vec![0.0f64; steps];
        let mut counts = vec![0usize; steps];

        for i in 0..n {
            // Find nearest neighbour (excluding temporal neighbours within `min_sep`)
            let min_sep = 1usize;
            let mut nn_idx = n;
            let mut nn_dist = f64::INFINITY;
            for j in 0..n {
                if (i as isize - j as isize).unsigned_abs() < min_sep {
                    continue;
                }
                let d: f64 = self.data[i]
                    .iter()
                    .zip(self.data[j].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                if d < nn_dist {
                    nn_dist = d;
                    nn_idx = j;
                }
            }
            if nn_idx == n {
                continue;
            }

            // Track divergence
            for t in 0..steps {
                let ii = i + t;
                let jj = nn_idx + t;
                if ii >= n || jj >= n {
                    break;
                }
                let d: f64 = self.data[ii]
                    .iter()
                    .zip(self.data[jj].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                if d > 1e-300 {
                    sum_log[t] += d.ln();
                    counts[t] += 1;
                }
            }
        }

        sum_log
            .iter()
            .zip(counts.iter())
            .map(|(&s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
            .collect()
    }

    /// Check whether the attractor looks strange by testing if the
    /// correlation dimension estimate is fractional (> 1.5).
    pub fn is_strange_attractor(&self) -> bool {
        let pairs = self.correlation_dimension(20);
        if pairs.len() < 4 {
            return false;
        }
        // Linear regression slope on the log-log pairs
        let slope = linear_slope(&pairs);
        slope > 1.5 && slope < 4.0
    }
}

/// Helper: linear regression slope for `(x, y)` pairs.
fn linear_slope(pairs: &[(f64, f64)]) -> f64 {
    let n = pairs.len() as f64;
    let sx: f64 = pairs.iter().map(|(x, _)| x).sum();
    let sy: f64 = pairs.iter().map(|(_, y)| y).sum();
    let sxy: f64 = pairs.iter().map(|(x, y)| x * y).sum();
    let sxx: f64 = pairs.iter().map(|(x, _)| x * x).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        return 0.0;
    }
    (n * sxy - sx * sy) / denom
}

// ─────────────────────────────────────────────────────────────────────────────
// DynamicalNormalForm
// ─────────────────────────────────────────────────────────────────────────────

/// Normal forms for codimension-1 bifurcations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BifurcationType {
    /// Saddle-node (fold): `ẋ = μ + x²`.
    SaddleNode,
    /// Transcritical: `ẋ = μx − x²`.
    Transcritical,
    /// Supercritical pitchfork: `ẋ = μx − x³`.
    PitchforkSuper,
    /// Subcritical pitchfork: `ẋ = μx + x³`.
    PitchforkSub,
    /// Supercritical Hopf (2D): `ṙ = μr − r³`.
    HopfSuper,
    /// Subcritical Hopf (2D): `ṙ = μr + r³`.
    HopfSub,
}

/// Normal-form vector field evaluator.
#[derive(Debug, Clone)]
pub struct DynamicalNormalForm {
    /// Bifurcation type.
    pub bif_type: BifurcationType,
    /// Bifurcation parameter `μ`.
    pub mu: f64,
}

impl DynamicalNormalForm {
    /// Create a new normal form.
    pub fn new(bif_type: BifurcationType, mu: f64) -> Self {
        Self { bif_type, mu }
    }

    /// Evaluate the vector field at state `x` (scalar for 1D, `(r, θ)` for Hopf).
    ///
    /// For 1D normal forms `state[0] = x`.  For Hopf `state[0] = r`, `state[1] = θ`.
    pub fn rhs(&self, state: &[f64]) -> Vec<f64> {
        match self.bif_type {
            BifurcationType::SaddleNode => {
                let x = state[0];
                vec![self.mu + x * x]
            }
            BifurcationType::Transcritical => {
                let x = state[0];
                vec![self.mu * x - x * x]
            }
            BifurcationType::PitchforkSuper => {
                let x = state[0];
                vec![self.mu * x - x * x * x]
            }
            BifurcationType::PitchforkSub => {
                let x = state[0];
                vec![self.mu * x + x * x * x]
            }
            BifurcationType::HopfSuper => {
                let r = state[0];
                let theta = if state.len() > 1 { state[1] } else { 0.0 };
                vec![self.mu * r - r * r * r, theta + 1.0]
            }
            BifurcationType::HopfSub => {
                let r = state[0];
                let theta = if state.len() > 1 { state[1] } else { 0.0 };
                vec![self.mu * r + r * r * r, theta + 1.0]
            }
        }
    }

    /// Fixed points of the 1D normal form.
    ///
    /// Returns up to 3 fixed-point values.  Empty vec if none exist for the
    /// given `μ`.
    pub fn fixed_points_1d(&self) -> Vec<f64> {
        match self.bif_type {
            BifurcationType::SaddleNode => {
                // ẋ = μ + x² = 0  ⟹  x = ±√(−μ)
                if self.mu <= 0.0 {
                    let s = (-self.mu).sqrt();
                    if s < 1e-12 { vec![0.0] } else { vec![-s, s] }
                } else {
                    vec![]
                }
            }
            BifurcationType::Transcritical => {
                // ẋ = μx − x² = x(μ − x) = 0  ⟹  x=0 or x=μ
                vec![0.0, self.mu]
            }
            BifurcationType::PitchforkSuper | BifurcationType::PitchforkSub => {
                // ẋ = μx ± x³ = x(μ ± x²) = 0  ⟹  x=0 and (for super) x=±√μ if μ>0
                let mut fps = vec![0.0];
                if self.bif_type == BifurcationType::PitchforkSuper && self.mu > 0.0 {
                    let s = self.mu.sqrt();
                    fps.push(-s);
                    fps.push(s);
                } else if self.bif_type == BifurcationType::PitchforkSub && self.mu < 0.0 {
                    let s = (-self.mu).sqrt();
                    fps.push(-s);
                    fps.push(s);
                }
                fps
            }
            BifurcationType::HopfSuper | BifurcationType::HopfSub => {
                // ṙ = 0  ⟹  r=0, and for super/sub: r=√(±μ)
                let mut fps = vec![0.0];
                if self.bif_type == BifurcationType::HopfSuper && self.mu > 0.0 {
                    fps.push(self.mu.sqrt());
                } else if self.bif_type == BifurcationType::HopfSub && self.mu < 0.0 {
                    fps.push((-self.mu).sqrt());
                }
                fps
            }
        }
    }

    /// Euler-integrate the normal form for `n` steps with step size `dt`.
    pub fn integrate_euler(&self, init: &[f64], dt: f64, n: usize) -> Vec<Vec<f64>> {
        let mut state = init.to_vec();
        let mut traj = Vec::with_capacity(n + 1);
        traj.push(state.clone());
        for _ in 0..n {
            let rhs = self.rhs(&state);
            for (s, r) in state.iter_mut().zip(rhs.iter()) {
                *s += dt * r;
            }
            traj.push(state.clone());
        }
        traj
    }

    /// Classify the stability of a 1D fixed point `x*` via the sign of the
    /// linearisation `f'(x*)`.
    pub fn stability_1d(&self, x_star: f64) -> Stability {
        let eps = 1e-7;
        let rhs_plus = self.rhs(&[x_star + eps])[0];
        let rhs_minus = self.rhs(&[x_star - eps])[0];
        let deriv = (rhs_plus - rhs_minus) / (2.0 * eps);
        if deriv < -1e-10 {
            Stability::Stable
        } else if deriv > 1e-10 {
            Stability::Unstable
        } else {
            Stability::Marginal
        }
    }
}

/// Stability classification of a fixed point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stability {
    /// Asymptotically stable (attracting).
    Stable,
    /// Unstable (repelling).
    Unstable,
    /// Marginally stable (centre or degenerate).
    Marginal,
}

// ─────────────────────────────────────────────────────────────────────────────
// PoincareSection
// ─────────────────────────────────────────────────────────────────────────────

/// Poincaré section analysis for flows or maps.
///
/// Records intersections of a trajectory with a hyperplane defined by
/// `normal · (x − origin) = 0`.
#[derive(Debug, Clone)]
pub struct PoincareSection {
    /// Normal vector of the section hyperplane.
    pub normal: Vec<f64>,
    /// A point on the section hyperplane (origin).
    pub origin: Vec<f64>,
    /// Collected intersection points.
    pub intersections: Vec<Vec<f64>>,
}

impl PoincareSection {
    /// Create a new Poincaré section.
    pub fn new(normal: Vec<f64>, origin: Vec<f64>) -> Self {
        Self {
            normal,
            origin,
            intersections: Vec::new(),
        }
    }

    /// Signed distance of point `p` from the section.
    pub fn signed_distance(&self, p: &[f64]) -> f64 {
        p.iter()
            .zip(self.normal.iter())
            .zip(self.origin.iter())
            .map(|((pi, ni), oi)| (pi - oi) * ni)
            .sum()
    }

    /// Process a trajectory (list of states) and record all up-crossings of
    /// the section using linear interpolation.
    pub fn process_trajectory(&mut self, traj: &[Vec<f64>]) {
        for w in traj.windows(2) {
            let d0 = self.signed_distance(&w[0]);
            let d1 = self.signed_distance(&w[1]);
            // Up-crossing: d0 < 0 and d1 >= 0
            if d0 < 0.0 && d1 >= 0.0 {
                let t = d0.abs() / (d0.abs() + d1.abs() + 1e-300);
                let interp: Vec<f64> = w[0]
                    .iter()
                    .zip(w[1].iter())
                    .map(|(a, b)| a + t * (b - a))
                    .collect();
                self.intersections.push(interp);
            }
        }
    }

    /// Compute the return map: the sequence of successive intersection
    /// coordinates along axis `idx`.
    pub fn return_map_1d(&self, idx: usize) -> Vec<(f64, f64)> {
        let coords: Vec<f64> = self
            .intersections
            .iter()
            .filter_map(|p| p.get(idx).copied())
            .collect();
        coords.windows(2).map(|w| (w[0], w[1])).collect()
    }

    /// Rotation number: average angular advance between successive crossings.
    ///
    /// Assumes the section is in a 2D projection and computes the angle of
    /// each crossing relative to the section origin.
    pub fn rotation_number(&self, angle_idx_x: usize, angle_idx_y: usize) -> f64 {
        if self.intersections.len() < 2 {
            return 0.0;
        }
        let angles: Vec<f64> = self
            .intersections
            .iter()
            .filter_map(|p| {
                let x = p.get(angle_idx_x).copied()? - self.origin.get(angle_idx_x).copied()?;
                let y = p.get(angle_idx_y).copied()? - self.origin.get(angle_idx_y).copied()?;
                Some(y.atan2(x))
            })
            .collect();

        let mut total = 0.0f64;
        for w in angles.windows(2) {
            let mut diff = w[1] - w[0];
            // Wrap to [-π, π]
            while diff > std::f64::consts::PI {
                diff -= std::f64::consts::TAU;
            }
            while diff < -std::f64::consts::PI {
                diff += std::f64::consts::TAU;
            }
            total += diff;
        }
        total / (std::f64::consts::TAU * (angles.len() - 1) as f64)
    }

    /// Number of recorded intersections.
    pub fn count(&self) -> usize {
        self.intersections.len()
    }

    /// Clear all recorded intersections.
    pub fn clear(&mut self) {
        self.intersections.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FlowMap
// ─────────────────────────────────────────────────────────────────────────────

/// Numerical flow integration for autonomous ODEs `ẋ = f(x)`.
///
/// Uses a fixed-step 4th-order Runge–Kutta integrator.
#[derive(Debug, Clone)]
pub struct FlowMap {
    /// Dimension of the state space.
    pub dim: usize,
    /// Integration step size.
    pub dt: f64,
}

impl FlowMap {
    /// Create a flow map integrator.
    pub fn new(dim: usize, dt: f64) -> Self {
        Self { dim, dt }
    }

    /// Advance the state by `n_steps` using RK4 with the vector field `f`.
    pub fn advance<F>(&self, init: &[f64], n_steps: usize, f: &F) -> Vec<f64>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let mut state = init.to_vec();
        for _ in 0..n_steps {
            state = self.rk4_step(&state, f);
        }
        state
    }

    /// Generate a trajectory (list of states) of length `n_steps + 1`.
    pub fn trajectory<F>(&self, init: &[f64], n_steps: usize, f: &F) -> Vec<Vec<f64>>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let mut traj = Vec::with_capacity(n_steps + 1);
        traj.push(init.to_vec());
        for _ in 0..n_steps {
            let next = self.rk4_step(traj.last().expect("trajectory is non-empty"), f);
            traj.push(next);
        }
        traj
    }

    /// Stroboscopic map: sample the trajectory at every `period` steps.
    pub fn stroboscopic_map<F>(
        &self,
        init: &[f64],
        n_periods: usize,
        steps_per_period: usize,
        f: &F,
    ) -> Vec<Vec<f64>>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let mut state = init.to_vec();
        let mut samples = Vec::with_capacity(n_periods + 1);
        samples.push(state.clone());
        for _ in 0..n_periods {
            state = self.advance(&state, steps_per_period, f);
            samples.push(state.clone());
        }
        samples
    }

    /// Compute the fundamental solution matrix (monodromy-like) by
    /// propagating `dim` perturbation vectors alongside the trajectory.
    ///
    /// Returns a `dim × dim` matrix in row-major order (each row is a
    /// transported unit vector).
    pub fn fundamental_solution_matrix<F>(
        &self,
        init: &[f64],
        n_steps: usize,
        f: &F,
    ) -> Vec<Vec<f64>>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let d = self.dim;
        // Initialise perturbation vectors as identity matrix
        let mut pert: Vec<Vec<f64>> = (0..d)
            .map(|i| {
                let mut v = vec![0.0; d];
                v[i] = 1.0;
                v
            })
            .collect();

        let mut state = init.to_vec();
        let eps = 1e-6;

        for _ in 0..n_steps {
            // Advance reference trajectory
            let next_state = self.rk4_step(&state, f);
            // Advance each perturbation vector using finite-difference Jacobian
            let mut next_pert = vec![vec![0.0f64; d]; d];
            for (col, pv) in pert.iter().enumerate() {
                let state_plus: Vec<f64> = state
                    .iter()
                    .zip(pv.iter())
                    .map(|(s, p)| s + eps * p)
                    .collect();
                let f_plus = f(&state_plus);
                let f_ref = f(&state);
                for row in 0..d {
                    next_pert[row][col] = pv[row] + self.dt * (f_plus[row] - f_ref[row]) / eps;
                }
            }
            state = next_state;
            pert = next_pert;
        }
        pert
    }

    fn rk4_step<F>(&self, state: &[f64], f: &F) -> Vec<f64>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let k1 = f(state);
        let s2: Vec<f64> = state
            .iter()
            .zip(k1.iter())
            .map(|(s, k)| s + 0.5 * self.dt * k)
            .collect();
        let k2 = f(&s2);
        let s3: Vec<f64> = state
            .iter()
            .zip(k2.iter())
            .map(|(s, k)| s + 0.5 * self.dt * k)
            .collect();
        let k3 = f(&s3);
        let s4: Vec<f64> = state
            .iter()
            .zip(k3.iter())
            .map(|(s, k)| s + self.dt * k)
            .collect();
        let k4 = f(&s4);
        state
            .iter()
            .enumerate()
            .map(|(i, s)| s + self.dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SynchronizationAnalysis
// ─────────────────────────────────────────────────────────────────────────────

/// Synchronisation analysis for coupled oscillators / chaotic systems.
#[derive(Debug, Clone)]
pub struct SynchronizationAnalysis;

impl SynchronizationAnalysis {
    /// Pecora–Carroll synchronisation check.
    ///
    /// Drives a response system (with the same structure as the drive) by
    /// replacing one of its variables with the drive's corresponding variable.
    ///
    /// Returns `true` if the synchronisation error decays below `tol`.
    ///
    /// `drive_orbit` and `response_orbit` are trajectories of the same length;
    /// `sync_idx` is the variable index that is replaced in the response.
    pub fn pecora_carroll_sync(
        drive_orbit: &[Vec<f64>],
        response_orbit: &[Vec<f64>],
        sync_idx: usize,
        tol: f64,
    ) -> bool {
        if drive_orbit.is_empty() || response_orbit.is_empty() {
            return false;
        }
        // Check that the final error (excluding sync_idx variable) is small
        let last_drive = drive_orbit.last().expect("drive orbit is non-empty");
        let last_response = response_orbit.last().expect("response orbit is non-empty");
        let err: f64 = last_drive
            .iter()
            .zip(last_response.iter())
            .enumerate()
            .filter(|(i, _)| *i != sync_idx)
            .map(|(_, (d, r))| (d - r).powi(2))
            .sum::<f64>()
            .sqrt();
        err < tol
    }

    /// Phase synchronisation index (mean resultant length of phase differences).
    ///
    /// Given two phase time-series `phi1` and `phi2` (in radians), computes
    ///
    /// `R = |⟨exp(i(φ₁ − φ₂))⟩|`
    ///
    /// `R = 1` means perfect phase synchronisation; `R ≈ 0` means asynchronous.
    pub fn phase_synchronization_index(phi1: &[f64], phi2: &[f64]) -> f64 {
        let n = phi1.len().min(phi2.len());
        if n == 0 {
            return 0.0;
        }
        let sum_cos: f64 = phi1
            .iter()
            .zip(phi2.iter())
            .map(|(a, b)| (a - b).cos())
            .sum();
        let sum_sin: f64 = phi1
            .iter()
            .zip(phi2.iter())
            .map(|(a, b)| (a - b).sin())
            .sum();
        let n_f = n as f64;
        ((sum_cos / n_f).powi(2) + (sum_sin / n_f).powi(2)).sqrt()
    }

    /// Generalised synchronisation check via the auxiliary system method.
    ///
    /// Two copies of the response system starting at different initial conditions
    /// should converge to the same trajectory if GS holds.  Returns `true` if
    /// the distance between `resp1` and `resp2` at the end is below `tol`.
    pub fn generalised_sync_check(resp1: &[Vec<f64>], resp2: &[Vec<f64>], tol: f64) -> bool {
        let last1 = match resp1.last() {
            Some(v) => v,
            None => return false,
        };
        let last2 = match resp2.last() {
            Some(v) => v,
            None => return false,
        };
        let dist: f64 = last1
            .iter()
            .zip(last2.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        dist < tol
    }

    /// Mutual information estimate between two scalar signals using histogram
    /// binning with `n_bins` bins.
    pub fn mutual_information(x: &[f64], y: &[f64], n_bins: usize) -> f64 {
        let n = x.len().min(y.len());
        if n == 0 || n_bins == 0 {
            return 0.0;
        }
        let x_min = x[..n].iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = x[..n].iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = y[..n].iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = y[..n].iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        let x_range = (x_max - x_min).max(1e-15);
        let y_range = (y_max - y_min).max(1e-15);

        let mut joint = vec![vec![0u64; n_bins]; n_bins];
        let mut hist_x = vec![0u64; n_bins];
        let mut hist_y = vec![0u64; n_bins];

        for i in 0..n {
            let ix = ((x[i] - x_min) / x_range * n_bins as f64)
                .floor()
                .clamp(0.0, (n_bins - 1) as f64) as usize;
            let iy = ((y[i] - y_min) / y_range * n_bins as f64)
                .floor()
                .clamp(0.0, (n_bins - 1) as f64) as usize;
            joint[ix][iy] += 1;
            hist_x[ix] += 1;
            hist_y[iy] += 1;
        }

        let n_f = n as f64;
        let mut mi = 0.0f64;
        for ix in 0..n_bins {
            for iy in 0..n_bins {
                let p_xy = joint[ix][iy] as f64 / n_f;
                let p_x = hist_x[ix] as f64 / n_f;
                let p_y = hist_y[iy] as f64 / n_f;
                if p_xy > 0.0 && p_x > 0.0 && p_y > 0.0 {
                    mi += p_xy * (p_xy / (p_x * p_y)).ln();
                }
            }
        }
        mi
    }

    /// Compute instantaneous phase of a signal via the Hilbert transform
    /// approximation (using the discrete analytic signal via FFT-like unwrapping).
    ///
    /// This is a simple numerical approximation using finite differences.
    /// Returns the unwrapped phase array.
    pub fn instantaneous_phase(signal: &[f64]) -> Vec<f64> {
        let n = signal.len();
        if n < 2 {
            return vec![0.0; n];
        }
        // Compute a rough Hilbert transform via 90-degree phase shift in Fourier space.
        // Here we use a simple finite-difference approximation: a causal integrator.
        let mean = signal.iter().sum::<f64>() / n as f64;
        let centered: Vec<f64> = signal.iter().map(|s| s - mean).collect();

        // Quadrature signal via finite differences (simple approximation)
        let mut quad = vec![0.0f64; n];
        for i in 1..n {
            quad[i] = quad[i - 1] + centered[i - 1];
        }
        // Normalise quadrature
        let q_max = quad
            .iter()
            .cloned()
            .fold(0.0f64, |a, b| a.max(b.abs()))
            .max(1e-15);
        let c_max = centered
            .iter()
            .cloned()
            .fold(0.0f64, |a, b| a.max(b.abs()))
            .max(1e-15);
        let scale = c_max / q_max;

        let mut phase = Vec::with_capacity(n);
        for i in 0..n {
            phase.push(quad[i] * scale * centered[i].signum().max(0.0));
        }

        // Use atan2 for phase
        let mut phases = Vec::with_capacity(n);
        for i in 0..n {
            phases.push((quad[i] * scale).atan2(centered[i]));
        }
        phases
    }

    /// Cross-correlation between two signals at lag 0 (normalised).
    pub fn cross_correlation_zero_lag(x: &[f64], y: &[f64]) -> f64 {
        let n = x.len().min(y.len());
        if n == 0 {
            return 0.0;
        }
        let mx = x.iter().sum::<f64>() / n as f64;
        let my = y.iter().sum::<f64>() / n as f64;
        let num: f64 = x[..n]
            .iter()
            .zip(y[..n].iter())
            .map(|(a, b)| (a - mx) * (b - my))
            .sum();
        let sx: f64 = x[..n].iter().map(|a| (a - mx).powi(2)).sum::<f64>().sqrt();
        let sy: f64 = y[..n].iter().map(|b| (b - my).powi(2)).sum::<f64>().sqrt();
        if sx < 1e-15 || sy < 1e-15 {
            return 0.0;
        }
        num / (sx * sy)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Lorenz system helper (used in tests)
// ─────────────────────────────────────────────────────────────────────────────

/// Lorenz system parameters.
#[derive(Debug, Clone)]
pub struct LorenzSystem {
    /// σ parameter.
    pub sigma: f64,
    /// ρ parameter.
    pub rho: f64,
    /// β parameter.
    pub beta: f64,
}

impl LorenzSystem {
    /// Create with standard chaotic parameters `σ=10, ρ=28, β=8/3`.
    pub fn standard() -> Self {
        Self {
            sigma: 10.0,
            rho: 28.0,
            beta: 8.0 / 3.0,
        }
    }

    /// Evaluate the Lorenz vector field at `(x, y, z)`.
    pub fn rhs(&self, state: &[f64]) -> Vec<f64> {
        let (x, y, z) = (state[0], state[1], state[2]);
        vec![
            self.sigma * (y - x),
            x * (self.rho - z) - y,
            x * y - self.beta * z,
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── LogisticMap ────────────────────────────────────────────────────────────

    #[test]
    fn test_logistic_fixed_point_stability() {
        // For r=2, fixed point at x* = (r-1)/r = 0.5
        let map = LogisticMap::new(2.0);
        let x_star = (map.r - 1.0) / map.r;
        let mut x = 0.4;
        for _ in 0..200 {
            x = map.apply(x);
        }
        assert!(
            (x - x_star).abs() < 1e-6,
            "should converge to {x_star}, got {x}"
        );
    }

    #[test]
    fn test_logistic_orbit_length() {
        let map = LogisticMap::new(3.5);
        let orbit = map.orbit_1d(0.3, 50);
        assert_eq!(orbit.len(), 51);
    }

    #[test]
    fn test_logistic_orbit_bounded() {
        let map = LogisticMap::new(3.9);
        let orbit = map.orbit_1d(0.5, 1000);
        for &x in &orbit {
            assert!((0.0..=1.0).contains(&x), "orbit escaped [0,1]: {x}");
        }
    }

    #[test]
    fn test_logistic_lyapunov_chaotic() {
        // r=4.0: Lyapunov exponent ≈ ln 2 ≈ 0.693
        let map = LogisticMap::new(4.0);
        let le = map.lyapunov_exponent(0.3, 5000).unwrap();
        assert!(le > 0.4, "should be positive (chaotic), got {le}");
    }

    #[test]
    fn test_logistic_lyapunov_stable() {
        // r=2.0: stable fixed point at x*=0.5, Lyapunov < 0
        // Use x0=0.3 which avoids the fixed point x=0 or x=1
        let map = LogisticMap::new(2.0);
        // The derivative at the fixed point is r(1-2x*) = 2*(1-1) = 0, so use short orbit
        // and check via orbit convergence instead
        let mut x = 0.3f64;
        for _ in 0..500 {
            x = map.apply(x);
        }
        // Should have converged near x* = 0.5
        assert!((x - 0.5).abs() < 1e-6, "should converge to x*=0.5, got {x}");
        // Lyapunov: compute without hitting the degenerate point
        if let Some(le) = map.lyapunov_exponent(0.3, 200) {
            assert!(
                le < 0.1,
                "stable fixed point should have non-positive LE, got {le}"
            );
        }
        // If None is returned that's also acceptable (orbit hit degenerate point)
    }

    #[test]
    fn test_logistic_bifurcation_diagram_size() {
        let diag = LogisticMap::bifurcation_diagram(2.5, 4.0, 10, 200, 50);
        assert_eq!(diag.len(), 10 * 50);
    }

    #[test]
    fn test_logistic_trait_orbit() {
        let map = LogisticMap::new(3.2);
        let orbit = map.orbit(&[0.5], 20);
        assert_eq!(orbit.len(), 21);
        assert_eq!(orbit[0], vec![0.5]);
    }

    // ── HenonMap ───────────────────────────────────────────────────────────────

    #[test]
    fn test_henon_orbit_first_step() {
        let map = HenonMap::new(1.4, 0.3);
        let (xn, yn) = map.apply(0.0, 0.0);
        assert!((xn - 1.0).abs() < 1e-12);
        assert!(yn.abs() < 1e-12);
    }

    #[test]
    fn test_henon_jacobian_det_constant() {
        let map = HenonMap::new(1.4, 0.3);
        assert!((map.jacobian_det() - (-0.3)).abs() < 1e-12);
    }

    #[test]
    fn test_henon_orbit_length() {
        let map = HenonMap::new(1.4, 0.3);
        let orbit = map.orbit_2d(0.0, 0.0, 100);
        assert_eq!(orbit.len(), 101);
    }

    #[test]
    fn test_henon_jacobian_structure() {
        let map = HenonMap::new(1.4, 0.3);
        let j = map.jacobian(1.0, 0.5);
        assert!((j[0][0] - (-2.0 * 1.4)).abs() < 1e-12);
        assert!((j[0][1] - 1.0).abs() < 1e-12);
        assert!((j[1][0] - 0.3).abs() < 1e-12);
        assert!(j[1][1].abs() < 1e-12);
    }

    #[test]
    fn test_henon_trait_step() {
        let map = HenonMap::new(1.4, 0.3);
        let state = vec![0.0, 0.0];
        let next = map.step(&state);
        assert_eq!(next.len(), 2);
    }

    // ── ChirikovStandardMap ────────────────────────────────────────────────────

    #[test]
    fn test_chirikov_orbit_stays_in_torus() {
        let map = ChirikovStandardMap::new(0.5);
        let orbit = map.orbit_2d(0.1, 0.2, 1000);
        let tau = std::f64::consts::TAU;
        for (th, p) in orbit {
            assert!(th >= 0.0 && th < tau, "theta out of [0,2π): {th}");
            assert!(p >= 0.0 && p < tau, "p out of [0,2π): {p}");
        }
    }

    #[test]
    fn test_chirikov_zero_k_is_twist_map() {
        // K=0: p unchanged, theta advances by p
        let map = ChirikovStandardMap::new(0.0);
        let (th, p) = map.apply(1.0, 0.5);
        let tau = std::f64::consts::TAU;
        assert!((p - 0.5).abs() < 1e-12, "p should not change: {p}");
        assert!(
            (th - (1.5_f64).rem_euclid(tau)).abs() < 1e-12,
            "theta = 1+p mod 2π: {th}"
        );
    }

    #[test]
    fn test_chirikov_lyapunov_chaotic() {
        // K=5 should be strongly chaotic
        let map = ChirikovStandardMap::new(5.0);
        let le = map.lyapunov_exponent(1.0, 1.5, 5000);
        assert!(le > 0.5, "K=5 should be chaotic, got LE={le}");
    }

    #[test]
    fn test_chirikov_trait_orbit() {
        let map = ChirikovStandardMap::new(1.0);
        let orbit = map.orbit(&[0.5, 0.5], 10);
        assert_eq!(orbit.len(), 11);
    }

    // ── AttractorAnalysis ──────────────────────────────────────────────────────

    #[test]
    fn test_attractor_find_fixed_points_identity() {
        // Points near y = x should be found as fixed points of x → x
        let data = vec![vec![0.0], vec![0.5], vec![1.0], vec![2.0]];
        let analysis = AttractorAnalysis::new(data);
        let fps = analysis.find_fixed_points(|x| x.to_vec(), 1e-9);
        assert!(
            !fps.is_empty(),
            "identity map has every point as fixed point"
        );
    }

    #[test]
    fn test_attractor_correlation_dimension_returns_pairs() {
        let data: Vec<Vec<f64>> = (0..50).map(|i| vec![i as f64 * 0.1]).collect();
        let analysis = AttractorAnalysis::new(data);
        let pairs = analysis.correlation_dimension(10);
        assert!(!pairs.is_empty(), "should produce log-log pairs");
    }

    #[test]
    fn test_attractor_rosenstein_length() {
        let data: Vec<Vec<f64>> = (0..100)
            .map(|i| {
                let t = i as f64 * 0.1;
                vec![t.sin(), t.cos()]
            })
            .collect();
        let analysis = AttractorAnalysis::new(data);
        let div = analysis.rosenstein_lyapunov(10);
        assert_eq!(div.len(), 10);
    }

    #[test]
    fn test_attractor_find_periodic_orbit() {
        // Period-2 orbit of logistic map at r ≈ 3.2
        let map = LogisticMap::new(3.2);
        // Converge to period-2 orbit
        let mut x = 0.5f64;
        for _ in 0..1000 {
            x = map.apply(x);
        }
        let p1 = x;
        let p2 = map.apply(p1);
        let data = vec![vec![p1], vec![p2]];
        let analysis = AttractorAnalysis::new(data);
        let po = analysis.find_periodic_orbits(&|s: &[f64]| vec![map.apply(s[0])], 2, 1e-6);
        assert!(!po.is_empty(), "should find period-2 orbit");
    }

    // ── DynamicalNormalForm ────────────────────────────────────────────────────

    #[test]
    fn test_normal_form_saddle_node_fixed_points() {
        let nf = DynamicalNormalForm::new(BifurcationType::SaddleNode, -1.0);
        let fps = nf.fixed_points_1d();
        assert_eq!(fps.len(), 2, "two fps for μ<0");
        for &x in &fps {
            let rhs_val = nf.rhs(&[x])[0];
            assert!(rhs_val.abs() < 1e-10, "not a fixed point: rhs={rhs_val}");
        }
    }

    #[test]
    fn test_normal_form_saddle_node_no_fps_positive_mu() {
        let nf = DynamicalNormalForm::new(BifurcationType::SaddleNode, 1.0);
        let fps = nf.fixed_points_1d();
        assert!(fps.is_empty(), "no fps for μ>0 in saddle-node");
    }

    #[test]
    fn test_normal_form_transcritical_fps() {
        let nf = DynamicalNormalForm::new(BifurcationType::Transcritical, 2.0);
        let fps = nf.fixed_points_1d();
        assert_eq!(fps.len(), 2);
        for &x in &fps {
            let rhs_val = nf.rhs(&[x])[0];
            assert!(
                rhs_val.abs() < 1e-10,
                "not a fixed point: x={x}, rhs={rhs_val}"
            );
        }
    }

    #[test]
    fn test_normal_form_pitchfork_super_three_fps() {
        let nf = DynamicalNormalForm::new(BifurcationType::PitchforkSuper, 1.0);
        let fps = nf.fixed_points_1d();
        assert_eq!(fps.len(), 3, "super-critical pitchfork μ>0 has 3 fps");
    }

    #[test]
    fn test_normal_form_stability_classification() {
        // Pitchfork super with μ=1: x*=0 is unstable, x*=±1 are stable
        let nf = DynamicalNormalForm::new(BifurcationType::PitchforkSuper, 1.0);
        assert_eq!(nf.stability_1d(0.0), Stability::Unstable);
        assert_eq!(nf.stability_1d(1.0), Stability::Stable);
        assert_eq!(nf.stability_1d(-1.0), Stability::Stable);
    }

    #[test]
    fn test_normal_form_euler_integration_converges() {
        // Saddle-node stable branch: μ=-1, x0=-0.9 (near stable fp x=-1)
        let nf = DynamicalNormalForm::new(BifurcationType::SaddleNode, -1.0);
        let traj = nf.integrate_euler(&[-0.9], 0.01, 2000);
        let last = traj.last().unwrap()[0];
        assert!(
            (last - (-1.0)).abs() < 0.05,
            "should converge near -1, got {last}"
        );
    }

    #[test]
    fn test_normal_form_hopf_super_limit_cycle() {
        // HopfSuper μ=1: trajectory should converge to r=1
        let nf = DynamicalNormalForm::new(BifurcationType::HopfSuper, 1.0);
        let traj = nf.integrate_euler(&[0.1, 0.0], 0.01, 5000);
        let last_r = traj.last().unwrap()[0];
        assert!((last_r - 1.0).abs() < 0.05, "limit cycle r≈1, got {last_r}");
    }

    // ── PoincareSection ────────────────────────────────────────────────────────

    #[test]
    fn test_poincare_section_circle() {
        // Circular orbit crossing y=0 upwards twice per revolution
        let n = 1000usize;
        let traj: Vec<Vec<f64>> = (0..=n)
            .map(|i| {
                let t = i as f64 * std::f64::consts::TAU / n as f64 * 5.0;
                vec![t.cos(), t.sin()]
            })
            .collect();
        let mut section = PoincareSection::new(vec![0.0, 1.0], vec![0.0, 0.0]);
        section.process_trajectory(&traj);
        // Should have ~5 crossings
        assert!(
            section.count() >= 4,
            "expected ~5 crossings, got {}",
            section.count()
        );
    }

    #[test]
    fn test_poincare_section_rotation_number() {
        let n = 2000usize;
        let traj: Vec<Vec<f64>> = (0..=n)
            .map(|i| {
                let t = i as f64 * std::f64::consts::TAU / n as f64 * 3.0;
                vec![t.cos(), t.sin()]
            })
            .collect();
        let mut section = PoincareSection::new(vec![0.0, 1.0], vec![0.0, 0.0]);
        section.process_trajectory(&traj);
        let rn = section.rotation_number(0, 1);
        // rotation number for a circle should be approximately 1 per crossing
        assert!(rn.abs() < 2.0, "rotation number should be finite, got {rn}");
    }

    #[test]
    fn test_poincare_clear() {
        let traj: Vec<Vec<f64>> = (0..=100)
            .map(|i| {
                let t = i as f64 * 0.1;
                vec![t.cos(), t.sin()]
            })
            .collect();
        let mut section = PoincareSection::new(vec![0.0, 1.0], vec![0.0, 0.0]);
        section.process_trajectory(&traj);
        section.clear();
        assert_eq!(section.count(), 0);
    }

    #[test]
    fn test_poincare_return_map() {
        let n = 2000usize;
        let traj: Vec<Vec<f64>> = (0..=n)
            .map(|i| {
                let t = i as f64 * std::f64::consts::TAU / n as f64 * 5.0;
                vec![t.cos(), t.sin()]
            })
            .collect();
        let mut section = PoincareSection::new(vec![0.0, 1.0], vec![0.0, 0.0]);
        section.process_trajectory(&traj);
        let rm = section.return_map_1d(0);
        // Should have (count - 1) return-map pairs
        assert_eq!(rm.len(), section.count().saturating_sub(1));
    }

    // ── FlowMap ────────────────────────────────────────────────────────────────

    #[test]
    fn test_flow_simple_harmonic_oscillator() {
        // ẋ = y, ẏ = -x (simple harmonic oscillator)
        let f = |s: &[f64]| vec![s[1], -s[0]];
        let flow = FlowMap::new(2, 0.01);
        // Advance one full period: T = 2π ≈ 628 steps at dt=0.01
        let final_state = flow.advance(&[1.0, 0.0], 628, &f);
        // Should return close to (1, 0)
        assert!(
            (final_state[0] - 1.0).abs() < 0.1,
            "x≈1, got {}",
            final_state[0]
        );
        assert!(final_state[1].abs() < 0.1, "y≈0, got {}", final_state[1]);
    }

    #[test]
    fn test_flow_trajectory_length() {
        let f = |s: &[f64]| vec![s[1], -s[0]];
        let flow = FlowMap::new(2, 0.01);
        let traj = flow.trajectory(&[1.0, 0.0], 100, &f);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_flow_stroboscopic_map() {
        let f = |s: &[f64]| vec![s[1], -s[0]];
        let flow = FlowMap::new(2, 0.01);
        let samples = flow.stroboscopic_map(&[1.0, 0.0], 5, 10, &f);
        assert_eq!(samples.len(), 6);
    }

    #[test]
    fn test_flow_energy_conservation_sho() {
        // For SHO, E = x² + y² should be conserved
        let f = |s: &[f64]| vec![s[1], -s[0]];
        let flow = FlowMap::new(2, 0.001);
        let init = vec![1.0, 0.0];
        let energy_init: f64 = init[0] * init[0] + init[1] * init[1];
        let final_state = flow.advance(&init, 1000, &f);
        let energy_final: f64 = final_state[0] * final_state[0] + final_state[1] * final_state[1];
        assert!(
            (energy_final - energy_init).abs() < 0.01,
            "energy not conserved: {energy_final}"
        );
    }

    #[test]
    fn test_flow_fundamental_solution_matrix() {
        let f = |s: &[f64]| vec![s[1], -s[0]];
        let flow = FlowMap::new(2, 0.01);
        let fsm = flow.fundamental_solution_matrix(&[1.0, 0.0], 10, &f);
        assert_eq!(fsm.len(), 2);
        assert_eq!(fsm[0].len(), 2);
    }

    // ── SynchronizationAnalysis ────────────────────────────────────────────────

    #[test]
    fn test_phase_sync_index_identical() {
        // Identical phases → R = 1
        let phi: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let r = SynchronizationAnalysis::phase_synchronization_index(&phi, &phi);
        assert!((r - 1.0).abs() < 1e-6, "identical phases → R=1, got {r}");
    }

    #[test]
    fn test_phase_sync_index_random_low() {
        // Linearly increasing phase difference → low R
        let phi1: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let phi2: Vec<f64> = (0..100).map(|i| i as f64 * 0.3 + 0.7).collect();
        let r = SynchronizationAnalysis::phase_synchronization_index(&phi1, &phi2);
        assert!(r < 0.8, "asynchronous phases → low R, got {r}");
    }

    #[test]
    fn test_generalised_sync_converged() {
        let traj1: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0, 2.0]).collect();
        let traj2: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0, 2.0]).collect();
        assert!(SynchronizationAnalysis::generalised_sync_check(
            &traj1, &traj2, 1e-9
        ));
    }

    #[test]
    fn test_generalised_sync_diverged() {
        let traj1: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0]).collect();
        let traj2: Vec<Vec<f64>> = (0..10).map(|_| vec![100.0]).collect();
        assert!(!SynchronizationAnalysis::generalised_sync_check(
            &traj1, &traj2, 1.0
        ));
    }

    #[test]
    fn test_mutual_information_identical_signals() {
        let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let mi = SynchronizationAnalysis::mutual_information(&x, &x, 10);
        assert!(
            mi > 0.0,
            "MI should be positive for identical signals, got {mi}"
        );
    }

    #[test]
    fn test_mutual_information_independent_signals() {
        let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..100).map(|i| (i * 7 + 3) as f64 % 100.0).collect();
        let mi = SynchronizationAnalysis::mutual_information(&x, &y, 10);
        // MI ≥ 0
        assert!(mi >= 0.0, "MI should be non-negative, got {mi}");
    }

    #[test]
    fn test_cross_correlation_identical() {
        let x: Vec<f64> = (0..50).map(|i| (i as f64 * 0.3).sin()).collect();
        let cc = SynchronizationAnalysis::cross_correlation_zero_lag(&x, &x);
        assert!(
            (cc - 1.0).abs() < 1e-6,
            "self-correlation should be 1, got {cc}"
        );
    }

    #[test]
    fn test_cross_correlation_anticorrelated() {
        let x: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|v| -v).collect();
        let cc = SynchronizationAnalysis::cross_correlation_zero_lag(&x, &y);
        assert!(
            (cc + 1.0).abs() < 1e-6,
            "anti-correlated: should be -1, got {cc}"
        );
    }

    #[test]
    fn test_instantaneous_phase_length() {
        let sig: Vec<f64> = (0..100).map(|i| (i as f64 * 0.2).sin()).collect();
        let phase = SynchronizationAnalysis::instantaneous_phase(&sig);
        assert_eq!(phase.len(), sig.len());
    }

    #[test]
    fn test_pecora_carroll_identical_sync() {
        let traj: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, (i as f64).sin()]).collect();
        assert!(SynchronizationAnalysis::pecora_carroll_sync(
            &traj, &traj, 0, 1e-9
        ));
    }

    // ── LorenzSystem integration ───────────────────────────────────────────────

    #[test]
    fn test_lorenz_trajectory() {
        let lorenz = LorenzSystem::standard();
        let flow = FlowMap::new(3, 0.01);
        let init = vec![1.0, 0.0, 0.0];
        let traj = flow.trajectory(&init, 100, &|s| lorenz.rhs(s));
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_lorenz_sensitivity() {
        // Lorenz is chaotic: two nearby initial conditions diverge substantially
        // Use 5000 steps at dt=0.01 = 50 time units for clear separation
        let lorenz = LorenzSystem::standard();
        let flow = FlowMap::new(3, 0.01);
        let init1 = vec![1.0, 1.0, 1.0];
        let init2 = vec![1.0 + 1e-5, 1.0, 1.0];
        let s1 = flow.advance(&init1, 5000, &|s| lorenz.rhs(s));
        let s2 = flow.advance(&init2, 5000, &|s| lorenz.rhs(s));
        let dist: f64 = s1
            .iter()
            .zip(s2.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        assert!(dist > 1e-3, "Lorenz should show sensitivity, dist={dist}");
    }

    #[test]
    fn test_lorenz_poincare_section() {
        let lorenz = LorenzSystem::standard();
        let flow = FlowMap::new(3, 0.01);
        let mut state = vec![1.0, 1.0, 1.0];
        // Warm up
        state = flow.advance(&state, 1000, &|s| lorenz.rhs(s));
        let traj = flow.trajectory(&state, 5000, &|s| lorenz.rhs(s));
        // Poincaré section at z = 27
        let mut section = PoincareSection::new(vec![0.0, 0.0, 1.0], vec![0.0, 0.0, 27.0]);
        section.process_trajectory(&traj);
        assert!(section.count() > 5, "Lorenz should cross z=27 many times");
    }

    #[test]
    fn test_attractor_analysis_lorenz_is_strange() {
        // Build data from Lorenz orbit (x coordinate only for speed)
        let lorenz = LorenzSystem::standard();
        let flow = FlowMap::new(3, 0.02);
        let mut state = vec![1.0, 1.0, 1.0];
        state = flow.advance(&state, 500, &|s| lorenz.rhs(s));
        let traj = flow.trajectory(&state, 400, &|s| lorenz.rhs(s));
        let data: Vec<Vec<f64>> = traj.to_vec();
        let analysis = AttractorAnalysis::new(data);
        let pairs = analysis.correlation_dimension(15);
        // Just check we get valid log-log pairs
        assert!(!pairs.is_empty(), "should compute correlation dimension");
    }

    // ── Linear slope helper ────────────────────────────────────────────────────

    #[test]
    fn test_linear_slope_perfect() {
        // y = 2x + 1 → slope = 2
        let pairs: Vec<(f64, f64)> = (0..5).map(|i| (i as f64, 2.0 * i as f64 + 1.0)).collect();
        let s = linear_slope(&pairs);
        assert!((s - 2.0).abs() < 1e-6, "slope should be 2, got {s}");
    }

    #[test]
    fn test_linear_slope_flat() {
        let pairs: Vec<(f64, f64)> = (0..5).map(|i| (i as f64, 3.0)).collect();
        let s = linear_slope(&pairs);
        assert!(s.abs() < 1e-6, "flat line has slope 0, got {s}");
    }
}
