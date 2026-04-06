// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Scientific plotting data structures for dynamical-systems analysis.
//!
//! All types are pure Rust, renderer-agnostic, and dependency-free (no nalgebra).
//! They represent *what* to draw; a downstream renderer converts them to pixels.
//!
//! # Structures
//! - [`PlotPoint`] — 2-D point `(x, y)`
//! - [`PlotSeries`] — named series of [`PlotPoint`]s with an RGBA colour
//! - [`PhasePortrait`] — 2-D vector-field grid, fixed-point detection, stability
//! - [`BifurcationDiagram`] — attractor vs. parameter sweep, period detection
//! - [`PoincareSection`] — trajectory → section crossings → 2-D return map
//! - [`LyapunovSpectrum`] — maximal Lyapunov exponent (Benettin algorithm)
//! - [`AttractorData`] — 3-D trajectory storage + correlation-dimension estimate
//! - [`TimeFrequencyData`] — short-time Fourier transform (STFT) spectrogram
//! - [`ScientificAxis`] — linear/log axis with Wilkinson-style tick generation

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// PlotPoint
// ---------------------------------------------------------------------------

/// A 2-D data point with `x` and `y` coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlotPoint {
    /// Horizontal coordinate.
    pub x: f64,
    /// Vertical coordinate.
    pub y: f64,
}

impl PlotPoint {
    /// Construct a new point.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Euclidean distance to another point.
    pub fn distance_to(&self, other: &PlotPoint) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        (dx * dx + dy * dy).sqrt()
    }
}

// ---------------------------------------------------------------------------
// PlotSeries
// ---------------------------------------------------------------------------

/// A named series of [`PlotPoint`]s with an associated RGBA colour.
#[derive(Debug, Clone)]
pub struct PlotSeries {
    /// Human-readable series name (used for legends).
    pub name: String,
    /// Ordered data points.
    pub points: Vec<PlotPoint>,
    /// RGBA colour — each component in `[0.0, 1.0]`.
    pub color: [f32; 4],
}

impl PlotSeries {
    /// Create a new, empty series.
    pub fn new(name: impl Into<String>, color: [f32; 4]) -> Self {
        Self {
            name: name.into(),
            points: Vec::new(),
            color,
        }
    }

    /// Append a point.
    pub fn push(&mut self, x: f64, y: f64) {
        self.points.push(PlotPoint::new(x, y));
    }

    /// Number of data points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Returns `true` if there are no data points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Compute the axis-aligned bounding box `(x_min, x_max, y_min, y_max)`.
    ///
    /// Returns `None` when the series is empty.
    pub fn bounding_box(&self) -> Option<(f64, f64, f64, f64)> {
        if self.points.is_empty() {
            return None;
        }
        let mut xmin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        let mut ymin = f64::INFINITY;
        let mut ymax = f64::NEG_INFINITY;
        for p in &self.points {
            if p.x < xmin {
                xmin = p.x;
            }
            if p.x > xmax {
                xmax = p.x;
            }
            if p.y < ymin {
                ymin = p.y;
            }
            if p.y > ymax {
                ymax = p.y;
            }
        }
        Some((xmin, xmax, ymin, ymax))
    }
}

// ---------------------------------------------------------------------------
// Stability classification
// ---------------------------------------------------------------------------

/// Stability classification of a fixed point in a 2-D ODE system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixedPointStability {
    /// Both eigenvalues have negative real parts with zero imaginary part.
    StableNode,
    /// Both eigenvalues have negative real parts with non-zero imaginary part.
    StableSpiral,
    /// Eigenvalues are real with opposite signs.
    Saddle,
    /// Both eigenvalues have positive real parts with zero imaginary part.
    UnstableNode,
    /// Both eigenvalues have positive real parts with non-zero imaginary part.
    UnstableSpiral,
    /// Both eigenvalues are purely imaginary (centre).
    Centre,
    /// Cannot be determined (degenerate or non-hyperbolic).
    Unknown,
}

/// A detected fixed point together with its stability.
#[derive(Debug, Clone)]
pub struct FixedPoint {
    /// The `(x, y)` location of the fixed point.
    pub location: PlotPoint,
    /// Stability classification based on Jacobian eigenvalues.
    pub stability: FixedPointStability,
    /// Jacobian matrix at the fixed point `[[a, b\], [c, d]]`.
    pub jacobian: [[f64; 2]; 2],
}

// ---------------------------------------------------------------------------
// PhasePortrait
// ---------------------------------------------------------------------------

/// A 2-D phase portrait: vector-field arrows on a regular grid, plus detected
/// fixed points with stability classification.
///
/// The user supplies two functions `fx(x, y)` and `fy(x, y)` representing
/// `dx/dt` and `dy/dt`.  After construction the grid is evaluated once and
/// fixed points can be found with [`PhasePortrait::find_fixed_points`].
#[derive(Debug, Clone)]
pub struct PhasePortrait {
    /// x-values of grid columns.
    pub x_values: Vec<f64>,
    /// y-values of grid rows.
    pub y_values: Vec<f64>,
    /// `dx/dt` at each `(col, row)` — indexed `[row * nx + col]`.
    pub dx: Vec<f64>,
    /// `dy/dt` at each `(col, row)` — indexed `[row * nx + col]`.
    pub dy: Vec<f64>,
    /// Fixed points found by Newton's method.
    pub fixed_points: Vec<FixedPoint>,
}

impl PhasePortrait {
    /// Build a phase portrait by evaluating `fx` and `fy` on a regular grid.
    ///
    /// - `x_range` — `(x_min, x_max)`
    /// - `y_range` — `(y_min, y_max)`
    /// - `nx`, `ny` — number of grid points per axis (clamped to ≥ 2)
    /// - `fx`, `fy` — ODE right-hand sides
    pub fn new(
        x_range: (f64, f64),
        y_range: (f64, f64),
        nx: usize,
        ny: usize,
        fx: impl Fn(f64, f64) -> f64,
        fy: impl Fn(f64, f64) -> f64,
    ) -> Self {
        let nx = nx.max(2);
        let ny = ny.max(2);
        let xs: Vec<f64> = linspace(x_range.0, x_range.1, nx);
        let ys: Vec<f64> = linspace(y_range.0, y_range.1, ny);
        let n = nx * ny;
        let mut dx_vec = Vec::with_capacity(n);
        let mut dy_vec = Vec::with_capacity(n);
        for &y in &ys {
            for &x in &xs {
                dx_vec.push(fx(x, y));
                dy_vec.push(fy(x, y));
            }
        }
        Self {
            x_values: xs,
            y_values: ys,
            dx: dx_vec,
            dy: dy_vec,
            fixed_points: Vec::new(),
        }
    }

    /// Search for fixed points near a set of seed guesses using Newton's method.
    ///
    /// Each seed is refined until `‖(fx, fy)‖ < tol` or `max_iter` is reached.
    /// Duplicates (within `dedup_eps`) are removed.  Fixed points are stored in
    /// `self.fixed_points` and also returned.
    pub fn find_fixed_points(
        &mut self,
        seeds: &[PlotPoint],
        fx: impl Fn(f64, f64) -> f64,
        fy: impl Fn(f64, f64) -> f64,
        tol: f64,
        max_iter: usize,
        dedup_eps: f64,
    ) -> &[FixedPoint] {
        let mut found: Vec<FixedPoint> = Vec::new();

        for seed in seeds {
            if let Some(fp) = newton_fixed_point(seed.x, seed.y, &fx, &fy, tol, max_iter) {
                // deduplication
                let is_dup = found
                    .iter()
                    .any(|e| e.location.distance_to(&fp.location) < dedup_eps);
                if !is_dup {
                    found.push(fp);
                }
            }
        }

        self.fixed_points = found;
        &self.fixed_points
    }

    /// Arrow vectors at grid cell `(col, row)`.  Returns `None` for out-of-bounds indices.
    pub fn arrow_at(&self, col: usize, row: usize) -> Option<(f64, f64)> {
        let nx = self.x_values.len();
        let ny = self.y_values.len();
        if col >= nx || row >= ny {
            return None;
        }
        let idx = row * nx + col;
        Some((self.dx[idx], self.dy[idx]))
    }
}

// ---------------------------------------------------------------------------
// BifurcationDiagram
// ---------------------------------------------------------------------------

/// A bifurcation diagram produced by sweeping a parameter `r` and recording
/// the long-term attractor of a scalar map `f(r, x)`.
///
/// Each entry in `attractors` corresponds to one value of `r` and contains the
/// last `tail_n` iterates after `transient` warm-up steps.
#[derive(Debug, Clone)]
pub struct BifurcationDiagram {
    /// Parameter values used in the sweep.
    pub r_values: Vec<f64>,
    /// For each `r_values[i]`, the attractor points `x` after transient removal.
    pub attractors: Vec<Vec<f64>>,
    /// Detected period at each `r` (unique values rounded to `period_tol`).
    pub periods: Vec<usize>,
}

impl BifurcationDiagram {
    /// Compute a bifurcation diagram for the 1-D map `f(r, x)`.
    ///
    /// - `r_range`  — `(r_min, r_max)`
    /// - `nr`       — number of parameter steps
    /// - `x0`       — initial condition
    /// - `transient`— number of iterates discarded as warm-up
    /// - `tail_n`   — number of iterates retained as attractor
    /// - `period_tol` — rounding tolerance for counting unique periods
    /// - `f`        — map: `f(r, x) → x_next`
    pub fn compute(
        r_range: (f64, f64),
        nr: usize,
        x0: f64,
        transient: usize,
        tail_n: usize,
        period_tol: f64,
        f: impl Fn(f64, f64) -> f64,
    ) -> Self {
        let nr = nr.max(1);
        let r_values: Vec<f64> = linspace(r_range.0, r_range.1, nr);
        let mut attractors = Vec::with_capacity(nr);
        let mut periods = Vec::with_capacity(nr);

        for &r in &r_values {
            let mut x = x0;
            for _ in 0..transient {
                x = f(r, x);
            }
            let mut tail = Vec::with_capacity(tail_n);
            for _ in 0..tail_n {
                x = f(r, x);
                tail.push(x);
            }
            let p = count_unique_periods(&tail, period_tol);
            periods.push(p);
            attractors.push(tail);
        }

        Self {
            r_values,
            attractors,
            periods,
        }
    }

    /// Return the flat list of `(r, x)` pairs suitable for scatter plotting.
    pub fn scatter_points(&self) -> Vec<PlotPoint> {
        let mut pts = Vec::new();
        for (i, r) in self.r_values.iter().enumerate() {
            for &x in &self.attractors[i] {
                pts.push(PlotPoint::new(*r, x));
            }
        }
        pts
    }

    /// Return the index of the first parameter value at which the period
    /// doubles (period 1 → 2, 2 → 4, etc.), or `None` if no doubling is found.
    pub fn first_period_doubling(&self) -> Option<usize> {
        for i in 1..self.periods.len() {
            let prev = self.periods[i - 1];
            let curr = self.periods[i];
            if prev > 0 && curr == 2 * prev {
                return Some(i);
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// PoincareSection
// ---------------------------------------------------------------------------

/// Definition of a hyperplane used to collect Poincaré-section crossings.
#[derive(Debug, Clone, Copy)]
pub struct SectionPlane {
    /// Which axis the section is taken along (0 = x, 1 = y, 2 = z).
    pub axis: usize,
    /// Coordinate value of the section on that axis.
    pub value: f64,
    /// If `true`, only up-crossings (negative → positive) are recorded.
    pub positive_crossing_only: bool,
}

/// A Poincaré section: intersections of a trajectory with a hyper-plane,
/// projected to the two axes orthogonal to the section axis.
#[derive(Debug, Clone)]
pub struct PoincareSection {
    /// The defining plane.
    pub plane: SectionPlane,
    /// Crossing points in 2-D coordinates `(u, v)` where `u` and `v` are the
    /// two axes orthogonal to `plane.axis`.
    pub crossings: Vec<PlotPoint>,
}

impl PoincareSection {
    /// Construct an empty Poincaré section for the given plane.
    pub fn new(plane: SectionPlane) -> Self {
        Self {
            plane,
            crossings: Vec::new(),
        }
    }

    /// Process a 3-D trajectory `[(x, y, z)]` and collect crossings.
    ///
    /// A crossing is recorded when the trajectory passes through `plane.value`
    /// on `plane.axis`.  The crossing location is linearly interpolated.
    pub fn collect_crossings(&mut self, trajectory: &[[f64; 3]]) {
        if trajectory.len() < 2 {
            return;
        }
        let ax = self.plane.axis.min(2);
        let val = self.plane.value;

        for window in trajectory.windows(2) {
            let a = window[0];
            let b = window[1];
            let da = a[ax] - val;
            let db = b[ax] - val;
            // Detect sign change
            if da * db >= 0.0 {
                continue;
            }
            if self.plane.positive_crossing_only && da >= 0.0 {
                // Only up-crossings (da < 0 → db > 0)
                continue;
            }
            // Linear interpolation parameter
            let t = da.abs() / (da.abs() + db.abs());
            let interp: [f64; 3] = [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ];
            // Project to the two non-section axes
            let (u, v) = orthogonal_coords(interp, ax);
            self.crossings.push(PlotPoint::new(u, v));
        }
    }

    /// Number of recorded crossings.
    pub fn count(&self) -> usize {
        self.crossings.len()
    }
}

// ---------------------------------------------------------------------------
// LyapunovSpectrum
// ---------------------------------------------------------------------------

/// Maximal Lyapunov exponent computed via the Benettin orbit-divergence method.
///
/// This implementation supports:
/// - Arbitrary 1-D maps via [`LyapunovSpectrum::from_map_1d`]
/// - The logistic map `x_{n+1} = r x (1 - x)` via [`LyapunovSpectrum::logistic`]
#[derive(Debug, Clone)]
pub struct LyapunovSpectrum {
    /// The estimated maximal Lyapunov exponent λ₁.
    pub lambda_max: f64,
    /// Running sum of `ln|f'(x)|` at each step (for diagnostics).
    pub log_sum_history: Vec<f64>,
}

impl LyapunovSpectrum {
    /// Estimate the maximal Lyapunov exponent for a 1-D map `f(x)` using
    /// the sum-of-log-derivatives method.
    ///
    /// λ₁ ≈ (1/N) Σ ln|f'(xₙ)|
    ///
    /// The derivative `df` is evaluated numerically with step `h` when not
    /// supplied analytically.
    ///
    /// - `x0`       — initial condition
    /// - `transient`— warm-up iterates (discarded)
    /// - `n_iter`   — iterates used for the average
    /// - `f`        — map `xₙ₊₁ = f(xₙ)`
    /// - `df`       — derivative `f'(xₙ)`; if `None`, finite difference with `h=1e-7`
    pub fn from_map_1d(
        x0: f64,
        transient: usize,
        n_iter: usize,
        f: impl Fn(f64) -> f64,
        df: Option<impl Fn(f64) -> f64>,
    ) -> Self {
        let h = 1e-7_f64;
        let has_df = df.is_some();
        let deriv: Box<dyn Fn(f64) -> f64> = match df {
            Some(d) => Box::new(d),
            None => Box::new(move |_x: f64| 0.0), // placeholder, won't be called
        };

        let mut x = x0;
        for _ in 0..transient {
            x = f(x);
        }

        let mut log_sum = 0.0;
        let mut history = Vec::with_capacity(n_iter);
        for _ in 0..n_iter {
            let slope = if has_df {
                deriv(x).abs()
            } else {
                ((f(x + h) - f(x - h)) / (2.0 * h)).abs()
            };
            if slope > 1e-300 {
                log_sum += slope.ln();
            }
            history.push(log_sum);
            x = f(x);
        }

        let lambda_max = if n_iter > 0 {
            log_sum / n_iter as f64
        } else {
            0.0
        };

        Self {
            lambda_max,
            log_sum_history: history,
        }
    }

    /// Compute the maximal Lyapunov exponent for the logistic map
    /// `f(x) = r·x·(1−x)` with derivative `f'(x) = r·(1−2x)`.
    ///
    /// - `r`        — growth parameter
    /// - `x0`       — initial condition (should be in `(0, 1)`)
    /// - `transient`— warm-up iterates
    /// - `n_iter`   — iterates used for the average
    pub fn logistic(r: f64, x0: f64, transient: usize, n_iter: usize) -> Self {
        Self::from_map_1d(
            x0,
            transient,
            n_iter,
            move |x| r * x * (1.0 - x),
            Some(move |x| r * (1.0 - 2.0 * x)),
        )
    }

    /// Return `true` when the system is likely chaotic (λ₁ > `threshold`).
    pub fn is_chaotic(&self, threshold: f64) -> bool {
        self.lambda_max > threshold
    }
}

// ---------------------------------------------------------------------------
// AttractorData
// ---------------------------------------------------------------------------

/// Storage for a 3-D trajectory (e.g. Lorenz or Rössler attractor) plus a
/// crude correlation-dimension estimate D₂.
#[derive(Debug, Clone)]
pub struct AttractorData {
    /// Trajectory points `[x, y, z]`.
    pub points: Vec<[f64; 3]>,
}

impl AttractorData {
    /// Create an empty attractor container.
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    /// Append a trajectory point.
    pub fn push(&mut self, x: f64, y: f64, z: f64) {
        self.points.push([x, y, z]);
    }

    /// Integrate the Lorenz system
    /// `dx/dt = σ(y−x)`, `dy/dt = x(ρ−z)−y`, `dz/dt = xy−βz`
    /// using RK4 and store the trajectory.
    ///
    /// - `sigma`, `rho`, `beta` — Lorenz parameters
    /// - `x0`, `y0`, `z0`      — initial condition
    /// - `dt`                  — time step
    /// - `n_steps`             — number of steps
    pub fn lorenz(
        sigma: f64,
        rho: f64,
        beta: f64,
        x0: f64,
        y0: f64,
        z0: f64,
        dt: f64,
        n_steps: usize,
    ) -> Self {
        let mut data = Self::new();
        let mut state = [x0, y0, z0];
        data.push(state[0], state[1], state[2]);

        let lorenz_rhs = |s: [f64; 3]| -> [f64; 3] {
            [
                sigma * (s[1] - s[0]),
                s[0] * (rho - s[2]) - s[1],
                s[0] * s[1] - beta * s[2],
            ]
        };

        for _ in 0..n_steps {
            state = rk4_step_3d(state, dt, lorenz_rhs);
            data.push(state[0], state[1], state[2]);
        }
        data
    }

    /// Integrate the Rössler system
    /// `dx/dt = −y−z`, `dy/dt = x+ay`, `dz/dt = b+z(x−c)`
    /// using RK4 and store the trajectory.
    pub fn rossler(
        a: f64,
        b: f64,
        c: f64,
        x0: f64,
        y0: f64,
        z0: f64,
        dt: f64,
        n_steps: usize,
    ) -> Self {
        let mut data = Self::new();
        let mut state = [x0, y0, z0];
        data.push(state[0], state[1], state[2]);

        let rhs = move |s: [f64; 3]| -> [f64; 3] {
            [-s[1] - s[2], s[0] + a * s[1], b + s[2] * (s[0] - c)]
        };

        for _ in 0..n_steps {
            state = rk4_step_3d(state, dt, rhs);
            data.push(state[0], state[1], state[2]);
        }
        data
    }

    /// Estimate the correlation dimension D₂ using the Grassberger–Procaccia
    /// method on a random sub-sample of `sample` point pairs.
    ///
    /// D₂ ≈ d ln C(ε) / d ln ε evaluated between `eps_small` and `eps_large`.
    ///
    /// Returns `None` if fewer than 10 points are stored.
    pub fn correlation_dimension(
        &self,
        eps_small: f64,
        eps_large: f64,
        sample: usize,
    ) -> Option<f64> {
        let n = self.points.len();
        if n < 10 {
            return None;
        }
        // Use evenly spaced indices as a deterministic sub-sample
        let step = (n / sample.min(n)).max(1);
        let pts: Vec<[f64; 3]> = self.points.iter().step_by(step).copied().collect();
        let m = pts.len();
        if m < 4 {
            return None;
        }

        let count_pairs = |eps: f64| -> usize {
            let eps2 = eps * eps;
            let mut cnt = 0usize;
            for i in 0..m {
                for j in (i + 1)..m {
                    let d2 = dist2_3d(pts[i], pts[j]);
                    if d2 < eps2 {
                        cnt += 1;
                    }
                }
            }
            cnt
        };

        let c_small = count_pairs(eps_small) as f64;
        let c_large = count_pairs(eps_large) as f64;

        if c_small < 1.0 || c_large < 1.0 || c_large <= c_small {
            return None;
        }

        let d2 = (c_large.ln() - c_small.ln()) / (eps_large.ln() - eps_small.ln());
        Some(d2)
    }

    /// Number of stored trajectory points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Returns `true` when the trajectory is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

impl Default for AttractorData {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TimeFrequencyData  (STFT)
// ---------------------------------------------------------------------------

/// A short-time Fourier transform (STFT) spectrogram.
///
/// Computed with a Hann window, configurable hop size and FFT length.
/// Magnitudes (not power) are stored: `|X[k]|`.
#[derive(Debug, Clone)]
pub struct TimeFrequencyData {
    /// Center-times of each frame (seconds, given `sample_rate`).
    pub times: Vec<f64>,
    /// Frequency bins (Hz) — length `fft_size / 2 + 1`.
    pub frequencies: Vec<f64>,
    /// Magnitude spectrogram: `magnitude[frame][bin]`.
    pub magnitude: Vec<Vec<f64>>,
    /// FFT size used.
    pub fft_size: usize,
    /// Sample rate (Hz).
    pub sample_rate: f64,
}

impl TimeFrequencyData {
    /// Compute a STFT magnitude spectrogram.
    ///
    /// - `signal`      — input time-domain samples
    /// - `fft_size`    — DFT length (also window length; clamped to ≥ 2)
    /// - `hop_size`    — hop between successive frames (clamped to ≥ 1)
    /// - `sample_rate` — sample rate in Hz (used only for labelling)
    pub fn compute(signal: &[f64], fft_size: usize, hop_size: usize, sample_rate: f64) -> Self {
        let fft_size = fft_size.max(2);
        let hop_size = hop_size.max(1);
        let n_bins = fft_size / 2 + 1;

        // Hann window
        let window: Vec<f64> = hann_window(fft_size);

        // Frequency axis
        let frequencies: Vec<f64> = (0..n_bins)
            .map(|k| k as f64 * sample_rate / fft_size as f64)
            .collect();

        let mut times = Vec::new();
        let mut magnitude = Vec::new();

        let mut frame_start = 0usize;
        while frame_start + fft_size <= signal.len() {
            let frame_center_sample = frame_start + fft_size / 2;
            times.push(frame_center_sample as f64 / sample_rate);

            // Apply window
            let windowed: Vec<f64> = (0..fft_size)
                .map(|i| signal[frame_start + i] * window[i])
                .collect();

            // DFT (naive O(N²) — acceptable for test/demo sizes)
            let frame_mag = naive_dft_magnitude(&windowed, n_bins);
            magnitude.push(frame_mag);

            frame_start += hop_size;
        }

        Self {
            times,
            frequencies,
            magnitude,
            fft_size,
            sample_rate,
        }
    }

    /// Number of time frames.
    pub fn n_frames(&self) -> usize {
        self.times.len()
    }

    /// Number of frequency bins.
    pub fn n_bins(&self) -> usize {
        self.frequencies.len()
    }

    /// Peak frequency (Hz) in frame `frame_idx`.
    ///
    /// Returns `None` if the frame index is out of range.
    pub fn peak_frequency(&self, frame_idx: usize) -> Option<f64> {
        let mags = self.magnitude.get(frame_idx)?;
        let (peak_bin, _) = mags
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        self.frequencies.get(peak_bin).copied()
    }
}

// ---------------------------------------------------------------------------
// ScientificAxis
// ---------------------------------------------------------------------------

/// Scale type for a plot axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisScale {
    /// Evenly-spaced linear ticks.
    Linear,
    /// Logarithmically-spaced ticks (data must be > 0).
    Log,
}

/// A formatted tick label.
#[derive(Debug, Clone)]
pub struct TickLabel {
    /// The data-space value at which the tick is placed.
    pub value: f64,
    /// Human-readable label string.
    pub label: String,
}

/// A plot axis with linear or log scale, and Wilkinson-style tick generation.
#[derive(Debug, Clone)]
pub struct ScientificAxis {
    /// Minimum data value.
    pub data_min: f64,
    /// Maximum data value.
    pub data_max: f64,
    /// Scale type.
    pub scale: AxisScale,
    /// Axis label (e.g. "Time (s)").
    pub label: String,
}

impl ScientificAxis {
    /// Construct a new axis.
    pub fn new(data_min: f64, data_max: f64, scale: AxisScale, label: impl Into<String>) -> Self {
        Self {
            data_min,
            data_max,
            scale,
            label: label.into(),
        }
    }

    /// Generate tick positions and labels.
    ///
    /// Uses a simplified Wilkinson extended algorithm for linear axes.
    /// Log axes produce one tick per decade.
    ///
    /// - `target_count` — desired number of ticks (hint only; clamped to ≥ 2)
    pub fn generate_ticks(&self, target_count: usize) -> Vec<TickLabel> {
        let target_count = target_count.max(2);
        match self.scale {
            AxisScale::Linear => self.linear_ticks(target_count),
            AxisScale::Log => self.log_ticks(),
        }
    }

    fn linear_ticks(&self, target: usize) -> Vec<TickLabel> {
        let range = self.data_max - self.data_min;
        if range <= 0.0 || !range.is_finite() {
            return vec![TickLabel {
                value: self.data_min,
                label: format_tick(self.data_min),
            }];
        }

        // Wilkinson "nice number" step
        let raw_step = range / (target - 1) as f64;
        let step = nice_step(raw_step);

        let first = (self.data_min / step).floor() * step;
        let last = (self.data_max / step).ceil() * step;

        let mut ticks = Vec::new();
        let mut v = first;
        let eps = step * 1e-9;
        while v <= last + eps {
            if v >= self.data_min - eps && v <= self.data_max + eps {
                ticks.push(TickLabel {
                    value: v,
                    label: format_tick(v),
                });
            }
            v += step;
            // Guard against floating-point drift creating infinite loops
            if ticks.len() > 1000 {
                break;
            }
        }
        ticks
    }

    fn log_ticks(&self) -> Vec<TickLabel> {
        if self.data_min <= 0.0 || self.data_max <= 0.0 {
            return Vec::new();
        }
        let log_min = self.data_min.log10().floor() as i32;
        let log_max = self.data_max.log10().ceil() as i32;
        (log_min..=log_max)
            .map(|exp| {
                let v = 10f64.powi(exp);
                TickLabel {
                    value: v,
                    label: format!("10^{exp}"),
                }
            })
            .filter(|t| t.value >= self.data_min && t.value <= self.data_max)
            .collect()
    }

    /// Map a data value to normalised display coordinate `[0, 1]`.
    pub fn normalise(&self, value: f64) -> f64 {
        match self.scale {
            AxisScale::Linear => (value - self.data_min) / (self.data_max - self.data_min),
            AxisScale::Log => {
                if value <= 0.0 || self.data_min <= 0.0 {
                    return 0.0;
                }
                (value.ln() - self.data_min.ln()) / (self.data_max.ln() - self.data_min.ln())
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Evenly-spaced vector from `start` to `end` with `n` points.
fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    if n == 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (end - start) * i as f64 / (n - 1) as f64)
        .collect()
}

/// Attempt Newton's method starting from `(x0, y0)` for `fx(x,y) = 0`, `fy(x,y) = 0`.
/// Returns a [`FixedPoint`] if converged within `tol` and `max_iter` iterations.
fn newton_fixed_point(
    x0: f64,
    y0: f64,
    fx: &impl Fn(f64, f64) -> f64,
    fy: &impl Fn(f64, f64) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<FixedPoint> {
    let h = 1e-6_f64;
    let mut x = x0;
    let mut y = y0;

    for _ in 0..max_iter {
        let fxv = fx(x, y);
        let fyv = fy(x, y);
        if (fxv * fxv + fyv * fyv).sqrt() < tol {
            let j = numerical_jacobian(x, y, fx, fy, h);
            let stab = classify_stability(j);
            return Some(FixedPoint {
                location: PlotPoint::new(x, y),
                stability: stab,
                jacobian: j,
            });
        }
        // Jacobian for Newton step
        let j = numerical_jacobian(x, y, fx, fy, h);
        // Solve J * delta = -F  (2×2 system)
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() < 1e-14 {
            break;
        }
        let inv = [
            [j[1][1] / det, -j[0][1] / det],
            [-j[1][0] / det, j[0][0] / det],
        ];
        let dx = -(inv[0][0] * fxv + inv[0][1] * fyv);
        let dy = -(inv[1][0] * fxv + inv[1][1] * fyv);
        x += dx;
        y += dy;
    }
    // Check convergence at the last iterate
    let fxv = fx(x, y);
    let fyv = fy(x, y);
    if (fxv * fxv + fyv * fyv).sqrt() < tol {
        let j = numerical_jacobian(x, y, fx, fy, h);
        let stab = classify_stability(j);
        Some(FixedPoint {
            location: PlotPoint::new(x, y),
            stability: stab,
            jacobian: j,
        })
    } else {
        None
    }
}

/// Compute the 2×2 Jacobian `[[∂fx/∂x, ∂fx/∂y\], [∂fy/∂x, ∂fy/∂y]]` numerically.
fn numerical_jacobian(
    x: f64,
    y: f64,
    fx: &impl Fn(f64, f64) -> f64,
    fy: &impl Fn(f64, f64) -> f64,
    h: f64,
) -> [[f64; 2]; 2] {
    let dfdx_x = (fx(x + h, y) - fx(x - h, y)) / (2.0 * h);
    let dfdx_y = (fx(x, y + h) - fx(x, y - h)) / (2.0 * h);
    let dfdy_x = (fy(x + h, y) - fy(x - h, y)) / (2.0 * h);
    let dfdy_y = (fy(x, y + h) - fy(x, y - h)) / (2.0 * h);
    [[dfdx_x, dfdx_y], [dfdy_x, dfdy_y]]
}

/// Classify stability from a 2×2 Jacobian matrix.
fn classify_stability(j: [[f64; 2]; 2]) -> FixedPointStability {
    // Characteristic equation: λ² − tr·λ + det = 0
    let tr = j[0][0] + j[1][1];
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    let disc = tr * tr - 4.0 * det;

    if disc >= 0.0 {
        // Real eigenvalues
        let sqrt_disc = disc.sqrt();
        let l1 = (tr + sqrt_disc) / 2.0;
        let l2 = (tr - sqrt_disc) / 2.0;
        if l1 * l2 < 0.0 {
            FixedPointStability::Saddle
        } else if l1 < 0.0 && l2 < 0.0 {
            FixedPointStability::StableNode
        } else if l1 > 0.0 && l2 > 0.0 {
            FixedPointStability::UnstableNode
        } else {
            FixedPointStability::Unknown
        }
    } else {
        // Complex eigenvalues; real part = tr / 2
        let real_part = tr / 2.0;
        if real_part.abs() < 1e-10 {
            FixedPointStability::Centre
        } else if real_part < 0.0 {
            FixedPointStability::StableSpiral
        } else {
            FixedPointStability::UnstableSpiral
        }
    }
}

/// Count approximate unique periods in a tail by rounding values.
fn count_unique_periods(tail: &[f64], tol: f64) -> usize {
    if tail.is_empty() {
        return 0;
    }
    let mut unique: Vec<f64> = Vec::new();
    for &v in tail {
        let is_new = unique.iter().all(|&u| (v - u).abs() > tol);
        if is_new {
            unique.push(v);
        }
    }
    unique.len()
}

/// Return the two coordinates orthogonal to `axis` in a 3-D point.
fn orthogonal_coords(p: [f64; 3], axis: usize) -> (f64, f64) {
    match axis {
        0 => (p[1], p[2]),
        1 => (p[0], p[2]),
        _ => (p[0], p[1]),
    }
}

/// Single RK4 step for a 3-D ODE.
fn rk4_step_3d(state: [f64; 3], dt: f64, rhs: impl Fn([f64; 3]) -> [f64; 3]) -> [f64; 3] {
    let k1 = rhs(state);
    let k2 = rhs(add3(state, scale3(k1, dt * 0.5)));
    let k3 = rhs(add3(state, scale3(k2, dt * 0.5)));
    let k4 = rhs(add3(state, scale3(k3, dt)));
    [
        state[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        state[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        state[2] + dt / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
    ]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dist2_3d(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Hann window of length `n`.
fn hann_window(n: usize) -> Vec<f64> {
    use std::f64::consts::PI;
    (0..n)
        .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos()))
        .collect()
}

/// Naive DFT computing `|X[k]|` for `k = 0..n_bins`.
fn naive_dft_magnitude(signal: &[f64], n_bins: usize) -> Vec<f64> {
    let n = signal.len();
    (0..n_bins)
        .map(|k| {
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            for (i, &s) in signal.iter().enumerate() {
                let angle = -2.0 * PI * k as f64 * i as f64 / n as f64;
                re += s * angle.cos();
                im += s * angle.sin();
            }
            (re * re + im * im).sqrt()
        })
        .collect()
}

/// Round a step to the nearest "nice" number: 1, 2, 5 × 10^k.
fn nice_step(raw: f64) -> f64 {
    if raw <= 0.0 || !raw.is_finite() {
        return 1.0;
    }
    let exp = raw.log10().floor();
    let pow = 10f64.powf(exp);
    let frac = raw / pow;
    let nice_frac = if frac < 1.5 {
        1.0
    } else if frac < 3.5 {
        2.0
    } else if frac < 7.5 {
        5.0
    } else {
        10.0
    };
    nice_frac * pow
}

/// Format a tick value as a short decimal string.
fn format_tick(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let mag = v.abs().log10().floor() as i32;
    if (-3..=4).contains(&mag) {
        // Fixed notation with enough decimal places
        let decimals = ((-mag + 1).max(0)) as usize;
        format!("{:.prec$}", v, prec = decimals)
    } else {
        format!("{:.2e}", v)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- PlotPoint ----

    #[test]
    fn test_plot_point_distance() {
        let a = PlotPoint::new(0.0, 0.0);
        let b = PlotPoint::new(3.0, 4.0);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_plot_point_same_is_zero() {
        let a = PlotPoint::new(1.0, 2.0);
        assert!(a.distance_to(&a) < 1e-15);
    }

    // ---- PlotSeries ----

    #[test]
    fn test_plot_series_empty() {
        let s = PlotSeries::new("test", [1.0, 0.0, 0.0, 1.0]);
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
        assert!(s.bounding_box().is_none());
    }

    #[test]
    fn test_plot_series_bounding_box() {
        let mut s = PlotSeries::new("s", [1.0, 1.0, 1.0, 1.0]);
        s.push(-1.0, 0.0);
        s.push(2.0, 3.0);
        s.push(0.5, -0.5);
        let bb = s.bounding_box().unwrap();
        assert!((bb.0 - (-1.0)).abs() < 1e-12); // xmin
        assert!((bb.1 - 2.0).abs() < 1e-12); // xmax
        assert!((bb.2 - (-0.5)).abs() < 1e-12); // ymin
        assert!((bb.3 - 3.0).abs() < 1e-12); // ymax
    }

    // ---- PhasePortrait / fixed points ----

    #[test]
    fn test_phase_portrait_grid_size() {
        let pp = PhasePortrait::new((-1.0, 1.0), (-1.0, 1.0), 5, 5, |x, _y| x, |_x, y| y);
        assert_eq!(pp.x_values.len(), 5);
        assert_eq!(pp.y_values.len(), 5);
        assert_eq!(pp.dx.len(), 25);
    }

    #[test]
    fn test_phase_portrait_arrow_at_bounds() {
        let pp = PhasePortrait::new((0.0, 1.0), (0.0, 1.0), 3, 3, |x, _y| x, |_x, y| y);
        assert!(pp.arrow_at(0, 0).is_some());
        assert!(pp.arrow_at(5, 0).is_none());
    }

    #[test]
    fn test_fixed_point_detection_logistic_growth() {
        // ODE: dx/dt = x(1-x), dy/dt = -y
        // Fixed points: (0,0) and (1,0)
        let fx = |x: f64, _y: f64| x * (1.0 - x);
        let fy = |_x: f64, y: f64| -y;

        let mut pp = PhasePortrait::new((-0.2, 1.2), (-0.5, 0.5), 8, 8, fx, fy);
        let seeds = vec![PlotPoint::new(0.1, 0.1), PlotPoint::new(0.9, 0.1)];
        let fps = pp.find_fixed_points(&seeds, fx, fy, 1e-10, 100, 1e-4);
        assert_eq!(fps.len(), 2);
        let xs: Vec<f64> = fps.iter().map(|fp| fp.location.x).collect();
        assert!(xs.iter().any(|&x| x.abs() < 1e-6), "should find FP at x≈0");
        assert!(
            xs.iter().any(|&x| (x - 1.0).abs() < 1e-6),
            "should find FP at x≈1"
        );
    }

    #[test]
    fn test_fixed_point_stability_saddle() {
        // Linear system: dx/dt = x, dy/dt = -y → saddle at origin
        let fx = |x: f64, _y: f64| x;
        let fy = |_x: f64, y: f64| -y;
        let j = numerical_jacobian(0.0, 0.0, &fx, &fy, 1e-6);
        let stab = classify_stability(j);
        assert_eq!(stab, FixedPointStability::Saddle);
    }

    #[test]
    fn test_fixed_point_stability_stable_node() {
        // Linear: dx/dt = -2x, dy/dt = -3y → stable node
        let fx = |x: f64, _y: f64| -2.0 * x;
        let fy = |_x: f64, y: f64| -3.0 * y;
        let j = numerical_jacobian(0.0, 0.0, &fx, &fy, 1e-6);
        let stab = classify_stability(j);
        assert_eq!(stab, FixedPointStability::StableNode);
    }

    #[test]
    fn test_fixed_point_stability_unstable_node() {
        // Linear: dx/dt = 2x, dy/dt = 3y → unstable node
        let fx = |x: f64, _y: f64| 2.0 * x;
        let fy = |_x: f64, y: f64| 3.0 * y;
        let j = numerical_jacobian(0.0, 0.0, &fx, &fy, 1e-6);
        let stab = classify_stability(j);
        assert_eq!(stab, FixedPointStability::UnstableNode);
    }

    #[test]
    fn test_fixed_point_stability_stable_spiral() {
        // dx/dt = -x - y, dy/dt = x - y → stable spiral (eigenvalues -1 ± i)
        let fx = |x: f64, y: f64| -x - y;
        let fy = |x: f64, y: f64| x - y;
        let j = numerical_jacobian(0.0, 0.0, &fx, &fy, 1e-6);
        let stab = classify_stability(j);
        assert_eq!(stab, FixedPointStability::StableSpiral);
    }

    // ---- Lyapunov spectrum ----

    #[test]
    fn test_lyapunov_chaotic_logistic_r39() {
        // r=3.9: logistic map is chaotic, λ > 0
        let spec = LyapunovSpectrum::logistic(3.9, 0.5, 1000, 5000);
        assert!(
            spec.lambda_max > 0.0,
            "r=3.9 should be chaotic, got λ={}",
            spec.lambda_max
        );
        assert!(spec.is_chaotic(0.0));
    }

    #[test]
    fn test_lyapunov_stable_logistic_r25() {
        // r=2.5: converges to fixed point x*=(r-1)/r, λ < 0
        let spec = LyapunovSpectrum::logistic(2.5, 0.5, 1000, 5000);
        assert!(
            spec.lambda_max < 0.0,
            "r=2.5 should be stable, got λ={}",
            spec.lambda_max
        );
        assert!(!spec.is_chaotic(0.0));
    }

    #[test]
    fn test_lyapunov_history_length() {
        let n = 200usize;
        let spec = LyapunovSpectrum::logistic(3.9, 0.5, 100, n);
        assert_eq!(spec.log_sum_history.len(), n);
    }

    #[test]
    fn test_lyapunov_stable_r1() {
        // r=1.0: all orbits go to 0 → strongly stable
        let spec = LyapunovSpectrum::logistic(1.0, 0.3, 100, 1000);
        assert!(spec.lambda_max < 0.0, "r=1.0 should be stable");
    }

    // ---- Poincaré section ----

    #[test]
    fn test_poincare_section_crossings_z_plane() {
        // Helix: t in [0, 4π], x=cos(t), y=sin(t), z=t/(2π)
        // Section at z=0.5 (z axis=2): 2 crossings in [0, 4π]
        let n = 2000usize;
        let traj: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = 4.0 * PI * i as f64 / (n - 1) as f64;
                [t.cos(), t.sin(), t / (2.0 * PI)]
            })
            .collect();

        let plane = SectionPlane {
            axis: 2,
            value: 0.5,
            positive_crossing_only: false,
        };
        let mut ps = PoincareSection::new(plane);
        ps.collect_crossings(&traj);
        // Should detect 2 crossings (z goes 0→2, crosses 0.5 twice: once up, once... still up)
        // Actually z is monotone increasing, so 1 crossing
        assert!(ps.count() >= 1, "expected at least 1 crossing");
    }

    #[test]
    fn test_poincare_section_oscillating_trajectory() {
        // z(t) = sin(2πt/T): 4 zero-crossings in 2 full cycles
        let n = 1000usize;
        let traj: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = 2.0 * i as f64 / (n - 1) as f64; // t in [0, 2]
                [0.0, 0.0, (2.0 * PI * t).sin()]
            })
            .collect();

        let plane = SectionPlane {
            axis: 2,
            value: 0.0,
            positive_crossing_only: false,
        };
        let mut ps = PoincareSection::new(plane);
        ps.collect_crossings(&traj);
        // sin has 4 sign changes in 2 full periods
        assert!(
            ps.count() >= 2,
            "expected multiple crossings, got {}",
            ps.count()
        );
    }

    #[test]
    fn test_poincare_positive_crossing_only() {
        // Only record up-crossings (neg → pos)
        let traj: Vec<[f64; 3]> = vec![
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],  // up-crossing at z=0
            [0.0, 0.0, -1.0], // down-crossing — should be skipped
            [0.0, 0.0, 1.0],  // up-crossing
        ];
        let plane = SectionPlane {
            axis: 2,
            value: 0.0,
            positive_crossing_only: true,
        };
        let mut ps = PoincareSection::new(plane);
        ps.collect_crossings(&traj);
        assert_eq!(ps.count(), 2, "expected 2 up-crossings, got {}", ps.count());
    }

    // ---- BifurcationDiagram ----

    #[test]
    fn test_bifurcation_diagram_period_one_at_r25() {
        // r=2.5: single stable fixed point → period 1
        let diag = BifurcationDiagram::compute((2.5, 2.5), 1, 0.5, 2000, 64, 1e-6, |r, x| {
            r * x * (1.0 - x)
        });
        assert_eq!(diag.periods[0], 1, "r=2.5 should have period 1");
    }

    #[test]
    fn test_bifurcation_period_doubling_detected() {
        // Sweep r=2.5 → 4.0: period doubling should be detected somewhere
        let diag = BifurcationDiagram::compute((2.5, 4.0), 50, 0.5, 2000, 128, 1e-5, |r, x| {
            r * x * (1.0 - x)
        });
        // Should find at least one period doubling (period 1 → 2)
        let idx = diag.first_period_doubling();
        assert!(
            idx.is_some(),
            "period doubling should be detected in r=[2.5, 4.0]"
        );
    }

    #[test]
    fn test_bifurcation_scatter_points_count() {
        let tail_n = 32usize;
        let nr = 10usize;
        let diag = BifurcationDiagram::compute((3.0, 3.5), nr, 0.5, 500, tail_n, 1e-6, |r, x| {
            r * x * (1.0 - x)
        });
        let pts = diag.scatter_points();
        assert_eq!(pts.len(), nr * tail_n);
    }

    // ---- AttractorData ----

    #[test]
    fn test_lorenz_trajectory_length() {
        let n_steps = 500usize;
        let data = AttractorData::lorenz(10.0, 28.0, 8.0 / 3.0, 0.1, 0.0, 0.0, 0.01, n_steps);
        assert_eq!(data.len(), n_steps + 1);
    }

    #[test]
    fn test_rossler_trajectory_grows() {
        let data = AttractorData::rossler(0.2, 0.2, 5.7, 0.1, 0.0, 0.0, 0.01, 500);
        assert!(!data.is_empty());
        // Check the trajectory has diverged from initial condition
        let last = data.points.last().unwrap();
        let first = data.points[0];
        let d2 = dist2_3d(*last, first);
        assert!(
            d2 > 1e-6,
            "Rössler trajectory should evolve away from initial point"
        );
    }

    #[test]
    fn test_attractor_correlation_dimension_lorenz() {
        let data = AttractorData::lorenz(10.0, 28.0, 8.0 / 3.0, 1.0, 0.0, 0.0, 0.01, 2000);
        let d2 = data.correlation_dimension(0.5, 5.0, 200);
        assert!(d2.is_some(), "correlation dimension should be computable");
        let d = d2.unwrap();
        // Lorenz D2 ≈ 2.05; allow wide range for small sample
        assert!(d > 1.0 && d < 4.0, "unexpected D2={}", d);
    }

    #[test]
    fn test_attractor_correlation_dimension_needs_10_points() {
        let mut data = AttractorData::new();
        for i in 0..5 {
            data.push(i as f64, 0.0, 0.0);
        }
        assert!(data.correlation_dimension(0.1, 1.0, 50).is_none());
    }

    // ---- TimeFrequencyData / STFT ----

    #[test]
    fn test_stft_pure_tone_peak_frequency() {
        // Generate a pure 100 Hz sine at 1 kHz sample rate
        let fs = 1000.0_f64;
        let f0 = 100.0_f64;
        let n = 512usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * f0 * i as f64 / fs).sin())
            .collect();
        let tf = TimeFrequencyData::compute(&signal, 256, 128, fs);
        let peak = tf.peak_frequency(0).unwrap();
        // Allow ±1 bin tolerance
        let bin_width = fs / 256.0;
        assert!(
            (peak - f0).abs() <= bin_width + 1e-9,
            "peak freq {} should be near {}; bin_width={}",
            peak,
            f0,
            bin_width
        );
    }

    #[test]
    fn test_stft_dc_signal_peak_at_zero() {
        let signal = vec![1.0_f64; 256];
        let tf = TimeFrequencyData::compute(&signal, 64, 32, 1000.0);
        let peak = tf.peak_frequency(0).unwrap();
        assert!(peak < 1.0, "DC signal peak should be near 0 Hz, got {peak}");
    }

    #[test]
    fn test_stft_frame_count() {
        let signal = vec![0.0_f64; 512];
        let fft = 64usize;
        let hop = 32usize;
        let tf = TimeFrequencyData::compute(&signal, fft, hop, 1.0);
        let expected = (512 - fft) / hop + 1;
        assert_eq!(tf.n_frames(), expected);
    }

    #[test]
    fn test_stft_bin_count() {
        let signal = vec![0.0_f64; 256];
        let fft = 64usize;
        let tf = TimeFrequencyData::compute(&signal, fft, 32, 1000.0);
        assert_eq!(tf.n_bins(), fft / 2 + 1);
    }

    // ---- ScientificAxis ----

    #[test]
    fn test_linear_axis_ticks_in_range() {
        let ax = ScientificAxis::new(0.0, 10.0, AxisScale::Linear, "x");
        let ticks = ax.generate_ticks(6);
        assert!(!ticks.is_empty());
        for t in &ticks {
            assert!(
                t.value >= 0.0 - 1e-9 && t.value <= 10.0 + 1e-9,
                "tick {} out of [0, 10]",
                t.value
            );
        }
    }

    #[test]
    fn test_log_axis_ticks_decades() {
        let ax = ScientificAxis::new(1.0, 1000.0, AxisScale::Log, "freq");
        let ticks = ax.generate_ticks(4);
        // Should have 10^0=1, 10^1=10, 10^2=100, 10^3=1000
        assert_eq!(
            ticks.len(),
            4,
            "expected 4 decade ticks, got {:?}",
            ticks.iter().map(|t| t.value).collect::<Vec<_>>()
        );
        assert!((ticks[0].value - 1.0).abs() < 1e-10);
        assert!((ticks[3].value - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_log_axis_labels_format() {
        let ax = ScientificAxis::new(1.0, 100.0, AxisScale::Log, "y");
        let ticks = ax.generate_ticks(3);
        assert!(ticks.iter().any(|t| t.label.contains("10^")));
    }

    #[test]
    fn test_axis_normalise_linear() {
        let ax = ScientificAxis::new(0.0, 10.0, AxisScale::Linear, "x");
        assert!((ax.normalise(0.0)).abs() < 1e-12);
        assert!((ax.normalise(10.0) - 1.0).abs() < 1e-12);
        assert!((ax.normalise(5.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_axis_normalise_log() {
        let ax = ScientificAxis::new(1.0, 100.0, AxisScale::Log, "x");
        assert!((ax.normalise(1.0)).abs() < 1e-12);
        assert!((ax.normalise(100.0) - 1.0).abs() < 1e-12);
        // log(10)/log(100) = 0.5
        assert!((ax.normalise(10.0) - 0.5).abs() < 1e-9);
    }

    // ---- Helpers ----

    #[test]
    fn test_linspace_endpoints() {
        let v = linspace(0.0, 1.0, 5);
        assert_eq!(v.len(), 5);
        assert!((v[0]).abs() < 1e-15);
        assert!((v[4] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_hann_window_endpoints_near_zero() {
        let w = hann_window(64);
        assert!(w[0].abs() < 1e-10, "Hann window should start near 0");
        assert!(w[63].abs() < 1e-3, "Hann window should end near 0");
    }

    #[test]
    fn test_nice_step_round_numbers() {
        assert!((nice_step(0.9) - 1.0).abs() < 1e-10);
        assert!((nice_step(1.5) - 2.0).abs() < 1e-10);
        assert!((nice_step(4.0) - 5.0).abs() < 1e-10);
        assert!((nice_step(9.0) - 10.0).abs() < 1e-10);
    }
}
