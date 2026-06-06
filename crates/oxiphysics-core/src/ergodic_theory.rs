// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ergodic theory and dynamical systems analysis.
//!
//! Provides tools for analysing long-time statistical behaviour of dynamical
//! systems, including Lyapunov exponent estimation, Poincaré section extraction,
//! recurrence-plot analysis, mixing-rate estimation, Birkhoff ergodic averages,
//! invariant measure estimation, Kolmogorov–Sinai entropy, and topological entropy.

// ─── Invariant Measure / Ergodic Measure ─────────────────────────────────────

/// Empirical approximation of an invariant measure for a dynamical system.
///
/// The measure is represented as a histogram over a bounded interval.  Given a
/// long trajectory the histogram converges (by the ergodic theorem) to the
/// natural invariant measure of the system.
#[derive(Debug, Clone)]
pub struct ErgodicMeasure {
    /// Left endpoint of the histogram domain.
    pub x_min: f64,
    /// Right endpoint of the histogram domain.
    pub x_max: f64,
    /// Number of bins.
    pub num_bins: usize,
    /// Raw counts per bin.
    pub counts: Vec<u64>,
    /// Total samples ingested.
    pub total: u64,
}

impl ErgodicMeasure {
    /// Construct an empty measure over `[x_min, x_max]` with `num_bins` bins.
    ///
    /// # Panics
    /// Panics if `num_bins == 0` or if `x_min >= x_max`.
    pub fn new(x_min: f64, x_max: f64, num_bins: usize) -> Self {
        assert!(num_bins > 0, "num_bins must be > 0");
        assert!(x_min < x_max, "x_min must be < x_max");
        Self {
            x_min,
            x_max,
            num_bins,
            counts: vec![0u64; num_bins],
            total: 0,
        }
    }

    /// Add a single observation `x` to the histogram.
    ///
    /// Values outside `[x_min, x_max)` are silently dropped.
    pub fn update(&mut self, x: f64) {
        if x < self.x_min || x >= self.x_max {
            return;
        }
        let frac = (x - self.x_min) / (self.x_max - self.x_min);
        let bin = ((frac * self.num_bins as f64) as usize).min(self.num_bins - 1);
        self.counts[bin] += 1;
        self.total += 1;
    }

    /// Return the normalised density in each bin (counts / total / bin_width).
    ///
    /// Returns a zero vector if no samples have been ingested.
    pub fn density(&self) -> Vec<f64> {
        if self.total == 0 {
            return vec![0.0; self.num_bins];
        }
        let bin_width = (self.x_max - self.x_min) / self.num_bins as f64;
        self.counts
            .iter()
            .map(|&c| c as f64 / (self.total as f64 * bin_width))
            .collect()
    }

    /// Compute the empirical (time) average of an observable `f` under this measure.
    ///
    /// Uses the midpoint of each bin as the representative x value.
    /// Returns 0 if no samples have been ingested.
    pub fn space_average<F>(&self, f: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        if self.total == 0 {
            return 0.0;
        }
        let bin_width = (self.x_max - self.x_min) / self.num_bins as f64;
        let mut acc = 0.0_f64;
        for (k, &c) in self.counts.iter().enumerate() {
            let mid = self.x_min + (k as f64 + 0.5_f64) * bin_width;
            acc += f(mid) * c as f64;
        }
        acc / self.total as f64
    }

    /// Compare time average and space average of an observable `f` along a trajectory.
    ///
    /// Returns `(time_average, space_average, absolute_difference)`.
    pub fn compare_averages<F>(&self, trajectory: &[f64], f: &F) -> (f64, f64, f64)
    where
        F: Fn(f64) -> f64,
    {
        let time_avg = if trajectory.is_empty() {
            0.0_f64
        } else {
            trajectory.iter().map(|&x| f(x)).sum::<f64>() / trajectory.len() as f64
        };
        let space_avg = self.space_average(f);
        let diff = (time_avg - space_avg).abs();
        (time_avg, space_avg, diff)
    }
}

// ─── Birkhoff Average ────────────────────────────────────────────────────────

/// Running Birkhoff time-average of an observable along a trajectory.
///
/// Tracks how the time average converges as more of the trajectory is consumed.
/// For an ergodic system, the running average converges to the space average
/// under the invariant measure.
#[derive(Debug, Clone)]
pub struct BirkhoffAverage {
    /// The accumulated sum of `f(x_t)` observed so far.
    pub sum: f64,
    /// Number of samples consumed.
    pub count: usize,
    /// History of running averages (one entry per `update` call).
    pub history: Vec<f64>,
}

impl BirkhoffAverage {
    /// Create a new, empty Birkhoff accumulator.
    pub fn new() -> Self {
        Self {
            sum: 0.0,
            count: 0,
            history: Vec::new(),
        }
    }

    /// Consume a single value `fx = f(x_t)` and append to history.
    pub fn update(&mut self, fx: f64) {
        self.sum += fx;
        self.count += 1;
        self.history.push(self.sum / self.count as f64);
    }

    /// Return the current time average (returns 0 if no samples consumed).
    pub fn current(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.sum / self.count as f64
        }
    }

    /// Ingest an entire trajectory slice, applying observable `f`.
    pub fn consume_trajectory<F>(&mut self, trajectory: &[f64], f: &F)
    where
        F: Fn(f64) -> f64,
    {
        for &x in trajectory {
            self.update(f(x));
        }
    }
}

impl Default for BirkhoffAverage {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute the Birkhoff ergodic (time) average of a discrete signal.
///
/// Equivalent to `signal.iter().sum::`f64`() / signal.len()` but returns 0
/// for an empty signal.
pub fn birkhoff_average(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    signal.iter().sum::<f64>() / signal.len() as f64
}

// ─── Lyapunov Spectrum ────────────────────────────────────────────────────────

/// Analyser for the full Lyapunov spectrum of a dynamical system via QR
/// decomposition (Gram-Schmidt orthogonalisation) of the accumulated Jacobian.
///
/// Lyapunov exponents characterise the average rate of divergence (or
/// convergence) of nearby trajectories.  A positive leading exponent is a
/// necessary condition for chaos.
pub struct LyapunovSpectrum {
    /// Dimension of the phase space.
    pub dim: usize,
    /// Accumulated log-stretching rates for each mode.
    pub log_sums: Vec<f64>,
    /// Number of QR steps taken.
    pub steps: usize,
}

impl LyapunovSpectrum {
    /// Construct a new spectrum analyser for a phase space of dimension `dim`.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            log_sums: vec![0.0_f64; dim],
            steps: 0,
        }
    }

    /// Perform one QR step given the current `dim × dim` Jacobian matrix
    /// stored in row-major order.
    ///
    /// The Gram-Schmidt QR decomposition extracts the R diagonal, whose
    /// log-absolute-values are accumulated into the Lyapunov sums.
    ///
    /// # Panics
    /// Panics if `jacobian.len() != dim * dim`.
    pub fn update_qr(&mut self, jacobian: &[f64]) {
        let d = self.dim;
        assert_eq!(jacobian.len(), d * d, "jacobian size mismatch");

        // Extract column vectors from row-major storage.
        let mut cols: Vec<Vec<f64>> = (0..d)
            .map(|j| (0..d).map(|i| jacobian[i * d + j]).collect())
            .collect();

        // Gram-Schmidt with R diagonal extraction.
        for j in 0..d {
            // Compute R[j,j] = ||cols[j]||.
            let norm = cols[j].iter().map(|&v| v * v).sum::<f64>().sqrt();
            if norm > 1e-300 {
                self.log_sums[j] += norm.ln();
                // Normalise column j.
                for v in cols[j].iter_mut() {
                    *v /= norm;
                }
            }
            // Project out component along cols[j] from remaining columns.
            for k in (j + 1)..d {
                let dot: f64 = cols[j]
                    .iter()
                    .zip(cols[k].iter())
                    .map(|(&a, &b)| a * b)
                    .sum();
                let cj = cols[j].clone();
                for (val, cj_val) in cols[k].iter_mut().zip(cj.iter()) {
                    *val -= dot * cj_val;
                }
            }
        }
        self.steps += 1;
    }

    /// Return the current Lyapunov exponents (per unit time `dt`).
    ///
    /// Divides each accumulated sum by total elapsed time `steps * dt`.
    /// Returns zeros if no steps have been taken.
    pub fn exponents(&self, dt: f64) -> Vec<f64> {
        let total = self.steps as f64 * dt;
        if total < 1e-300 {
            return vec![0.0_f64; self.dim];
        }
        self.log_sums.iter().map(|&s| s / total).collect()
    }
}

/// Legacy analyser wrapper kept for backwards compatibility.
pub struct LyapunovAnalyzer {
    /// Dimension of the phase space.
    pub dim: usize,
}

impl LyapunovAnalyzer {
    /// Construct a new analyser for a phase space of dimension `dim`.
    pub fn new(dim: usize) -> Self {
        Self { dim }
    }
}

/// Estimate Lyapunov exponents from a finite-time trajectory segment.
///
/// The trajectory is interpreted as a sequence of state vectors of equal
/// length.  Successive displacement vectors are accumulated into log-stretching
/// rates.  Returns exponents sorted in descending order.
///
/// # Panics
/// Panics if `trajectory` is empty or if `dt` ≤ 0.
pub fn lyapunov_exponents(trajectory: &[Vec<f64>], dt: f64) -> Vec<f64> {
    assert!(!trajectory.is_empty(), "trajectory must not be empty");
    assert!(dt > 0.0, "dt must be positive");

    let n = trajectory[0].len();
    if n == 0 || trajectory.len() < 2 {
        return vec![0.0; n];
    }

    let steps = trajectory.len() - 1;
    let mut sums = vec![0.0_f64; n];

    for t in 0..steps {
        let prev = &trajectory[t];
        let next = &trajectory[t + 1];
        for ((s, nxt), prv) in sums.iter_mut().zip(next.iter()).zip(prev.iter()) {
            let delta = (nxt - prv).abs();
            if delta > 1e-300 {
                *s += delta.ln();
            }
        }
    }

    let total_time = steps as f64 * dt;
    let mut exponents: Vec<f64> = sums.iter().map(|&s| s / total_time).collect();
    exponents.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    exponents
}

/// Estimate the maximal Lyapunov exponent of the logistic map `x → r·x·(1-x)`.
///
/// Uses the analytic formula: `λ = (1/N) Σ ln|r - 2·r·x_t|`.
/// For `r = 4` the value converges to `ln 2 ≈ 0.693`.
pub fn logistic_lyapunov(r: f64, x0: f64, n: usize) -> f64 {
    let mut x = x0;
    let mut sum = 0.0_f64;
    for _ in 0..n {
        let deriv = (r - 2.0_f64 * r * x).abs();
        if deriv > 1e-300 {
            sum += deriv.ln();
        }
        x = r * x * (1.0_f64 - x);
    }
    if n == 0 { 0.0_f64 } else { sum / n as f64 }
}

// ─── Poincaré Section ────────────────────────────────────────────────────────

/// Poincaré section recorder.
///
/// Accumulates trajectory crossings through a hyperplane defined by a normal
/// vector and a point on the plane.
#[derive(Debug, Clone)]
pub struct PoincareSection {
    /// Outward normal of the section plane (need not be unit-length).
    pub plane_normal: [f64; 3],
    /// A point lying on the section plane.
    pub plane_point: [f64; 3],
    /// Recorded crossing points (positions where the trajectory pierced the plane).
    pub crossings: Vec<[f64; 3]>,
}

impl PoincareSection {
    /// Construct a new (empty) Poincaré section from a normal and a point.
    pub fn new(normal: [f64; 3], point: [f64; 3]) -> Self {
        Self {
            plane_normal: normal,
            plane_point: point,
            crossings: Vec::new(),
        }
    }
}

/// Find all upward crossings of a plane by a 3-D trajectory.
///
/// A crossing occurs when the signed distance to the plane changes from
/// negative to positive between consecutive samples.  The crossing position is
/// linearly interpolated.
pub fn poincare_section(
    trajectory: &[[f64; 3]],
    normal: [f64; 3],
    point: [f64; 3],
) -> Vec<[f64; 3]> {
    if trajectory.len() < 2 {
        return Vec::new();
    }

    let dot3 = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let sub3 = |a: &[f64; 3], b: &[f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let signed_dist = |p: &[f64; 3]| dot3(&sub3(p, &point), &normal);

    let mut crossings = Vec::new();
    let mut prev_dist = signed_dist(&trajectory[0]);

    for i in 1..trajectory.len() {
        let curr_dist = signed_dist(&trajectory[i]);
        if prev_dist < 0.0 && curr_dist >= 0.0 {
            let t = if (curr_dist - prev_dist).abs() < 1e-15 {
                0.5_f64
            } else {
                (-prev_dist) / (curr_dist - prev_dist)
            };
            let p = &trajectory[i - 1];
            let q = &trajectory[i];
            let crossing = [
                p[0] + t * (q[0] - p[0]),
                p[1] + t * (q[1] - p[1]),
                p[2] + t * (q[2] - p[2]),
            ];
            crossings.push(crossing);
        }
        prev_dist = curr_dist;
    }
    crossings
}

// ─── Recurrence Analysis ─────────────────────────────────────────────────────

/// A recurrence plot (RP) computed from a scalar time series.
///
/// Entry `matrix[i][j]` is `true` when `|x[i] - x[j]| <= threshold`.
#[derive(Debug, Clone)]
pub struct RecurrencePlot {
    /// Side length of the (square) recurrence matrix.
    pub size: usize,
    /// Distance threshold below which two states are considered recurrent.
    pub threshold: f64,
    /// Boolean recurrence matrix.
    pub matrix: Vec<Vec<bool>>,
}

impl RecurrencePlot {
    /// Build a recurrence plot from a scalar signal.
    pub fn new(signal: &[f64], threshold: f64) -> Self {
        let n = signal.len();
        let mut matrix = vec![vec![false; n]; n];
        for i in 0..n {
            for j in 0..n {
                matrix[i][j] = (signal[i] - signal[j]).abs() <= threshold;
            }
        }
        Self {
            size: n,
            threshold,
            matrix,
        }
    }
}

/// Recurrence quantification analysis (RQA) results.
///
/// Summarises the statistical structure of a recurrence plot via standard
/// RQA measures.
#[derive(Debug, Clone)]
pub struct RecurrenceAnalysis {
    /// Recurrence rate (RR): fraction of recurrent points.
    pub recurrence_rate: f64,
    /// Determinism (DET): fraction of recurrent points on diagonal lines ≥ min_line.
    pub determinism: f64,
    /// Mean diagonal line length (L_mean).
    pub mean_diagonal_length: f64,
    /// Longest diagonal line length (L_max).
    pub max_diagonal_length: usize,
    /// Entropy of diagonal line-length distribution (Shannon).
    pub diagonal_entropy: f64,
    /// Laminarity (LAM): fraction of recurrent points on vertical lines ≥ min_line.
    pub laminarity: f64,
    /// Trapping time (TT): mean vertical line length.
    pub trapping_time: f64,
}

impl RecurrenceAnalysis {
    /// Compute full RQA from a recurrence plot with minimum line length `min_line`.
    ///
    /// `min_line` is the minimum length (≥ 1) for a line to be counted.
    pub fn compute(rp: &RecurrencePlot, min_line: usize) -> Self {
        let rr = recurrence_rate(rp);
        let det = determinism(rp);
        let total_recurrent = (rp.size * rp.size) as f64 * rr;

        // Diagonal line statistics.
        let diag_lengths = diagonal_line_lengths(rp, min_line);
        let total_on_diag: usize = diag_lengths.iter().copied().sum();
        let mean_diag = if diag_lengths.is_empty() {
            0.0_f64
        } else {
            total_on_diag as f64 / diag_lengths.len() as f64
        };
        let max_diag = diag_lengths.iter().cloned().max().unwrap_or(0);

        // Shannon entropy of diagonal line lengths.
        let diag_entropy = if diag_lengths.is_empty() {
            0.0_f64
        } else {
            let n_lines = diag_lengths.len() as f64;
            // Build histogram of lengths.
            let max_l = *diag_lengths.iter().max().unwrap_or(&1);
            let mut hist = vec![0usize; max_l + 1];
            for &l in &diag_lengths {
                hist[l] += 1;
            }
            hist.iter()
                .filter(|&&c| c > 0)
                .map(|&c| {
                    let p = c as f64 / n_lines;
                    -p * p.ln()
                })
                .sum()
        };

        // Vertical line statistics (laminarity / trapping time).
        let vert_lengths = vertical_line_lengths(rp, min_line);
        let total_on_vert: usize = vert_lengths.iter().copied().sum();
        let lam = if total_recurrent < 1.0 {
            0.0_f64
        } else {
            total_on_vert as f64 / total_recurrent
        };
        let trapping = if vert_lengths.is_empty() {
            0.0_f64
        } else {
            total_on_vert as f64 / vert_lengths.len() as f64
        };

        Self {
            recurrence_rate: rr,
            determinism: det,
            mean_diagonal_length: mean_diag,
            max_diagonal_length: max_diag,
            diagonal_entropy: diag_entropy,
            laminarity: lam,
            trapping_time: trapping,
        }
    }
}

/// Collect diagonal line lengths ≥ `min_line` from a recurrence plot.
fn diagonal_line_lengths(rp: &RecurrencePlot, min_line: usize) -> Vec<usize> {
    let n = rp.size;
    let mut lengths = Vec::new();
    for k in -(n as isize - 1)..=(n as isize - 1) {
        let mut run = 0usize;
        let i_start = if k < 0 { (-k) as usize } else { 0 };
        let i_end = if k < 0 {
            n
        } else {
            n.saturating_sub(k as usize)
        };
        for i in i_start..i_end {
            let j = (i as isize + k) as usize;
            if rp.matrix[i][j] {
                run += 1;
            } else {
                if run >= min_line {
                    lengths.push(run);
                }
                run = 0;
            }
        }
        if run >= min_line {
            lengths.push(run);
        }
    }
    lengths
}

/// Collect vertical line lengths ≥ `min_line` from a recurrence plot.
fn vertical_line_lengths(rp: &RecurrencePlot, min_line: usize) -> Vec<usize> {
    let n = rp.size;
    let mut lengths = Vec::new();
    for j in 0..n {
        let mut run = 0usize;
        for i in 0..n {
            if rp.matrix[i][j] {
                run += 1;
            } else {
                if run >= min_line {
                    lengths.push(run);
                }
                run = 0;
            }
        }
        if run >= min_line {
            lengths.push(run);
        }
    }
    lengths
}

/// Compute the recurrence rate: fraction of `true` entries in the RP.
///
/// Returns 0 if the plot is empty.
pub fn recurrence_rate(rp: &RecurrencePlot) -> f64 {
    if rp.size == 0 {
        return 0.0;
    }
    let total = (rp.size * rp.size) as f64;
    let recurrent: usize = rp
        .matrix
        .iter()
        .flat_map(|row| row.iter())
        .filter(|&&v| v)
        .count();
    recurrent as f64 / total
}

/// Compute determinism: fraction of recurrent points that lie on diagonal lines
/// of length ≥ 2.
///
/// Returns 0 if there are no recurrent points.
pub fn determinism(rp: &RecurrencePlot) -> f64 {
    if rp.size == 0 {
        return 0.0;
    }
    let n = rp.size;
    let mut on_diag = 0usize;
    let mut total_recurrent = 0usize;

    for i in 0..n {
        for j in 0..n {
            if rp.matrix[i][j] {
                total_recurrent += 1;
            }
        }
    }

    if total_recurrent == 0 {
        return 0.0;
    }

    for k in -(n as isize - 1)..=(n as isize - 1) {
        let mut run = 0usize;
        let i_start = if k < 0 { (-k) as usize } else { 0 };
        let i_end = if k < 0 { n } else { n - k as usize };
        for i in i_start..i_end {
            let j = (i as isize + k) as usize;
            if rp.matrix[i][j] {
                run += 1;
            } else {
                if run >= 2 {
                    on_diag += run;
                }
                run = 0;
            }
        }
        if run >= 2 {
            on_diag += run;
        }
    }

    on_diag as f64 / total_recurrent as f64
}

// ─── Entropic Analysis ───────────────────────────────────────────────────────

/// Entropy analysis of a dynamical system trajectory.
///
/// Provides estimates of Kolmogorov–Sinai (KS) entropy via the symbolic
/// partition method, and topological entropy via the growth rate of the
/// number of distinct symbolic words.
#[derive(Debug, Clone)]
pub struct EntropicAnalysis {
    /// Word length used for symbolic encoding.
    pub word_length: usize,
    /// Alphabet size (number of partition cells).
    pub alphabet_size: usize,
    /// Estimated KS entropy (bits per unit time).
    pub ks_entropy: f64,
    /// Estimated topological entropy.
    pub topological_entropy: f64,
}

impl EntropicAnalysis {
    /// Estimate KS entropy and topological entropy from a symbolic sequence.
    ///
    /// `symbols` is a sequence of partition-cell indices in `0..alphabet_size`.
    /// `word_length` is the window size for word counting (recommend ≥ 4).
    pub fn from_symbols(symbols: &[usize], word_length: usize, alphabet_size: usize) -> Self {
        let ks = kolmogorov_sinai_entropy(symbols, word_length);
        let topo = topological_entropy_estimate(symbols, word_length);
        Self {
            word_length,
            alphabet_size,
            ks_entropy: ks,
            topological_entropy: topo,
        }
    }
}

/// Estimate the Kolmogorov–Sinai entropy from a symbolic sequence.
///
/// Uses the approximation `h_KS ≈ H(L+1) - H(L)` where `H(L)` is the
/// block entropy of words of length `L`.
pub fn kolmogorov_sinai_entropy(symbols: &[usize], word_length: usize) -> f64 {
    if symbols.len() <= word_length + 1 || word_length == 0 {
        return 0.0_f64;
    }
    let h_l = block_entropy(symbols, word_length);
    let h_l1 = block_entropy(symbols, word_length + 1);
    (h_l1 - h_l).max(0.0_f64)
}

/// Compute the block entropy `H(L)` of a symbolic sequence.
///
/// `H(L) = -Σ p(w) log2(p(w))` where the sum is over all words `w` of length `L`.
pub fn block_entropy(symbols: &[usize], word_length: usize) -> f64 {
    use std::collections::HashMap;
    if symbols.len() < word_length || word_length == 0 {
        return 0.0_f64;
    }
    let mut counts: HashMap<Vec<usize>, usize> = HashMap::new();
    let num_words = symbols.len() - word_length + 1;
    for i in 0..num_words {
        let word = symbols[i..i + word_length].to_vec();
        *counts.entry(word).or_insert(0) += 1;
    }
    let total = num_words as f64;
    counts
        .values()
        .map(|&c| {
            let p = c as f64 / total;
            -p * p.log2()
        })
        .sum()
}

/// Estimate topological entropy as `log2(number of distinct words of length L) / L`.
pub fn topological_entropy_estimate(symbols: &[usize], word_length: usize) -> f64 {
    use std::collections::HashSet;
    if symbols.len() < word_length || word_length == 0 {
        return 0.0_f64;
    }
    let mut words: HashSet<Vec<usize>> = HashSet::new();
    let num_words = symbols.len() - word_length + 1;
    for i in 0..num_words {
        words.insert(symbols[i..i + word_length].to_vec());
    }
    let distinct = words.len() as f64;
    if distinct <= 1.0_f64 {
        0.0_f64
    } else {
        distinct.log2() / word_length as f64
    }
}

/// Encode a real-valued time series into symbols by uniform partitioning.
///
/// Partitions `[x_min, x_max)` into `num_symbols` cells and maps each
/// element of `signal` to its cell index.
pub fn symbolic_encode(signal: &[f64], x_min: f64, x_max: f64, num_symbols: usize) -> Vec<usize> {
    assert!(num_symbols > 0, "num_symbols must be > 0");
    let range = x_max - x_min;
    if range <= 0.0_f64 {
        return vec![0usize; signal.len()];
    }
    signal
        .iter()
        .map(|&x| {
            let frac = (x - x_min) / range;
            let idx = (frac * num_symbols as f64) as usize;
            idx.min(num_symbols - 1)
        })
        .collect()
}

// ─── Mixing Rate ─────────────────────────────────────────────────────────────

/// Estimates the mixing rate of a dynamical system from the decay of its
/// autocorrelation function.
///
/// An exponentially mixing system has autocorrelation `C(k) ~ exp(-k/τ)`;
/// this struct fits that exponential to estimate the mixing time `τ`.
#[derive(Debug, Clone)]
pub struct MixingRate {
    /// Autocorrelation values at successive lags (lag 0 = 1 by convention).
    pub correlation_decay: Vec<f64>,
    /// Estimated exponential decay rate `γ` such that `C(k) ≈ exp(-γ k)`.
    pub decay_rate: f64,
    /// Estimated mixing time `τ = 1/γ` (in units of lag).
    pub mixing_timescale: f64,
}

impl MixingRate {
    /// Construct from a pre-computed autocorrelation sequence and estimate
    /// the exponential decay rate via log-linear regression.
    pub fn from_autocorr(corr: Vec<f64>) -> Self {
        let rate = estimate_decay_rate(&corr);
        let timescale = if rate > 1e-300 {
            1.0_f64 / rate
        } else {
            f64::INFINITY
        };
        Self {
            correlation_decay: corr,
            decay_rate: rate,
            mixing_timescale: timescale,
        }
    }

    /// Estimate the mixing time as the first lag at which the autocorrelation
    /// drops below `threshold`.
    pub fn mixing_time(&self, threshold: f64) -> Option<usize> {
        mixing_time(&self.correlation_decay, threshold)
    }
}

/// Fit an exponential decay `C(k) = exp(-γ k)` to an autocorrelation sequence
/// by log-linear least squares regression.
///
/// Returns `γ ≥ 0`.
pub fn estimate_decay_rate(corr: &[f64]) -> f64 {
    if corr.len() < 2 {
        return 0.0_f64;
    }
    let mut sum_k = 0.0_f64;
    let mut sum_log_c = 0.0_f64;
    let mut sum_k2 = 0.0_f64;
    let mut sum_k_log_c = 0.0_f64;
    let mut count = 0.0_f64;
    for (k, &c) in corr.iter().enumerate() {
        if c > 0.0_f64 {
            let log_c = c.ln();
            let kf = k as f64;
            sum_k += kf;
            sum_log_c += log_c;
            sum_k2 += kf * kf;
            sum_k_log_c += kf * log_c;
            count += 1.0_f64;
        }
    }
    if count < 2.0_f64 {
        return 0.0_f64;
    }
    let slope = (count * sum_k_log_c - sum_k * sum_log_c) / (count * sum_k2 - sum_k * sum_k);
    (-slope).max(0.0_f64)
}

/// Legacy alias kept for API compatibility.
///
/// Estimates the mixing rate of a dynamical system from the decay of its
/// autocorrelation function.
#[derive(Debug, Clone)]
pub struct MixingAnalyzer {
    /// Autocorrelation values at successive lags (lag 0 = 1 by convention).
    pub correlation_decay: Vec<f64>,
}

impl MixingAnalyzer {
    /// Construct from a pre-computed autocorrelation sequence.
    pub fn new(correlation_decay: Vec<f64>) -> Self {
        Self { correlation_decay }
    }

    /// Estimate the mixing time as the first lag at which the autocorrelation
    /// drops below `threshold`.
    pub fn mixing_time(&self, threshold: f64) -> Option<usize> {
        mixing_time(&self.correlation_decay, threshold)
    }
}

/// Compute the normalised autocorrelation of `signal` up to `max_lag`.
///
/// The autocorrelation at lag `k` is `C(k) / C(0)` where
/// `C(k) = mean((x[t] - mu) * (x[t+k] - mu))`.  Returns a vector of length
/// `min(max_lag + 1, signal.len())`.
pub fn autocorrelation(signal: &[f64], max_lag: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return Vec::new();
    }
    let mu = signal.iter().sum::<f64>() / n as f64;
    let c0: f64 = signal.iter().map(|&x| (x - mu) * (x - mu)).sum::<f64>() / n as f64;
    if c0 < 1e-300 {
        return vec![1.0_f64; (max_lag + 1).min(n)];
    }

    let lags = (max_lag + 1).min(n);
    let mut result = Vec::with_capacity(lags);
    for k in 0..lags {
        let len = n - k;
        let ck: f64 = (0..len)
            .map(|t| (signal[t] - mu) * (signal[t + k] - mu))
            .sum::<f64>()
            / n as f64;
        result.push(ck / c0);
    }
    result
}

/// Return the first lag at which the autocorrelation drops below `threshold`.
///
/// Returns `None` if the autocorrelation never falls below the threshold.
pub fn mixing_time(autocorr: &[f64], threshold: f64) -> Option<usize> {
    autocorr.iter().position(|&v| v < threshold)
}

// ─── Ergodic Mean ────────────────────────────────────────────────────────────

/// Online estimator for the time-average ergodic mean and variance.
///
/// Accumulates samples one-by-one and computes the running mean and variance
/// using Welford's online algorithm.
#[derive(Debug, Clone)]
pub struct ErgodicMean {
    /// All accumulated sample values.
    pub values: Vec<f64>,
    /// Current running mean (Welford).
    mean_acc: f64,
    /// Current running M2 (for variance, Welford).
    m2_acc: f64,
}

impl ErgodicMean {
    /// Create a new, empty estimator.
    pub fn new() -> Self {
        Self {
            values: Vec::new(),
            mean_acc: 0.0,
            m2_acc: 0.0,
        }
    }

    /// Ingest a new sample, updating the running statistics.
    pub fn update(&mut self, x: f64) {
        self.values.push(x);
        let n = self.values.len() as f64;
        let delta = x - self.mean_acc;
        self.mean_acc += delta / n;
        let delta2 = x - self.mean_acc;
        self.m2_acc += delta * delta2;
    }

    /// Return the current time-average mean.
    ///
    /// Returns 0 if no samples have been ingested.
    pub fn mean(&self) -> f64 {
        self.mean_acc
    }

    /// Return the current sample variance (unbiased, Bessel-corrected).
    ///
    /// Returns 0 if fewer than two samples have been ingested.
    pub fn variance(&self) -> f64 {
        let n = self.values.len();
        if n < 2 {
            0.0
        } else {
            self.m2_acc / (n as f64 - 1.0)
        }
    }
}

impl Default for ErgodicMean {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Baker's Map ─────────────────────────────────────────────────────────────

/// Iterate the baker's map on the unit square `[0,1)²`.
///
/// The baker's map is: `(x, y) → (2x mod 1, (y + floor(2x)) / 2)`.
/// It is uniformly hyperbolic with Lyapunov exponents `±ln 2`.
pub fn bakers_map_step(x: f64, y: f64) -> (f64, f64) {
    let two_x = 2.0_f64 * x;
    let floor = two_x.floor();
    let xn = two_x - floor;
    let yn = (y + floor) / 2.0_f64;
    (xn, yn)
}

/// Generate a trajectory of the baker's map starting at `(x0, y0)`.
pub fn bakers_map_trajectory(x0: f64, y0: f64, n: usize) -> Vec<(f64, f64)> {
    let mut traj = Vec::with_capacity(n + 1);
    let mut x = x0;
    let mut y = y0;
    traj.push((x, y));
    for _ in 0..n {
        let (xn, yn) = bakers_map_step(x, y);
        x = xn;
        y = yn;
        traj.push((x, y));
    }
    traj
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- ErgodicMeasure ---

    #[test]
    fn test_ergodic_measure_update_and_density() {
        let mut em = ErgodicMeasure::new(0.0, 1.0, 10);
        for i in 0..100 {
            em.update(i as f64 / 100.0);
        }
        assert_eq!(em.total, 100);
        let dens = em.density();
        assert_eq!(dens.len(), 10);
        // Each bin should have roughly 10 samples → density ≈ 1.0 (uniform).
        for &d in &dens {
            assert!(d > 0.0, "density should be positive");
        }
    }

    #[test]
    fn test_ergodic_measure_out_of_bounds_dropped() {
        let mut em = ErgodicMeasure::new(0.0, 1.0, 5);
        em.update(-0.5);
        em.update(1.5);
        assert_eq!(em.total, 0);
    }

    #[test]
    fn test_ergodic_measure_density_zero_when_empty() {
        let em = ErgodicMeasure::new(0.0, 10.0, 5);
        let dens = em.density();
        for &d in &dens {
            assert_eq!(d, 0.0);
        }
    }

    #[test]
    fn test_ergodic_measure_space_average_constant() {
        let mut em = ErgodicMeasure::new(0.0, 1.0, 100);
        for i in 0..1000 {
            em.update(i as f64 / 1000.0);
        }
        // Space average of identity f(x) = x over uniform distribution ≈ 0.5.
        let avg = em.space_average(|x| x);
        assert!((avg - 0.5).abs() < 0.05, "space average ≈ 0.5, got {}", avg);
    }

    #[test]
    fn test_ergodic_measure_compare_averages() {
        let mut em = ErgodicMeasure::new(0.0, 1.0, 50);
        let traj: Vec<f64> = (0..500).map(|i| i as f64 / 500.0).collect();
        for &x in &traj {
            em.update(x);
        }
        let (_ta, _sa, diff) = em.compare_averages(&traj, &|x| x);
        assert!(
            diff < 0.1,
            "time and space averages should be close: diff = {}",
            diff
        );
    }

    #[test]
    fn test_ergodic_measure_panics_zero_bins() {
        let result = std::panic::catch_unwind(|| ErgodicMeasure::new(0.0, 1.0, 0));
        assert!(result.is_err());
    }

    // --- BirkhoffAverage ---

    #[test]
    fn test_birkhoff_average_empty_current() {
        let ba = BirkhoffAverage::new();
        assert_eq!(ba.current(), 0.0);
    }

    #[test]
    fn test_birkhoff_average_single_update() {
        let mut ba = BirkhoffAverage::new();
        ba.update(7.0);
        assert!((ba.current() - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_birkhoff_average_convergence() {
        let mut ba = BirkhoffAverage::new();
        for i in 1..=100 {
            ba.update(i as f64);
        }
        // Mean of 1..=100 = 50.5.
        assert!((ba.current() - 50.5).abs() < 1e-10);
    }

    #[test]
    fn test_birkhoff_average_history_length() {
        let mut ba = BirkhoffAverage::new();
        for _ in 0..20 {
            ba.update(1.0);
        }
        assert_eq!(ba.history.len(), 20);
    }

    #[test]
    fn test_birkhoff_average_consume_trajectory() {
        let mut ba = BirkhoffAverage::new();
        let traj: Vec<f64> = vec![2.0, 4.0, 6.0, 8.0];
        ba.consume_trajectory(&traj, &|x| x);
        // Mean = 5.0.
        assert!((ba.current() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_birkhoff_average_observable() {
        let mut ba = BirkhoffAverage::new();
        let traj: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        ba.consume_trajectory(&traj, &|x| x * x);
        // Mean of squares: (1+4+9+16)/4 = 7.5.
        assert!((ba.current() - 7.5).abs() < 1e-12);
    }

    #[test]
    fn test_birkhoff_average_default() {
        let ba = BirkhoffAverage::default();
        assert_eq!(ba.count, 0);
    }

    // --- LyapunovSpectrum ---

    #[test]
    fn test_lyapunov_spectrum_new() {
        let ls = LyapunovSpectrum::new(2);
        assert_eq!(ls.dim, 2);
        assert_eq!(ls.steps, 0);
    }

    #[test]
    fn test_lyapunov_spectrum_identity_jacobian() {
        let mut ls = LyapunovSpectrum::new(2);
        // Identity Jacobian: stretching rates = 0 for all modes.
        let identity = [1.0, 0.0, 0.0, 1.0];
        for _ in 0..10 {
            ls.update_qr(&identity);
        }
        let exp = ls.exponents(1.0);
        // log(1) / total_time = 0.
        for &e in &exp {
            assert!(
                e.abs() < 0.5,
                "exponent should be near 0 for identity: {}",
                e
            );
        }
    }

    #[test]
    fn test_lyapunov_spectrum_diagonal_stretching() {
        let mut ls = LyapunovSpectrum::new(2);
        // Diagonal Jacobian with eigenvalues 2 and 0.5.
        let diag = [2.0, 0.0, 0.0, 0.5];
        for _ in 0..100 {
            ls.update_qr(&diag);
        }
        let exp = ls.exponents(1.0);
        assert_eq!(exp.len(), 2);
        // Leading exponent should be ln(2) ≈ 0.693.
        assert!(
            exp[0] > 0.0,
            "first exponent should be positive: {}",
            exp[0]
        );
    }

    #[test]
    fn test_lyapunov_spectrum_exponents_zero_dt() {
        let ls = LyapunovSpectrum::new(3);
        let exp = ls.exponents(0.0);
        for &e in &exp {
            assert_eq!(e, 0.0);
        }
    }

    // --- logistic_lyapunov ---

    #[test]
    fn test_logistic_lyapunov_r4_positive() {
        // At r = 4 the Lyapunov exponent is ln(2) ≈ 0.693.
        let le = logistic_lyapunov(4.0, 0.3, 10_000);
        assert!(le > 0.5, "logistic LE at r=4 should be > 0.5, got {}", le);
        assert!(le < 1.0, "logistic LE at r=4 should be < 1.0, got {}", le);
    }

    #[test]
    fn test_logistic_lyapunov_r2_negative() {
        // At r = 2 the map has a stable fixed point → negative LE.
        let le = logistic_lyapunov(2.0, 0.4, 5_000);
        assert!(
            le < 0.0,
            "logistic LE at r=2 should be negative, got {}",
            le
        );
    }

    #[test]
    fn test_logistic_lyapunov_zero_steps() {
        assert_eq!(logistic_lyapunov(4.0, 0.3, 0), 0.0);
    }

    // --- LyapunovAnalyzer / lyapunov_exponents ---

    #[test]
    fn test_lyapunov_empty_trajectory_panics() {
        let result = std::panic::catch_unwind(|| lyapunov_exponents(&[], 0.01));
        assert!(result.is_err());
    }

    #[test]
    fn test_lyapunov_zero_dt_panics() {
        let traj = vec![vec![1.0, 0.0], vec![1.1, 0.0]];
        let result = std::panic::catch_unwind(|| lyapunov_exponents(&traj, 0.0));
        assert!(result.is_err());
    }

    #[test]
    fn test_lyapunov_single_step_returns_vector() {
        let traj = vec![vec![1.0, 0.0], vec![2.0, 0.0]];
        let exps = lyapunov_exponents(&traj, 1.0);
        assert_eq!(exps.len(), 2);
    }

    #[test]
    fn test_lyapunov_growing_trajectory_positive() {
        let n = 20;
        let dt = 1.0;
        let traj: Vec<Vec<f64>> = (0..=n).map(|t| vec![2.0_f64.powi(t)]).collect();
        let exps = lyapunov_exponents(&traj, dt);
        assert!(!exps.is_empty());
        assert!(
            exps[0] > 0.0,
            "leading exponent should be positive: {}",
            exps[0]
        );
    }

    #[test]
    fn test_lyapunov_constant_trajectory_near_zero() {
        let traj: Vec<Vec<f64>> = vec![vec![1.0, 2.0]; 10];
        let exps = lyapunov_exponents(&traj, 0.1);
        for e in &exps {
            assert!(e.abs() < 1.0, "exponent {} should be near zero", e);
        }
    }

    #[test]
    fn test_lyapunov_sorted_descending() {
        let traj: Vec<Vec<f64>> = (0..20)
            .map(|t| vec![2.0_f64.powi(t), 1.5_f64.powi(t), 0.9_f64.powi(t)])
            .collect();
        let exps = lyapunov_exponents(&traj, 1.0);
        for i in 1..exps.len() {
            assert!(
                exps[i - 1] >= exps[i],
                "exponents not sorted at index {}",
                i
            );
        }
    }

    #[test]
    fn test_lyapunov_analyzer_new() {
        let la = LyapunovAnalyzer::new(3);
        assert_eq!(la.dim, 3);
    }

    // --- PoincareSection / poincare_section ---

    #[test]
    fn test_poincare_empty_trajectory() {
        let crossings = poincare_section(&[], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(crossings.is_empty());
    }

    #[test]
    fn test_poincare_single_point_no_crossings() {
        let crossings = poincare_section(&[[0.0, 0.0, 1.0]], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(crossings.is_empty());
    }

    #[test]
    fn test_poincare_one_crossing_through_z_plane() {
        let traj = [[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]];
        let crossings = poincare_section(&traj, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0][2]).abs() < 1e-10, "crossing z should be ~0");
    }

    #[test]
    fn test_poincare_multiple_crossings() {
        let n = 200;
        let traj: Vec<[f64; 3]> = (0..n)
            .map(|i| [0.0, 0.0, (2.0 * PI * i as f64 / 100.0 + PI).sin()])
            .collect();
        let crossings = poincare_section(&traj, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(!crossings.is_empty(), "expected at least one crossing");
    }

    #[test]
    fn test_poincare_section_struct_new() {
        let ps = PoincareSection::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(ps.crossings.is_empty());
    }

    #[test]
    fn test_poincare_downward_not_counted() {
        let traj = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
        let crossings = poincare_section(&traj, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(
            crossings.is_empty(),
            "downward crossing should not be counted"
        );
    }

    // --- RecurrencePlot & RQA ---

    #[test]
    fn test_recurrence_rate_constant_signal() {
        let sig = vec![5.0; 10];
        let rp = RecurrencePlot::new(&sig, 0.1);
        let rate = recurrence_rate(&rp);
        assert!((rate - 1.0).abs() < 1e-10, "rate = {}", rate);
    }

    #[test]
    fn test_recurrence_rate_sparse() {
        let sig: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let rp = RecurrencePlot::new(&sig, 0.5);
        let rate = recurrence_rate(&rp);
        assert!((rate - 0.1).abs() < 1e-10, "rate = {}", rate);
    }

    #[test]
    fn test_recurrence_rate_bounds() {
        let sig: Vec<f64> = (0..20).map(|i| (i as f64).sin()).collect();
        let rp = RecurrencePlot::new(&sig, 0.3);
        let rate = recurrence_rate(&rp);
        assert!((0.0..=1.0).contains(&rate), "rate out of bounds: {}", rate);
    }

    #[test]
    fn test_recurrence_rate_empty() {
        let rp = RecurrencePlot {
            size: 0,
            threshold: 0.1,
            matrix: Vec::new(),
        };
        assert_eq!(recurrence_rate(&rp), 0.0);
    }

    #[test]
    fn test_determinism_constant_signal() {
        let sig = vec![1.0; 5];
        let rp = RecurrencePlot::new(&sig, 0.01);
        let det = determinism(&rp);
        assert!(det > 0.0, "determinism should be > 0 for constant signal");
    }

    #[test]
    fn test_determinism_bounds() {
        let sig: Vec<f64> = (0..15).map(|i| (i as f64 * 0.5).sin()).collect();
        let rp = RecurrencePlot::new(&sig, 0.3);
        let det = determinism(&rp);
        assert!(
            (0.0..=1.0).contains(&det),
            "determinism out of bounds: {}",
            det
        );
    }

    #[test]
    fn test_determinism_empty() {
        let rp = RecurrencePlot {
            size: 0,
            threshold: 0.1,
            matrix: Vec::new(),
        };
        assert_eq!(determinism(&rp), 0.0);
    }

    #[test]
    fn test_rqa_compute_sine_wave() {
        let sig: Vec<f64> = (0..50)
            .map(|i| (2.0 * PI * i as f64 / 25.0).sin())
            .collect();
        let rp = RecurrencePlot::new(&sig, 0.3);
        let rqa = RecurrenceAnalysis::compute(&rp, 2);
        assert!(rqa.recurrence_rate > 0.0);
        assert!(rqa.determinism >= 0.0 && rqa.determinism <= 1.0);
    }

    #[test]
    fn test_rqa_max_diagonal_length() {
        // Constant signal: full recurrence → longest diagonal = n.
        let sig = vec![1.0_f64; 8];
        let rp = RecurrencePlot::new(&sig, 0.01);
        let rqa = RecurrenceAnalysis::compute(&rp, 2);
        assert!(
            rqa.max_diagonal_length >= 6,
            "max diag = {}",
            rqa.max_diagonal_length
        );
    }

    // --- Entropic Analysis ---

    #[test]
    fn test_symbolic_encode_uniform() {
        let sig: Vec<f64> = (0..10).map(|i| i as f64 / 10.0).collect();
        let syms = symbolic_encode(&sig, 0.0, 1.0, 10);
        assert_eq!(syms.len(), 10);
        for (i, &s) in syms.iter().enumerate() {
            assert_eq!(s, i, "symbol mismatch at index {}", i);
        }
    }

    #[test]
    fn test_block_entropy_constant() {
        // Constant signal → single symbol → entropy = 0.
        let syms = vec![0usize; 20];
        let h = block_entropy(&syms, 3);
        assert!(h.abs() < 1e-10, "entropy of constant = {}", h);
    }

    #[test]
    fn test_block_entropy_uniform() {
        // Alternating symbols → non-zero entropy.
        let syms: Vec<usize> = (0..40).map(|i| i % 4).collect();
        let h = block_entropy(&syms, 2);
        assert!(h > 0.0, "entropy should be > 0 for varied signal");
    }

    #[test]
    fn test_ks_entropy_logistic_chaos() {
        // For the logistic map at r=4, we expect non-trivial KS entropy.
        let mut x = 0.3_f64;
        let mut traj = Vec::with_capacity(1000);
        for _ in 0..1000 {
            x = 4.0 * x * (1.0 - x);
            traj.push(x);
        }
        let syms = symbolic_encode(&traj, 0.0, 1.0, 8);
        let ks = kolmogorov_sinai_entropy(&syms, 4);
        assert!(ks >= 0.0, "KS entropy must be non-negative: {}", ks);
    }

    #[test]
    fn test_topological_entropy_constant() {
        let syms = vec![0usize; 30];
        let te = topological_entropy_estimate(&syms, 3);
        assert_eq!(te, 0.0, "topo entropy of constant signal should be 0");
    }

    #[test]
    fn test_entropic_analysis_from_symbols() {
        let syms: Vec<usize> = (0..100).map(|i| i % 4).collect();
        let ea = EntropicAnalysis::from_symbols(&syms, 3, 4);
        assert!(ea.ks_entropy >= 0.0);
        assert!(ea.topological_entropy >= 0.0);
    }

    // --- MixingRate ---

    #[test]
    fn test_autocorrelation_lag0_is_one() {
        let sig: Vec<f64> = (0..50).map(|i| (i as f64).sin()).collect();
        let ac = autocorrelation(&sig, 10);
        assert!(
            (ac[0] - 1.0).abs() < 1e-10,
            "lag-0 autocorrelation must be 1"
        );
    }

    #[test]
    fn test_autocorrelation_constant_signal() {
        let sig = vec![3.0; 20];
        let ac = autocorrelation(&sig, 5);
        for &v in &ac {
            assert!((v - 1.0).abs() < 1e-10, "constant signal AC should be 1");
        }
    }

    #[test]
    fn test_autocorrelation_sine_half_period() {
        let n = 200usize;
        let sig: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let ac = autocorrelation(&sig, n / 2);
        let ac_half = ac[n / 2];
        assert!(
            ac_half < -0.45,
            "AC at half-period should be strongly negative, got {}",
            ac_half
        );
    }

    #[test]
    fn test_autocorrelation_empty_signal() {
        let ac = autocorrelation(&[], 5);
        assert!(ac.is_empty());
    }

    #[test]
    fn test_autocorrelation_length() {
        let sig: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let ac = autocorrelation(&sig, 10);
        assert_eq!(ac.len(), 11);
    }

    #[test]
    fn test_mixing_time_found() {
        let ac = vec![1.0, 0.5, 0.2, 0.05, -0.01];
        assert_eq!(mixing_time(&ac, 0.1), Some(3));
    }

    #[test]
    fn test_mixing_time_not_found() {
        let ac = vec![1.0, 0.9, 0.8, 0.7];
        assert_eq!(mixing_time(&ac, 0.1), None);
    }

    #[test]
    fn test_mixing_time_immediately() {
        let ac = vec![0.0, 0.5, 1.0];
        assert_eq!(mixing_time(&ac, 0.1), Some(0));
    }

    #[test]
    fn test_mixing_analyzer_new() {
        let ma = MixingAnalyzer::new(vec![1.0, 0.5, 0.1]);
        assert_eq!(ma.correlation_decay.len(), 3);
        assert_eq!(ma.mixing_time(0.2), Some(2));
    }

    #[test]
    fn test_mixing_rate_from_exponential_decay() {
        // Exponential decay with rate γ = 0.5.
        let corr: Vec<f64> = (0..20).map(|k| (-0.5_f64 * k as f64).exp()).collect();
        let mr = MixingRate::from_autocorr(corr);
        // decay_rate should be close to 0.5.
        assert!(
            (mr.decay_rate - 0.5).abs() < 0.1,
            "decay rate = {}",
            mr.decay_rate
        );
        assert!(mr.mixing_timescale > 0.0);
    }

    // --- ErgodicMean ---

    #[test]
    fn test_ergodic_mean_empty() {
        let em = ErgodicMean::new();
        assert_eq!(em.mean(), 0.0);
        assert_eq!(em.variance(), 0.0);
    }

    #[test]
    fn test_ergodic_mean_single() {
        let mut em = ErgodicMean::new();
        em.update(5.0);
        assert!((em.mean() - 5.0).abs() < 1e-12);
        assert_eq!(em.variance(), 0.0);
    }

    #[test]
    fn test_ergodic_mean_convergence() {
        let mut em = ErgodicMean::new();
        for i in 1..=100 {
            em.update(i as f64);
        }
        assert!((em.mean() - 50.5).abs() < 1e-10, "mean = {}", em.mean());
    }

    #[test]
    fn test_ergodic_variance() {
        let mut em = ErgodicMean::new();
        for x in [1.0, 2.0, 3.0, 4.0, 5.0] {
            em.update(x);
        }
        assert!(
            (em.variance() - 2.5).abs() < 1e-10,
            "variance = {}",
            em.variance()
        );
    }

    #[test]
    fn test_ergodic_mean_default() {
        let em = ErgodicMean::default();
        assert_eq!(em.mean(), 0.0);
    }

    // --- birkhoff_average ---

    #[test]
    fn test_birkhoff_average_empty() {
        assert_eq!(birkhoff_average(&[]), 0.0);
    }

    #[test]
    fn test_birkhoff_average_constant() {
        let sig = vec![7.0; 100];
        assert!((birkhoff_average(&sig) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_birkhoff_average_integers() {
        let sig: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        assert!((birkhoff_average(&sig) - 5.5).abs() < 1e-12);
    }

    #[test]
    fn test_birkhoff_average_zero_mean() {
        let sig: Vec<f64> = (0..100)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        assert!(birkhoff_average(&sig).abs() < 1e-12);
    }

    // --- Baker's map ---

    #[test]
    fn test_bakers_map_step_in_unit_square() {
        let (xn, yn) = bakers_map_step(0.3, 0.7);
        assert!((0.0..1.0).contains(&xn), "x out of [0,1): {}", xn);
        assert!((0.0..1.0).contains(&yn), "y out of [0,1): {}", yn);
    }

    #[test]
    fn test_bakers_map_trajectory_length() {
        let traj = bakers_map_trajectory(0.2, 0.5, 50);
        assert_eq!(traj.len(), 51);
    }

    #[test]
    fn test_bakers_map_trajectory_first_point() {
        let traj = bakers_map_trajectory(0.3, 0.6, 10);
        assert!((traj[0].0 - 0.3).abs() < 1e-14);
        assert!((traj[0].1 - 0.6).abs() < 1e-14);
    }

    #[test]
    fn test_bakers_map_ergodic_mean_x() {
        // For the ergodic baker's map, all trajectory points must stay in [0,1).
        let traj = bakers_map_trajectory(0.123, 0.456, 5000);
        for (i, &(x, y)) in traj.iter().enumerate() {
            assert!(
                (0.0..1.0).contains(&x),
                "x out of [0,1) at step {}: {}",
                i,
                x
            );
            assert!(
                (0.0..1.0).contains(&y),
                "y out of [0,1) at step {}: {}",
                i,
                y
            );
        }
    }
}
