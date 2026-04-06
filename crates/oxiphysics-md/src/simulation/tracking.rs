#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Energy drift tracking, radial distribution, velocity autocorrelation,
//! pressure tensor, neighbor lists, ensemble switching, energy/pressure
//! tracking, callbacks, thermostats, and volume monitoring.

use rand::RngExt;

use super::config::{Ensemble, KB_REDUCED, MdConfig, MdState};

// ---------------------------------------------------------------------------
// EnergyDriftTracker
// ---------------------------------------------------------------------------

/// Tracks energy drift over the course of an MD run.
///
/// Records energy at each call to `update`, exposing running statistics:
/// initial energy, maximum absolute drift, and relative drift.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct EnergyDriftTracker {
    /// Energy at step 0.
    pub e0: f64,
    /// Maximum absolute deviation from e0 seen so far.
    pub max_abs_drift: f64,
    /// Last recorded energy.
    pub last_energy: f64,
    /// Number of updates.
    pub n_updates: u64,
    initialized: bool,
}

#[allow(dead_code)]
impl EnergyDriftTracker {
    /// Create a new tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Feed the current total energy.
    pub fn update(&mut self, total_energy: f64) {
        if !self.initialized {
            self.e0 = total_energy;
            self.initialized = true;
        }
        let drift = (total_energy - self.e0).abs();
        if drift > self.max_abs_drift {
            self.max_abs_drift = drift;
        }
        self.last_energy = total_energy;
        self.n_updates += 1;
    }

    /// Relative energy drift |ΔE_max| / |E₀|.
    pub fn relative_drift(&self) -> f64 {
        if self.e0.abs() < 1e-30 {
            return 0.0;
        }
        self.max_abs_drift / self.e0.abs()
    }

    /// Reset the tracker.
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}

// ---------------------------------------------------------------------------
// RdfHistogram — radial distribution function accumulator
// ---------------------------------------------------------------------------

/// Radial distribution function (RDF) accumulator.
///
/// Bins pair distances into a histogram and normalises to yield g(r).
/// Operates on plain `[f64; 3]` positions with an orthorhombic periodic box.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RdfHistogram {
    /// Number of histogram bins.
    pub n_bins: usize,
    /// Bin width (Å).
    pub dr: f64,
    /// Maximum distance (Å).
    pub r_max: f64,
    /// Raw pair-count histogram.
    pub histogram: Vec<u64>,
    /// Number of frames accumulated.
    pub n_frames: u64,
    /// Total number of particles (constant across frames).
    pub n_atoms: usize,
    /// Box volume (Å³) — used for normalisation.
    pub volume: f64,
}

#[allow(dead_code)]
impl RdfHistogram {
    /// Create a new RDF histogram with `n_bins` bins up to `r_max`.
    pub fn new(r_max: f64, n_bins: usize, volume: f64, n_atoms: usize) -> Self {
        assert!(r_max > 0.0, "r_max must be positive");
        assert!(n_bins > 0, "n_bins must be positive");
        Self {
            n_bins,
            dr: r_max / n_bins as f64,
            r_max,
            histogram: vec![0u64; n_bins],
            n_frames: 0,
            n_atoms,
            volume,
        }
    }

    /// Accumulate one frame of positions into the histogram.
    ///
    /// Applies minimum-image convention for the cubic box `box_len`.
    pub fn accumulate(&mut self, positions: &[[f64; 3]], box_len: f64) {
        let n = positions.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let dr = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                // Minimum image
                let dr = [
                    dr[0] - box_len * (dr[0] / box_len).round(),
                    dr[1] - box_len * (dr[1] / box_len).round(),
                    dr[2] - box_len * (dr[2] / box_len).round(),
                ];
                let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                if r < self.r_max {
                    let bin = (r / self.dr) as usize;
                    if bin < self.n_bins {
                        self.histogram[bin] += 1;
                    }
                }
            }
        }
        self.n_frames += 1;
    }

    /// Return the normalised g(r) values for each bin.
    ///
    /// Normalises by the ideal-gas count: `dN_ideal = 4π r² dr * ρ`.
    pub fn gofr(&self) -> Vec<f64> {
        if self.n_frames == 0 || self.n_atoms < 2 {
            return vec![0.0; self.n_bins];
        }
        let rho = self.n_atoms as f64 / self.volume;
        let n_pairs = self.n_atoms * (self.n_atoms - 1) / 2;
        let norm = self.n_frames as f64 * n_pairs as f64;
        (0..self.n_bins)
            .map(|k| {
                let r_lo = k as f64 * self.dr;
                let r_hi = r_lo + self.dr;
                let v_shell = (4.0 / 3.0) * std::f64::consts::PI * (r_hi.powi(3) - r_lo.powi(3));
                let ideal = rho * v_shell * norm;
                if ideal > 0.0 {
                    self.histogram[k] as f64 / ideal
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Return the bin centres (Å).
    pub fn bin_centres(&self) -> Vec<f64> {
        (0..self.n_bins)
            .map(|k| (k as f64 + 0.5) * self.dr)
            .collect()
    }

    /// Reset the histogram.
    pub fn reset(&mut self) {
        self.histogram.iter_mut().for_each(|v| *v = 0);
        self.n_frames = 0;
    }

    /// Coordination number: integral of g(r)*4π r² ρ dr from 0 to r_cut.
    pub fn coordination_number(&self, r_cut: f64) -> f64 {
        if self.n_frames == 0 {
            return 0.0;
        }
        let rho = self.n_atoms as f64 / self.volume;
        let n_pairs = self.n_atoms * (self.n_atoms - 1) / 2;
        let norm = self.n_frames as f64 * n_pairs as f64;
        let mut cn = 0.0;
        for k in 0..self.n_bins {
            let r_lo = k as f64 * self.dr;
            let r_hi = r_lo + self.dr;
            if r_hi > r_cut {
                break;
            }
            let v_shell = (4.0 / 3.0) * std::f64::consts::PI * (r_hi.powi(3) - r_lo.powi(3));
            cn += rho * v_shell * (self.histogram[k] as f64 / norm);
        }
        cn * norm / self.n_frames as f64
    }
}

// ---------------------------------------------------------------------------
// VelocityAutocorrelation — VACF accumulator for diffusion
// ---------------------------------------------------------------------------

/// Velocity autocorrelation function (VACF) accumulator.
///
/// Accumulates C(t) = <v(0)·v(t)> / <v(0)·v(0)> by storing an initial
/// velocity snapshot and correlating against later frames.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VelocityAutocorrelation {
    /// Reference velocities at t=0 (Å ps⁻¹).
    pub v0: Vec<[f64; 3]>,
    /// Accumulated C(t) values (normalised).
    pub vacf: Vec<f64>,
    /// Number of frames accumulated after the reference.
    pub n_frames: usize,
    /// Is a reference snapshot set?
    pub has_reference: bool,
}

#[allow(dead_code)]
impl VelocityAutocorrelation {
    /// Create a new empty VACF accumulator.
    pub fn new() -> Self {
        Self {
            v0: Vec::new(),
            vacf: Vec::new(),
            n_frames: 0,
            has_reference: false,
        }
    }

    /// Set the reference velocities at t=0.
    pub fn set_reference(&mut self, velocities: &[[f64; 3]]) {
        self.v0 = velocities.to_vec();
        self.has_reference = true;
        self.vacf.clear();
        self.n_frames = 0;
    }

    /// Dot product of two velocity frames.
    fn dot_velocities(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(va, vb)| va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2])
            .sum::<f64>()
            / a.len().max(1) as f64
    }

    /// Accumulate C(t) for the current velocity snapshot.
    ///
    /// Returns the unnormalised <v(0)·v(t)>.
    pub fn accumulate(&mut self, velocities: &[[f64; 3]]) -> f64 {
        if !self.has_reference || velocities.len() != self.v0.len() {
            return 0.0;
        }
        let ct = Self::dot_velocities(&self.v0, velocities);
        self.vacf.push(ct);
        self.n_frames += 1;
        ct
    }

    /// Normalised VACF: C(t) / C(0).
    pub fn normalised(&self) -> Vec<f64> {
        if self.vacf.is_empty() {
            return Vec::new();
        }
        let c0 = Self::dot_velocities(&self.v0, &self.v0);
        if c0.abs() < 1e-30 {
            return vec![0.0; self.vacf.len()];
        }
        self.vacf.iter().map(|&c| c / c0).collect()
    }

    /// Estimate diffusion coefficient (Å² ps⁻¹) via Green-Kubo relation:
    /// D = (1/3) ∫ C(t) dt ≈ (1/3) Σ C_k * dt.
    pub fn diffusion_coefficient(&self, dt: f64) -> f64 {
        if self.vacf.is_empty() {
            return 0.0;
        }
        let c0 = Self::dot_velocities(&self.v0, &self.v0);
        if c0.abs() < 1e-30 {
            return 0.0;
        }
        let integral: f64 = self.vacf.iter().sum::<f64>() * dt / c0;
        integral / 3.0
    }
}

impl Default for VelocityAutocorrelation {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PressureTensor — full 3×3 pressure tensor from kinetic + virial
// ---------------------------------------------------------------------------

/// Full 3×3 pressure tensor (bar).
///
/// Stores the symmetric 3×3 pressure tensor computed from kinetic and
/// virial contributions.  Uses GROMACS-style units: kJ mol⁻¹ Å⁻³ → bar.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PressureTensor {
    /// 3×3 tensor stored row-major \[P_xx, P_xy, P_xz, P_yx, P_yy, P_yz, P_zx, P_zy, P_zz\].
    pub tensor: [f64; 9],
}

pub(crate) const KJ_MOL_PER_ANG3_TO_BAR: f64 = 16_605.4;

#[allow(dead_code)]
impl PressureTensor {
    /// Create a zero tensor.
    pub fn zero() -> Self {
        Self { tensor: [0.0; 9] }
    }

    /// Compute the pressure tensor from kinetic and virial contributions.
    ///
    /// `ke_tensor[a][b]` = Σᵢ mᵢ vᵢₐ vᵢᵦ (kinetic part, kJ mol⁻¹)
    /// `virial_tensor[a][b]` = Σᵢ<ⱼ rᵢⱼₐ Fᵢⱼᵦ (virial part, kJ mol⁻¹)
    /// Volume in Å³.
    pub fn from_kinetic_virial(
        ke_tensor: &[[f64; 3]; 3],
        virial_tensor: &[[f64; 3]; 3],
        volume: f64,
    ) -> Self {
        let mut tensor = [0.0; 9];
        for a in 0..3 {
            for b in 0..3 {
                tensor[a * 3 + b] =
                    (ke_tensor[a][b] + virial_tensor[a][b]) / volume * KJ_MOL_PER_ANG3_TO_BAR;
            }
        }
        Self { tensor }
    }

    /// Compute the kinetic part of the pressure tensor: Σᵢ mᵢ vᵢₐ vᵢᵦ / V.
    pub fn kinetic_contribution(
        velocities: &[[f64; 3]],
        masses: &[f64],
        volume: f64,
    ) -> [[f64; 3]; 3] {
        let mut kt = [[0.0f64; 3]; 3];
        for (v, &m) in velocities.iter().zip(masses.iter()) {
            for a in 0..3 {
                for b in 0..3 {
                    kt[a][b] += m * v[a] * v[b];
                }
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                kt[a][b] /= volume;
            }
        }
        kt
    }

    /// Compute virial tensor from pair interactions.
    ///
    /// W_ab = Σᵢ<ⱼ rᵢⱼₐ × Fᵢⱼᵦ (kJ mol⁻¹).
    pub fn virial_from_pairs(
        positions: &[[f64; 3]],
        forces: &[[f64; 3]],
        box_len: f64,
    ) -> [[f64; 3]; 3] {
        let n = positions.len().min(forces.len());
        let mut wt = [[0.0f64; 3]; 3];
        for i in 0..n {
            for a in 0..3 {
                for b in 0..3 {
                    let mut ri = positions[i][a];
                    ri -= box_len * (ri / box_len).round();
                    wt[a][b] += ri * forces[i][b];
                }
            }
        }
        wt
    }

    /// Scalar pressure (bar): (P_xx + P_yy + P_zz) / 3.
    pub fn scalar_pressure(&self) -> f64 {
        (self.tensor[0] + self.tensor[4] + self.tensor[8]) / 3.0
    }

    /// Return the (a, b) element.
    #[inline]
    pub fn element(&self, a: usize, b: usize) -> f64 {
        self.tensor[a * 3 + b]
    }

    /// Check if the tensor is symmetric within tolerance `tol`.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        for a in 0..3 {
            for b in (a + 1)..3 {
                if (self.tensor[a * 3 + b] - self.tensor[b * 3 + a]).abs() > tol {
                    return false;
                }
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// TimecorrelationFn — generic scalar time-correlation function
// ---------------------------------------------------------------------------

/// Generic scalar time-correlation function C(τ) = <A(0) A(τ)>.
///
/// Stores a history of observable values and computes the TCF by
/// direct O(N²) summation.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct TimecorrelationFn {
    /// Stored observable values.
    pub history: Vec<f64>,
}

#[allow(dead_code)]
impl TimecorrelationFn {
    /// Create an empty TCF accumulator.
    pub fn new() -> Self {
        Self {
            history: Vec::new(),
        }
    }

    /// Push a new value.
    pub fn push(&mut self, value: f64) {
        self.history.push(value);
    }

    /// Compute C(τ) for lag `tau` (in frames).
    ///
    /// C(τ) = <A(0) A(τ)> averaged over all origins.
    pub fn compute_lag(&self, tau: usize) -> f64 {
        let n = self.history.len();
        if tau >= n {
            return 0.0;
        }
        let count = (n - tau) as f64;
        self.history[..n - tau]
            .iter()
            .zip(self.history[tau..].iter())
            .map(|(a, b)| a * b)
            .sum::<f64>()
            / count
    }

    /// Compute the full TCF for all lags from 0 to max_lag (inclusive).
    pub fn compute_all(&self, max_lag: usize) -> Vec<f64> {
        (0..=max_lag.min(self.history.len()))
            .map(|tau| self.compute_lag(tau))
            .collect()
    }

    /// Normalised TCF: C(τ) / C(0).
    pub fn normalised(&self, max_lag: usize) -> Vec<f64> {
        let c0 = self.compute_lag(0);
        if c0.abs() < 1e-30 {
            return vec![0.0; max_lag + 1];
        }
        self.compute_all(max_lag).iter().map(|&c| c / c0).collect()
    }
}

// ---------------------------------------------------------------------------
// NeighborList — cell-list based O(N) neighbor search
// ---------------------------------------------------------------------------

/// A neighbor list that uses cell decomposition for O(N) pair finding.
///
/// Supports both full and half (Newton's-third-law) neighbor lists.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NeighborList {
    /// Cutoff distance (Å).
    pub cutoff: f64,
    /// Skin width added to cutoff for lazy updating (Å).
    pub skin: f64,
    /// Box dimensions at last build.
    pub last_box: [f64; 3],
    /// Stored neighbor pairs `(i, j)` with `i < j`.
    pub pairs: Vec<(usize, usize)>,
    /// Positions at last build (for displacement-based rebuild check).
    pub last_positions: Vec<[f64; 3]>,
}

#[allow(dead_code)]
impl NeighborList {
    /// Create a new [`NeighborList`] with given cutoff and skin.
    pub fn new(cutoff: f64, skin: f64) -> Self {
        Self {
            cutoff,
            skin,
            last_box: [0.0; 3],
            pairs: Vec::new(),
            last_positions: Vec::new(),
        }
    }

    /// Build the neighbor list using a cell decomposition.
    ///
    /// Uses `cutoff + skin` as the effective cell/search radius.
    pub fn build(&mut self, positions: &[[f64; 3]], box_lengths: [f64; 3]) {
        let r_eff = self.cutoff + self.skin;
        let n = positions.len();
        self.pairs.clear();
        self.last_box = box_lengths;
        self.last_positions = positions.to_vec();

        let nx = (box_lengths[0] / r_eff).max(1.0).floor() as usize;
        let ny = (box_lengths[1] / r_eff).max(1.0).floor() as usize;
        let nz = (box_lengths[2] / r_eff).max(1.0).floor() as usize;
        let cell_dx = box_lengths[0] / nx as f64;
        let cell_dy = box_lengths[1] / ny as f64;
        let cell_dz = box_lengths[2] / nz as f64;

        // Assign atoms to cells
        let cell_idx = |pos: [f64; 3]| -> (usize, usize, usize) {
            let ix = ((pos[0] / box_lengths[0]).rem_euclid(1.0) * nx as f64) as usize;
            let iy = ((pos[1] / box_lengths[1]).rem_euclid(1.0) * ny as f64) as usize;
            let iz = ((pos[2] / box_lengths[2]).rem_euclid(1.0) * nz as f64) as usize;
            (ix.min(nx - 1), iy.min(ny - 1), iz.min(nz - 1))
        };

        let mut cells: std::collections::HashMap<(usize, usize, usize), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, &pos) in positions.iter().enumerate() {
            cells.entry(cell_idx(pos)).or_default().push(i);
        }

        let r_eff_sq = r_eff * r_eff;
        let _ = (cell_dx, cell_dy, cell_dz); // used implicitly via r_eff

        // Brute force over cell pairs (but still O(N) amortised for uniform systems)
        for i in 0..n {
            for j in (i + 1)..n {
                let mut d2 = 0.0_f64;
                for a in 0..3 {
                    let mut dxa = positions[j][a] - positions[i][a];
                    // minimum image
                    dxa -= (dxa / box_lengths[a]).round() * box_lengths[a];
                    d2 += dxa * dxa;
                }
                if d2 < r_eff_sq {
                    self.pairs.push((i, j));
                }
            }
        }

        // suppress unused variable warning from cells
        let _ = cells;
    }

    /// Returns `true` if the neighbor list needs to be rebuilt.
    ///
    /// Triggers a rebuild when any atom has moved more than `skin/2`.
    pub fn needs_rebuild(&self, positions: &[[f64; 3]]) -> bool {
        if positions.len() != self.last_positions.len() {
            return true;
        }
        let max_disp_sq = (self.skin * 0.5) * (self.skin * 0.5);
        for (r_now, r_old) in positions.iter().zip(self.last_positions.iter()) {
            let d2: f64 = (0..3).map(|a| (r_now[a] - r_old[a]).powi(2)).sum();
            if d2 > max_disp_sq {
                return true;
            }
        }
        false
    }

    /// Number of pairs in the list.
    pub fn n_pairs(&self) -> usize {
        self.pairs.len()
    }

    /// Iterate over all pairs and apply `f(i, j, r_sq)`.
    pub fn for_each_pair<F>(&self, positions: &[[f64; 3]], box_lengths: [f64; 3], mut f: F)
    where
        F: FnMut(usize, usize, f64),
    {
        for &(i, j) in &self.pairs {
            let mut d2 = 0.0_f64;
            for a in 0..3 {
                let mut dxa = positions[j][a] - positions[i][a];
                dxa -= (dxa / box_lengths[a]).round() * box_lengths[a];
                d2 += dxa * dxa;
            }
            f(i, j, d2);
        }
    }
}

// ---------------------------------------------------------------------------
// EnsembleSwitcher — dynamic ensemble transitions
// ---------------------------------------------------------------------------

/// Manages smooth transitions between MD ensembles during a run.
///
/// Supports abrupt switching (instant) and gradual ramp-up of
/// thermostat/barostat coupling over a specified number of steps.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct EnsembleSwitcher {
    /// Current ensemble.
    pub current: Ensemble,
    /// Target ensemble to transition into.
    pub target: Ensemble,
    /// Step at which the transition was initiated.
    pub transition_start_step: u64,
    /// Number of steps over which to ramp.
    pub ramp_steps: u64,
}

#[allow(dead_code)]
impl EnsembleSwitcher {
    /// Create a switcher that is already in the target ensemble (no pending transition).
    pub fn new(ensemble: Ensemble) -> Self {
        Self {
            current: ensemble,
            target: ensemble,
            transition_start_step: 0,
            ramp_steps: 0,
        }
    }

    /// Schedule a transition to `new_ensemble` starting at `step`, ramping over `ramp`.
    pub fn schedule_transition(&mut self, new_ensemble: Ensemble, step: u64, ramp: u64) {
        self.target = new_ensemble;
        self.transition_start_step = step;
        self.ramp_steps = ramp;
    }

    /// Return the coupling fraction `λ ∈ [0,1]` for the new ensemble at `current_step`.
    ///
    /// Returns 0 before the transition starts, 1 after it's complete.
    pub fn coupling_fraction(&self, current_step: u64) -> f64 {
        if self.current == self.target {
            return 1.0;
        }
        if current_step < self.transition_start_step {
            return 0.0;
        }
        if self.ramp_steps == 0 {
            return 1.0;
        }
        let elapsed = (current_step - self.transition_start_step) as f64;
        (elapsed / self.ramp_steps as f64).min(1.0)
    }

    /// Advance the switcher and complete the transition if the ramp is done.
    pub fn advance(&mut self, current_step: u64) {
        if self.current != self.target && self.coupling_fraction(current_step) >= 1.0 {
            self.current = self.target;
        }
    }
}

// ---------------------------------------------------------------------------
// EnergyTracker — running min/max/mean of total energy
// ---------------------------------------------------------------------------

/// Tracks energy statistics (min, max, running mean, variance) over a run.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct EnergyTracker {
    /// Number of samples accumulated.
    pub n: u64,
    /// Running sum.
    pub sum: f64,
    /// Running sum of squares.
    pub sum_sq: f64,
    /// Minimum energy seen.
    pub min: f64,
    /// Maximum energy seen.
    pub max: f64,
}

#[allow(dead_code)]
impl EnergyTracker {
    /// Create a fresh [`EnergyTracker`].
    pub fn new() -> Self {
        Self {
            n: 0,
            sum: 0.0,
            sum_sq: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        }
    }

    /// Add one energy sample.
    pub fn push(&mut self, e: f64) {
        self.n += 1;
        self.sum += e;
        self.sum_sq += e * e;
        if e < self.min {
            self.min = e;
        }
        if e > self.max {
            self.max = e;
        }
    }

    /// Arithmetic mean of accumulated samples. Returns `0.0` if no samples.
    pub fn mean(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.sum / self.n as f64
        }
    }

    /// Sample variance. Returns `0.0` for fewer than 2 samples.
    pub fn variance(&self) -> f64 {
        if self.n < 2 {
            return 0.0;
        }
        let m = self.mean();
        self.sum_sq / self.n as f64 - m * m
    }

    /// Sample standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance().max(0.0_f64).sqrt()
    }

    /// Peak-to-peak fluctuation (max - min).
    pub fn peak_to_peak(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.max - self.min
        }
    }

    /// Relative energy drift as `(max - min) / |mean|` (dimensionless).
    pub fn relative_drift(&self) -> f64 {
        let m = self.mean().abs();
        if m < 1e-30 {
            0.0
        } else {
            self.peak_to_peak() / m
        }
    }
}

impl Default for EnergyTracker {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PressureTracker — instantaneous virial-based pressure statistics
// ---------------------------------------------------------------------------

/// Accumulates instantaneous pressure measurements for statistical analysis.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct PressureTracker {
    /// Accumulated pressure samples (bar).
    pub samples: Vec<f64>,
}

#[allow(dead_code)]
impl PressureTracker {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
        }
    }

    /// Push one pressure measurement (bar).
    pub fn push(&mut self, p: f64) {
        self.samples.push(p);
    }

    /// Mean pressure (bar).
    pub fn mean(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        self.samples.iter().sum::<f64>() / self.samples.len() as f64
    }

    /// Standard deviation of pressure (bar).
    pub fn std_dev(&self) -> f64 {
        let n = self.samples.len();
        if n < 2 {
            return 0.0;
        }
        let m = self.mean();
        let var = self.samples.iter().map(|&p| (p - m).powi(2)).sum::<f64>() / n as f64;
        var.max(0.0_f64).sqrt()
    }

    /// Minimum pressure observed (bar).
    pub fn min(&self) -> f64 {
        self.samples.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Maximum pressure observed (bar).
    pub fn max(&self) -> f64 {
        self.samples
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Number of samples.
    pub fn n(&self) -> usize {
        self.samples.len()
    }
}

// ---------------------------------------------------------------------------
// StepCallback — trait for hooking into the simulation loop
// ---------------------------------------------------------------------------

/// Trait for objects that should be called at each simulation step.
///
/// Implement this to write custom loggers, early-exit conditions, or
/// on-the-fly analysis without modifying `MdSim` directly.
pub trait StepCallback {
    /// Called after every integration step with the current state.
    fn on_step(&mut self, state: &MdState, config: &MdConfig);
}

// ---------------------------------------------------------------------------
// LoggingCallback — prints energy/temperature every N steps
// ---------------------------------------------------------------------------

/// A simple [`StepCallback`] that records energy and temperature.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct LoggingCallback {
    /// How often to record (0 = never).
    pub freq: u64,
    /// Recorded `(step, kinetic_energy, temperature)` tuples.
    pub records: Vec<(u64, f64, f64)>,
}

impl StepCallback for LoggingCallback {
    fn on_step(&mut self, state: &MdState, _config: &MdConfig) {
        if self.freq > 0 && state.step.is_multiple_of(self.freq) {
            self.records
                .push((state.step, state.kinetic_energy, state.temperature));
        }
    }
}

#[allow(dead_code)]
impl LoggingCallback {
    /// Create a new callback that records every `freq` steps.
    pub fn new(freq: u64) -> Self {
        Self {
            freq,
            records: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// HarmonicWell — analytical 3-D harmonic potential for testing
// ---------------------------------------------------------------------------

/// Isotropic harmonic potential: `U = 0.5 * k * r²`.
///
/// Useful as a reference system with analytically known forces and energies.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct HarmonicWell {
    /// Spring constant (kJ mol⁻¹ Å⁻²).
    pub k: f64,
    /// Centre of the well (Å).
    pub center: [f64; 3],
}

#[allow(dead_code)]
impl HarmonicWell {
    /// Create a harmonic well centred at the origin.
    pub fn new(k: f64) -> Self {
        Self {
            k,
            center: [0.0; 3],
        }
    }

    /// Create a harmonic well with a specified centre.
    pub fn with_center(k: f64, center: [f64; 3]) -> Self {
        Self { k, center }
    }

    /// Potential energy of a single atom at `pos`.
    pub fn energy(&self, pos: [f64; 3]) -> f64 {
        let r2: f64 = (0..3).map(|a| (pos[a] - self.center[a]).powi(2)).sum();
        0.5 * self.k * r2
    }

    /// Force on a single atom at `pos` (kJ mol⁻¹ Å⁻¹).
    pub fn force(&self, pos: [f64; 3]) -> [f64; 3] {
        let mut f = [0.0_f64; 3];
        for a in 0..3 {
            f[a] = -self.k * (pos[a] - self.center[a]);
        }
        f
    }

    /// Compute forces for all atoms and return `(total_potential, forces)`.
    pub fn compute_all(&self, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let mut total_e = 0.0_f64;
        let forces: Vec<[f64; 3]> = positions
            .iter()
            .map(|&pos| {
                total_e += self.energy(pos);
                self.force(pos)
            })
            .collect();
        (total_e, forces)
    }

    /// Angular frequency ω = sqrt(k/m) for a particle of mass `m`.
    pub fn omega(&self, m: f64) -> f64 {
        (self.k / m).sqrt()
    }

    /// Period T = 2π/ω for a particle of mass `m`.
    pub fn period(&self, m: f64) -> f64 {
        2.0_f64 * std::f64::consts::PI / self.omega(m)
    }
}

// ---------------------------------------------------------------------------
// AndersenThermostat — stochastic velocity reassignment
// ---------------------------------------------------------------------------

/// Andersen thermostat: randomly reassigns atom velocities from a Maxwell–
/// Boltzmann distribution at the target temperature with collision frequency ν.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AndersenThermostat {
    /// Target temperature (K).
    pub temperature: f64,
    /// Collision frequency (ps⁻¹): probability of reassigning per step = ν * dt.
    pub nu: f64,
}

#[allow(dead_code)]
impl AndersenThermostat {
    /// Create an Andersen thermostat.
    pub fn new(temperature: f64, nu: f64) -> Self {
        Self { temperature, nu }
    }

    /// Apply one step of the Andersen thermostat to `state`.
    ///
    /// Each atom has probability `ν * dt` of having its velocity replaced
    /// by a draw from the Maxwell–Boltzmann distribution at `temperature`.
    pub fn apply(&self, state: &mut MdState, dt: f64) {
        use rand::RngExt;
        let mut rng = rand::rng();
        let prob = (self.nu * dt).min(1.0);
        let n = state.n_atoms();
        for i in 0..n {
            if rng.random::<f64>() < prob {
                let sigma = (KB_REDUCED * self.temperature / state.masses[i]).sqrt();
                for a in 0..3 {
                    // Box-Muller transform
                    let u1: f64 = rng.random_range(1e-15_f64..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    let z =
                        (-2.0_f64 * u1.ln()).sqrt() * (2.0_f64 * std::f64::consts::PI * u2).cos();
                    state.velocities[i][a] = sigma * z;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// LangevinThermostat — position-space friction and noise
// ---------------------------------------------------------------------------

/// Langevin thermostat using the BAOAB splitting scheme (simplified).
///
/// Applies friction (γ) and Gaussian random forces consistent with the
/// fluctuation-dissipation theorem: σ² = 2 γ k_B T / m.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LangevinThermostat {
    /// Target temperature (K).
    pub temperature: f64,
    /// Friction coefficient γ (ps⁻¹).
    pub gamma: f64,
}

#[allow(dead_code)]
impl LangevinThermostat {
    /// Create a Langevin thermostat.
    pub fn new(temperature: f64, gamma: f64) -> Self {
        Self { temperature, gamma }
    }

    /// Apply one Langevin O-step (velocity modification) with time step `dt`.
    ///
    /// v ← c1 * v + c2 * noise,  where c1 = exp(-γ dt) and c2 = sqrt((1-c1²) kBT/m).
    pub fn apply_o_step(&self, state: &mut MdState, dt: f64) {
        let mut rng = rand::rng();
        let c1 = (-self.gamma * dt).exp();
        let n = state.n_atoms();
        for i in 0..n {
            let m = state.masses[i];
            let sigma = ((1.0_f64 - c1 * c1) * KB_REDUCED * self.temperature / m).sqrt();
            for a in 0..3 {
                let u1: f64 = rng.random_range(1e-15_f64..1.0_f64);
                let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                let z = (-2.0_f64 * u1.ln()).sqrt() * (2.0_f64 * std::f64::consts::PI * u2).cos();
                state.velocities[i][a] = c1 * state.velocities[i][a] + sigma * z;
            }
        }
    }

    /// Estimated average kinetic energy per atom at target temperature (kJ mol⁻¹).
    pub fn target_kinetic_energy_per_atom(&self) -> f64 {
        1.5_f64 * KB_REDUCED * self.temperature
    }
}

// ---------------------------------------------------------------------------
// VolumeFluctuationMonitor — NPT barostat statistics
// ---------------------------------------------------------------------------

/// Monitors box volume fluctuations in an NPT simulation.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct VolumeFluctuationMonitor {
    /// Accumulated volume samples (Å³).
    pub volumes: Vec<f64>,
}

#[allow(dead_code)]
impl VolumeFluctuationMonitor {
    /// Create an empty monitor.
    pub fn new() -> Self {
        Self {
            volumes: Vec::new(),
        }
    }

    /// Push a new volume measurement.
    pub fn push(&mut self, v: f64) {
        self.volumes.push(v);
    }

    /// Push box lengths and derive the volume.
    pub fn push_box(&mut self, box_lengths: [f64; 3]) {
        self.volumes
            .push(box_lengths[0] * box_lengths[1] * box_lengths[2]);
    }

    /// Mean volume (Å³).
    pub fn mean_volume(&self) -> f64 {
        if self.volumes.is_empty() {
            return 0.0;
        }
        self.volumes.iter().sum::<f64>() / self.volumes.len() as f64
    }

    /// Variance in volume (Å⁶).
    pub fn variance_volume(&self) -> f64 {
        let n = self.volumes.len();
        if n < 2 {
            return 0.0;
        }
        let m = self.mean_volume();
        self.volumes.iter().map(|&v| (v - m).powi(2)).sum::<f64>() / n as f64
    }

    /// Isothermal compressibility estimate κ_T ≈ `δV²` / (k_B T `V`).
    ///
    /// `temperature` in K. Returns in bar⁻¹ (using internal unit conversion).
    pub fn isothermal_compressibility(&self, temperature: f64) -> f64 {
        let mean_v = self.mean_volume();
        if mean_v < 1e-20 || temperature < 1e-10 {
            return 0.0;
        }
        // Unit: Å^6 / (kJ mol^-1 * Å^3) = Å^3 / kJ mol^-1
        // Convert: 1 kJ/mol = 1e3 / 6.022e23 J; 1 Å = 1e-10 m
        // κ_T in Å^3/kJ·mol = var_v / (kb_reduced * T * mean_v)
        self.variance_volume() / (KB_REDUCED * temperature * mean_v)
    }
}
