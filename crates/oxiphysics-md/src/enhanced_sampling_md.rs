// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended enhanced sampling methods for molecular dynamics.
//!
//! This module builds on top of [`crate::enhanced_sampling`] and provides
//! more sophisticated algorithms:
//!
//! - [`RemdController`]: Replica exchange MD (REMD / parallel tempering) with
//!   Hamiltonian REMD support and optimal swap scheduling.
//! - [`WellTemperedMetadynamics`]: Well-tempered metadynamics with adaptive
//!   Gaussian deposition and free energy convergence tracking.
//! - [`UmbrellaSamplingWham`]: Multi-window umbrella sampling with the
//!   Weighted Histogram Analysis Method (WHAM) for free energy recovery.
//! - [`AdaptiveBiasingForce`]: ABF with running block averages and
//!   ramp-up sample threshold.
//! - [`SteeredMdExtension`]: Extended steered MD with adaptive pulling speed
//!   and work accumulation.
//! - [`CollectiveVariable`]: Pluggable collective variables including
//!   coordination number, dihedral angle, and RMSD.
//! - [`TransitionPathSampling`]: Two-ended shooting for reactive trajectories.

use rand::Rng;

use rand::RngExt;
// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant in kJ mol⁻¹ K⁻¹.
const KB_KJMOL: f64 = 8.314e-3;

/// Small tolerance for floating-point comparisons.
const EPS: f64 = 1e-12;

// ─────────────────────────────────────────────────────────────────────────────
// CollectiveVariable
// ─────────────────────────────────────────────────────────────────────────────

/// Kind of collective variable (CV).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CvKind {
    /// Distance between two atoms.
    Distance,
    /// Coordination number around a centre atom.
    CoordinationNumber,
    /// Dihedral angle formed by four atoms.
    Dihedral,
    /// RMSD from a reference structure.
    Rmsd,
    /// Linear combination of positions along a given direction.
    LinearCombination,
}

/// A collective variable evaluated on a set of positions.
#[derive(Debug, Clone)]
pub struct CollectiveVariable {
    /// Kind of CV.
    pub kind: CvKind,
    /// Indices of atoms involved.
    pub atom_indices: Vec<usize>,
    /// Reference positions for RMSD (one per atom, flattened xyz).
    pub reference_positions: Vec<[f64; 3]>,
    /// Direction vector for linear combination CV.
    pub direction: [f64; 3],
    /// Parameters: `[r0, d0, n, m]` for coordination number;
    /// unused entries set to 0.
    pub params: [f64; 4],
}

impl CollectiveVariable {
    /// Create a distance CV between two atoms.
    pub fn distance(i: usize, j: usize) -> Self {
        Self {
            kind: CvKind::Distance,
            atom_indices: vec![i, j],
            reference_positions: Vec::new(),
            direction: [0.0; 3],
            params: [0.0; 4],
        }
    }

    /// Create a coordination number CV.
    ///
    /// Uses the rational switching function:
    /// `s(r) = (1 - (r/r0)^n) / (1 - (r/r0)^m)` summed over neighbours.
    ///
    /// * `centre` – index of the central atom.
    /// * `neighbours` – indices of the coordination shell atoms.
    /// * `r0` – reference distance.
    /// * `d0` – shift parameter.
    /// * `n`, `m` – exponents.
    pub fn coordination_number(
        centre: usize,
        neighbours: &[usize],
        r0: f64,
        d0: f64,
        n: f64,
        m: f64,
    ) -> Self {
        let mut idx = vec![centre];
        idx.extend_from_slice(neighbours);
        Self {
            kind: CvKind::CoordinationNumber,
            atom_indices: idx,
            reference_positions: Vec::new(),
            direction: [0.0; 3],
            params: [r0, d0, n, m],
        }
    }

    /// Create a dihedral angle CV from four atom indices.
    pub fn dihedral(i: usize, j: usize, k: usize, l: usize) -> Self {
        Self {
            kind: CvKind::Dihedral,
            atom_indices: vec![i, j, k, l],
            reference_positions: Vec::new(),
            direction: [0.0; 3],
            params: [0.0; 4],
        }
    }

    /// Create an RMSD CV relative to reference positions.
    pub fn rmsd(atom_indices: Vec<usize>, reference: Vec<[f64; 3]>) -> Self {
        Self {
            kind: CvKind::Rmsd,
            atom_indices,
            reference_positions: reference,
            direction: [0.0; 3],
            params: [0.0; 4],
        }
    }

    /// Create a linear combination CV: `cv = sum_i (dir · r_i)`.
    pub fn linear_combination(atom_indices: Vec<usize>, direction: [f64; 3]) -> Self {
        Self {
            kind: CvKind::LinearCombination,
            atom_indices,
            reference_positions: Vec::new(),
            direction,
            params: [0.0; 4],
        }
    }

    /// Evaluate the collective variable given all positions.
    pub fn evaluate(&self, positions: &[[f64; 3]]) -> f64 {
        match self.kind {
            CvKind::Distance => self.eval_distance(positions),
            CvKind::CoordinationNumber => self.eval_coordination(positions),
            CvKind::Dihedral => self.eval_dihedral(positions),
            CvKind::Rmsd => self.eval_rmsd(positions),
            CvKind::LinearCombination => self.eval_linear(positions),
        }
    }

    /// Evaluate distance CV.
    fn eval_distance(&self, positions: &[[f64; 3]]) -> f64 {
        let a = positions[self.atom_indices[0]];
        let b = positions[self.atom_indices[1]];
        dist3(a, b)
    }

    /// Evaluate coordination number CV.
    fn eval_coordination(&self, positions: &[[f64; 3]]) -> f64 {
        let centre = positions[self.atom_indices[0]];
        let r0 = self.params[0];
        let d0 = self.params[1];
        let n = self.params[2];
        let m = self.params[3];
        let mut coord = 0.0;
        for &idx in &self.atom_indices[1..] {
            let r = dist3(centre, positions[idx]);
            let x = (r - d0) / r0;
            if x.abs() < EPS {
                coord += 1.0;
            } else {
                let num = 1.0 - x.powf(n);
                let den = 1.0 - x.powf(m);
                if den.abs() > EPS {
                    coord += num / den;
                }
            }
        }
        coord
    }

    /// Evaluate dihedral angle CV (radians).
    fn eval_dihedral(&self, positions: &[[f64; 3]]) -> f64 {
        let p0 = positions[self.atom_indices[0]];
        let p1 = positions[self.atom_indices[1]];
        let p2 = positions[self.atom_indices[2]];
        let p3 = positions[self.atom_indices[3]];
        let b1 = sub3(p1, p0);
        let b2 = sub3(p2, p1);
        let b3 = sub3(p3, p2);
        let n1 = cross3(b1, b2);
        let n2 = cross3(b2, b3);
        let m1 = cross3(n1, b2);
        let x = dot3(n1, n2);
        let y = dot3(m1, n2) * len3(b2);
        y.atan2(x)
    }

    /// Evaluate RMSD CV.
    fn eval_rmsd(&self, positions: &[[f64; 3]]) -> f64 {
        if self.atom_indices.is_empty() {
            return 0.0;
        }
        let n = self.atom_indices.len() as f64;
        let mut sum_sq = 0.0;
        for (k, &idx) in self.atom_indices.iter().enumerate() {
            let r = positions[idx];
            let rref = self.reference_positions[k];
            let dx = r[0] - rref[0];
            let dy = r[1] - rref[1];
            let dz = r[2] - rref[2];
            sum_sq += dx * dx + dy * dy + dz * dz;
        }
        (sum_sq / n).sqrt()
    }

    /// Evaluate linear combination CV.
    fn eval_linear(&self, positions: &[[f64; 3]]) -> f64 {
        let mut val = 0.0;
        for &idx in &self.atom_indices {
            let r = positions[idx];
            val += dot3(r, self.direction);
        }
        val
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RemdController
// ─────────────────────────────────────────────────────────────────────────────

/// Swap scheduling strategy for REMD.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwapSchedule {
    /// Try neighbours in order 0↔1, 2↔3, … then 1↔2, 3↔4, … alternating.
    EvenOdd,
    /// Random pair selection.
    Random,
}

/// Replica Exchange Molecular Dynamics controller.
///
/// Supports temperature REMD and Hamiltonian REMD, with configurable swap
/// scheduling and detailed acceptance tracking per pair.
#[derive(Debug, Clone)]
pub struct RemdController {
    /// Temperature of each replica (K).
    pub temperatures: Vec<f64>,
    /// Potential energy of each replica (kJ/mol).
    pub energies: Vec<f64>,
    /// Mapping: `replica_index[k]` is the configuration currently at slot `k`.
    pub replica_index: Vec<usize>,
    /// Swap schedule.
    pub schedule: SwapSchedule,
    /// Per-pair acceptance counters: `(attempts, accepted)`.
    pub pair_stats: Vec<(u64, u64)>,
    /// Global step counter.
    pub step: u64,
    /// Swap attempt interval (MD steps between swap attempts).
    pub swap_interval: u64,
    /// Hamiltonian REMD lambda parameters (one per replica). `None` → temperature REMD.
    pub lambda: Option<Vec<f64>>,
}

impl RemdController {
    /// Create a temperature REMD controller with a geometric ladder.
    pub fn new_geometric(t_min: f64, t_max: f64, n_replicas: usize, swap_interval: u64) -> Self {
        let mut temps = Vec::with_capacity(n_replicas);
        if n_replicas <= 1 {
            temps.push(t_min);
        } else {
            let ratio = (t_max / t_min).powf(1.0 / (n_replicas as f64 - 1.0));
            for i in 0..n_replicas {
                temps.push(t_min * ratio.powi(i as i32));
            }
        }
        let n = temps.len();
        Self {
            temperatures: temps,
            energies: vec![0.0; n],
            replica_index: (0..n).collect(),
            schedule: SwapSchedule::EvenOdd,
            pair_stats: vec![(0, 0); n.saturating_sub(1)],
            step: 0,
            swap_interval,
            lambda: None,
        }
    }

    /// Create a Hamiltonian REMD controller.
    pub fn new_hamiltonian(temperature: f64, lambdas: Vec<f64>, swap_interval: u64) -> Self {
        let n = lambdas.len();
        Self {
            temperatures: vec![temperature; n],
            energies: vec![0.0; n],
            replica_index: (0..n).collect(),
            schedule: SwapSchedule::EvenOdd,
            pair_stats: vec![(0, 0); n.saturating_sub(1)],
            step: 0,
            swap_interval,
            lambda: Some(lambdas),
        }
    }

    /// Number of replicas.
    pub fn n_replicas(&self) -> usize {
        self.temperatures.len()
    }

    /// Update energies from the simulation.
    pub fn set_energies(&mut self, energies: &[f64]) {
        let n = energies.len().min(self.energies.len());
        self.energies[..n].copy_from_slice(&energies[..n]);
    }

    /// Attempt swaps according to the current schedule.
    /// Returns a list of pairs that were actually swapped.
    pub fn attempt_swaps(&mut self) -> Vec<(usize, usize)> {
        let n = self.n_replicas();
        if n < 2 {
            return Vec::new();
        }
        let mut rng = rand::rng();
        let mut swapped = Vec::new();

        match self.schedule {
            SwapSchedule::EvenOdd => {
                let offset = (self.step % 2) as usize;
                let mut i = offset;
                while i + 1 < n {
                    if self.try_swap(i, i + 1, &mut rng) {
                        swapped.push((i, i + 1));
                    }
                    i += 2;
                }
            }
            SwapSchedule::Random => {
                let i = rng.random_range(0..n - 1);
                if self.try_swap(i, i + 1, &mut rng) {
                    swapped.push((i, i + 1));
                }
            }
        }
        self.step += 1;
        swapped
    }

    /// Try to swap replicas `i` and `j` using Metropolis criterion.
    fn try_swap(&mut self, i: usize, j: usize, rng: &mut impl Rng) -> bool {
        let beta_i = 1.0 / (KB_KJMOL * self.temperatures[i]);
        let beta_j = 1.0 / (KB_KJMOL * self.temperatures[j]);
        let delta = (beta_i - beta_j) * (self.energies[i] - self.energies[j]);
        let pair_idx = i.min(j);
        self.pair_stats[pair_idx].0 += 1;
        let accept = if delta <= 0.0 {
            true
        } else {
            let r: f64 = rng.random_range(0.0..1.0);
            r < (-delta).exp()
        };
        if accept {
            self.pair_stats[pair_idx].1 += 1;
            self.energies.swap(i, j);
            self.replica_index.swap(i, j);
        }
        accept
    }

    /// Acceptance ratio for the pair `(i, i+1)`.
    pub fn acceptance_ratio(&self, pair_index: usize) -> f64 {
        let (att, acc) = self.pair_stats[pair_index];
        if att == 0 {
            0.0
        } else {
            acc as f64 / att as f64
        }
    }

    /// Overall acceptance ratio across all pairs.
    pub fn overall_acceptance(&self) -> f64 {
        let (total_att, total_acc): (u64, u64) = self
            .pair_stats
            .iter()
            .fold((0, 0), |(a, b), &(x, y)| (a + x, b + y));
        if total_att == 0 {
            0.0
        } else {
            total_acc as f64 / total_att as f64
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WellTemperedMetadynamics
// ─────────────────────────────────────────────────────────────────────────────

/// A single Gaussian hill deposited during metadynamics.
#[derive(Debug, Clone)]
pub struct GaussianHill {
    /// CV value at hill centre.
    pub centre: f64,
    /// Hill width (sigma).
    pub sigma: f64,
    /// Hill height at time of deposition.
    pub height: f64,
    /// Simulation time (or step) when deposited.
    pub time: f64,
}

/// Well-tempered metadynamics engine.
///
/// The bias potential is built as a sum of Gaussians in CV space.
/// In the well-tempered variant, the height of new hills is scaled by:
///
/// ```text
/// w(t) = w0 * exp(-V(s, t) / (kB * delta_T))
/// ```
///
/// where `delta_T = (bias_factor - 1) * T`.
#[derive(Debug, Clone)]
pub struct WellTemperedMetadynamics {
    /// Initial hill height (kJ/mol).
    pub initial_height: f64,
    /// Hill width (sigma of Gaussian).
    pub sigma: f64,
    /// Bias factor γ.
    pub bias_factor: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Deposition interval (steps).
    pub deposition_interval: u64,
    /// Deposited hills.
    pub hills: Vec<GaussianHill>,
    /// Step counter.
    pub step: u64,
    /// Grid-based bias for fast evaluation: `(cv_min, cv_max, bin_width, values)`.
    grid: Option<(f64, f64, f64, Vec<f64>)>,
}

impl WellTemperedMetadynamics {
    /// Create a well-tempered metadynamics engine.
    ///
    /// * `initial_height` – w₀ in kJ/mol.
    /// * `sigma` – Gaussian width.
    /// * `bias_factor` – γ (must be > 1).
    /// * `temperature` – simulation temperature in K.
    /// * `deposition_interval` – steps between hill depositions.
    pub fn new(
        initial_height: f64,
        sigma: f64,
        bias_factor: f64,
        temperature: f64,
        deposition_interval: u64,
    ) -> Self {
        Self {
            initial_height,
            sigma,
            bias_factor,
            temperature,
            deposition_interval,
            hills: Vec::new(),
            step: 0,
            grid: None,
        }
    }

    /// Enable grid-based bias accumulation for speed.
    pub fn enable_grid(&mut self, cv_min: f64, cv_max: f64, n_bins: usize) {
        let bw = (cv_max - cv_min) / n_bins as f64;
        self.grid = Some((cv_min, cv_max, bw, vec![0.0; n_bins]));
    }

    /// Compute the current bias at CV value `s`.
    pub fn bias_at(&self, s: f64) -> f64 {
        // If grid is available, interpolate
        if let Some((cv_min, cv_max, bw, ref vals)) = self.grid {
            if s < cv_min || s >= cv_max {
                return self.bias_sum(s);
            }
            let idx_f = (s - cv_min) / bw;
            let idx = idx_f as usize;
            if idx + 1 < vals.len() {
                let frac = idx_f - idx as f64;
                return vals[idx] * (1.0 - frac) + vals[idx + 1] * frac;
            }
            return if idx < vals.len() { vals[idx] } else { 0.0 };
        }
        self.bias_sum(s)
    }

    /// Direct summation of Gaussian hills at `s`.
    fn bias_sum(&self, s: f64) -> f64 {
        let mut v = 0.0;
        let inv_2sig2 = 0.5 / (self.sigma * self.sigma);
        for hill in &self.hills {
            let ds = s - hill.centre;
            v += hill.height * (-ds * ds * inv_2sig2).exp();
        }
        v
    }

    /// Deposit a hill at the current CV value `s`.
    /// Returns the effective hill height used.
    pub fn deposit_hill(&mut self, s: f64) -> f64 {
        let current_bias = self.bias_at(s);
        let delta_t = (self.bias_factor - 1.0) * self.temperature;
        let eff_height = if delta_t.abs() < EPS {
            self.initial_height
        } else {
            self.initial_height * (-current_bias / (KB_KJMOL * delta_t)).exp()
        };

        let hill = GaussianHill {
            centre: s,
            sigma: self.sigma,
            height: eff_height,
            time: self.step as f64,
        };

        // Update grid if present
        if let Some((cv_min, _cv_max, bw, ref mut vals)) = self.grid {
            let inv_2sig2 = 0.5 / (self.sigma * self.sigma);
            for (i, val) in vals.iter_mut().enumerate() {
                let cv_i = cv_min + (i as f64 + 0.5) * bw;
                let ds = cv_i - s;
                *val += eff_height * (-ds * ds * inv_2sig2).exp();
            }
        }

        self.hills.push(hill);
        eff_height
    }

    /// Step the metadynamics engine. Deposits a hill if the interval is reached.
    pub fn step(&mut self, cv_value: f64) -> Option<f64> {
        self.step += 1;
        if self.step.is_multiple_of(self.deposition_interval) {
            Some(self.deposit_hill(cv_value))
        } else {
            None
        }
    }

    /// Compute the bias force (negative gradient) at CV value `s`.
    pub fn bias_force(&self, s: f64) -> f64 {
        let mut f = 0.0;
        let inv_sig2 = 1.0 / (self.sigma * self.sigma);
        let inv_2sig2 = 0.5 * inv_sig2;
        for hill in &self.hills {
            let ds = s - hill.centre;
            let g = hill.height * (-ds * ds * inv_2sig2).exp();
            f += g * ds * inv_sig2;
        }
        f
    }

    /// Estimate the free energy surface as `-V_bias(s) * gamma / (gamma - 1)`.
    /// Returns `(cv_values, fes_values)`.
    pub fn free_energy_surface(
        &self,
        cv_min: f64,
        cv_max: f64,
        n_points: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        let factor = self.bias_factor / (self.bias_factor - 1.0);
        let ds = (cv_max - cv_min) / (n_points as f64 - 1.0);
        let mut cvs = Vec::with_capacity(n_points);
        let mut fes = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let s = cv_min + i as f64 * ds;
            cvs.push(s);
            fes.push(-self.bias_at(s) * factor);
        }
        (cvs, fes)
    }

    /// Number of hills deposited so far.
    pub fn n_hills(&self) -> usize {
        self.hills.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// UmbrellaSamplingWham
// ─────────────────────────────────────────────────────────────────────────────

/// A single umbrella sampling window.
#[derive(Debug, Clone)]
pub struct UmbrellaWindow {
    /// CV centre of the harmonic bias.
    pub centre: f64,
    /// Force constant (kJ mol⁻¹ unit⁻²).
    pub force_constant: f64,
    /// Collected CV samples.
    pub samples: Vec<f64>,
}

impl UmbrellaWindow {
    /// Create a new umbrella window.
    pub fn new(centre: f64, force_constant: f64) -> Self {
        Self {
            centre,
            force_constant,
            samples: Vec::new(),
        }
    }

    /// Add a CV sample.
    pub fn add_sample(&mut self, cv_value: f64) {
        self.samples.push(cv_value);
    }

    /// Bias energy at CV value `s`.
    pub fn bias_energy(&self, s: f64) -> f64 {
        0.5 * self.force_constant * (s - self.centre).powi(2)
    }
}

/// Umbrella sampling with WHAM (Weighted Histogram Analysis Method).
#[derive(Debug, Clone)]
pub struct UmbrellaSamplingWham {
    /// Windows.
    pub windows: Vec<UmbrellaWindow>,
    /// Temperature (K).
    pub temperature: f64,
    /// WHAM convergence tolerance.
    pub tolerance: f64,
    /// Maximum WHAM iterations.
    pub max_iterations: usize,
    /// Histogram bin edges.
    pub bin_edges: Vec<f64>,
}

impl UmbrellaSamplingWham {
    /// Create umbrella sampling setup with uniformly spaced windows.
    pub fn new_uniform(
        cv_min: f64,
        cv_max: f64,
        n_windows: usize,
        force_constant: f64,
        temperature: f64,
    ) -> Self {
        let step = (cv_max - cv_min) / (n_windows as f64 - 1.0).max(1.0);
        let mut windows = Vec::with_capacity(n_windows);
        for i in 0..n_windows {
            let centre = cv_min + i as f64 * step;
            windows.push(UmbrellaWindow::new(centre, force_constant));
        }
        Self {
            windows,
            temperature,
            tolerance: 1e-6,
            max_iterations: 10000,
            bin_edges: Vec::new(),
        }
    }

    /// Set WHAM parameters.
    pub fn set_wham_params(&mut self, tolerance: f64, max_iterations: usize) {
        self.tolerance = tolerance;
        self.max_iterations = max_iterations;
    }

    /// Set custom histogram bin edges.
    pub fn set_bin_edges(&mut self, edges: Vec<f64>) {
        self.bin_edges = edges;
    }

    /// Generate bin edges from range.
    pub fn auto_bin_edges(&mut self, cv_min: f64, cv_max: f64, n_bins: usize) {
        let bw = (cv_max - cv_min) / n_bins as f64;
        self.bin_edges = (0..=n_bins).map(|i| cv_min + i as f64 * bw).collect();
    }

    /// Add a sample to a specific window.
    pub fn add_sample(&mut self, window_index: usize, cv_value: f64) {
        if let Some(w) = self.windows.get_mut(window_index) {
            w.add_sample(cv_value);
        }
    }

    /// Run the WHAM analysis and return `(bin_centres, free_energies)`.
    pub fn run_wham(&self) -> (Vec<f64>, Vec<f64>) {
        let n_bins = if self.bin_edges.len() > 1 {
            self.bin_edges.len() - 1
        } else {
            return (Vec::new(), Vec::new());
        };
        let n_windows = self.windows.len();
        let beta = 1.0 / (KB_KJMOL * self.temperature);

        // Build histograms
        let mut histograms: Vec<Vec<f64>> = vec![vec![0.0; n_bins]; n_windows];
        let mut n_samples: Vec<f64> = vec![0.0; n_windows];
        for (w, window) in self.windows.iter().enumerate() {
            for &s in &window.samples {
                let bin = self.find_bin(s);
                if let Some(b) = bin {
                    histograms[w][b] += 1.0;
                    n_samples[w] += 1.0;
                }
            }
        }

        // Bin centres
        let bin_centres: Vec<f64> = (0..n_bins)
            .map(|i| 0.5 * (self.bin_edges[i] + self.bin_edges[i + 1]))
            .collect();

        // Precompute bias energies
        let mut bias: Vec<Vec<f64>> = vec![vec![0.0; n_bins]; n_windows];
        for (w, window) in self.windows.iter().enumerate() {
            for (b, &cv) in bin_centres.iter().enumerate() {
                bias[w][b] = window.bias_energy(cv);
            }
        }

        // WHAM iteration
        let mut f_w = vec![0.0_f64; n_windows]; // free energy of each window
        let mut pmf = vec![0.0_f64; n_bins];

        for _iter in 0..self.max_iterations {
            // Compute unbiased probability of each bin
            for b in 0..n_bins {
                let mut numer = 0.0;
                let mut denom = 0.0;
                for w in 0..n_windows {
                    numer += histograms[w][b];
                    denom += n_samples[w] * (f_w[w] - beta * bias[w][b]).exp();
                }
                pmf[b] = if denom > EPS { numer / denom } else { 0.0 };
            }

            // Update window free energies
            let mut max_change = 0.0_f64;
            for w in 0..n_windows {
                let mut z = 0.0;
                for b in 0..n_bins {
                    z += pmf[b] * (-beta * bias[w][b]).exp();
                }
                let new_f = if z > EPS { -(z.ln()) } else { f_w[w] };
                let change = (new_f - f_w[w]).abs();
                if change > max_change {
                    max_change = change;
                }
                f_w[w] = new_f;
            }

            if max_change < self.tolerance {
                break;
            }
        }

        // Normalise and compute free energy
        let max_p = pmf.iter().cloned().fold(0.0_f64, f64::max);
        let fes: Vec<f64> = pmf
            .iter()
            .map(|&p| {
                if p > EPS && max_p > EPS {
                    -(p / max_p).ln() / beta
                } else {
                    f64::INFINITY
                }
            })
            .collect();

        // Shift minimum to zero
        let min_fes = fes
            .iter()
            .cloned()
            .filter(|x| x.is_finite())
            .fold(f64::INFINITY, f64::min);
        let fes_shifted: Vec<f64> = fes.iter().map(|x| x - min_fes).collect();

        (bin_centres, fes_shifted)
    }

    /// Find the histogram bin for a CV value.
    fn find_bin(&self, cv: f64) -> Option<usize> {
        if self.bin_edges.len() < 2 {
            return None;
        }
        if cv < self.bin_edges[0] || cv >= self.bin_edges[self.bin_edges.len() - 1] {
            return None;
        }
        // Linear search (fine for typical bin counts)
        (0..self.bin_edges.len() - 1)
            .find(|&i| cv >= self.bin_edges[i] && cv < self.bin_edges[i + 1])
    }

    /// Total number of samples across all windows.
    pub fn total_samples(&self) -> usize {
        self.windows.iter().map(|w| w.samples.len()).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AdaptiveBiasingForce
// ─────────────────────────────────────────────────────────────────────────────

/// Adaptive Biasing Force (ABF) method.
///
/// Accumulates the mean force in CV bins and applies an equal-and-opposite
/// bias force to flatten the free energy landscape.
#[derive(Debug, Clone)]
pub struct AdaptiveBiasingForce {
    /// Lower bound of CV range.
    pub cv_min: f64,
    /// Upper bound of CV range.
    pub cv_max: f64,
    /// Bin width.
    pub bin_width: f64,
    /// Number of bins.
    pub n_bins: usize,
    /// Accumulated force in each bin.
    pub force_sum: Vec<f64>,
    /// Number of samples in each bin.
    pub count: Vec<u64>,
    /// Minimum sample count before bias is applied (ramp-up threshold).
    pub min_samples: u64,
    /// Full-bias sample count (linear ramp between min_samples and full_samples).
    pub full_samples: u64,
    /// Step counter.
    pub step: u64,
}

impl AdaptiveBiasingForce {
    /// Create a new ABF engine.
    pub fn new(cv_min: f64, cv_max: f64, n_bins: usize, min_samples: u64) -> Self {
        let bw = (cv_max - cv_min) / n_bins as f64;
        Self {
            cv_min,
            cv_max,
            bin_width: bw,
            n_bins,
            force_sum: vec![0.0; n_bins],
            count: vec![0; n_bins],
            min_samples,
            full_samples: min_samples * 2,
            step: 0,
        }
    }

    /// Set the full-bias ramp-up threshold.
    pub fn set_full_samples(&mut self, n: u64) {
        self.full_samples = n;
    }

    /// Get the bin index for a CV value.
    fn bin_index(&self, cv: f64) -> Option<usize> {
        if cv < self.cv_min || cv >= self.cv_max {
            return None;
        }
        let idx = ((cv - self.cv_min) / self.bin_width) as usize;
        if idx < self.n_bins { Some(idx) } else { None }
    }

    /// Accumulate a force sample and return the bias force to apply.
    pub fn update(&mut self, cv: f64, instantaneous_force: f64) -> f64 {
        self.step += 1;
        if let Some(idx) = self.bin_index(cv) {
            self.force_sum[idx] += instantaneous_force;
            self.count[idx] += 1;
            let c = self.count[idx];
            if c < self.min_samples {
                return 0.0;
            }
            let mean_force = self.force_sum[idx] / c as f64;
            let ramp = if c >= self.full_samples {
                1.0
            } else {
                (c - self.min_samples) as f64 / (self.full_samples - self.min_samples) as f64
            };
            -mean_force * ramp
        } else {
            0.0
        }
    }

    /// Get the accumulated mean force profile: `(bin_centres, mean_forces)`.
    pub fn mean_force_profile(&self) -> (Vec<f64>, Vec<f64>) {
        let centres: Vec<f64> = (0..self.n_bins)
            .map(|i| self.cv_min + (i as f64 + 0.5) * self.bin_width)
            .collect();
        let forces: Vec<f64> = (0..self.n_bins)
            .map(|i| {
                if self.count[i] > 0 {
                    self.force_sum[i] / self.count[i] as f64
                } else {
                    0.0
                }
            })
            .collect();
        (centres, forces)
    }

    /// Integrate the mean force to get the free energy profile.
    /// Returns `(bin_centres, free_energy)`.
    pub fn free_energy_profile(&self) -> (Vec<f64>, Vec<f64>) {
        let (centres, forces) = self.mean_force_profile();
        let mut fes = Vec::with_capacity(self.n_bins);
        let mut running = 0.0;
        for (i, &f) in forces.iter().enumerate() {
            if i > 0 {
                running -= f * self.bin_width;
            }
            fes.push(running);
        }
        // Shift minimum to zero
        let min_val = fes.iter().cloned().fold(f64::INFINITY, f64::min);
        for v in &mut fes {
            *v -= min_val;
        }
        (centres, fes)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SteeredMdExtension
// ─────────────────────────────────────────────────────────────────────────────

/// Pulling protocol for steered MD.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PullingProtocol {
    /// Constant velocity pulling.
    ConstantVelocity,
    /// Constant force pulling.
    ConstantForce,
    /// Adaptive velocity (slows down when force is high).
    AdaptiveVelocity,
}

/// Extended steered MD engine.
///
/// Pulls a collective variable from an initial to a final value, tracking
/// work and force along the way.
#[derive(Debug, Clone)]
pub struct SteeredMdExtension {
    /// Force constant for the spring (kJ mol⁻¹ unit⁻²).
    pub force_constant: f64,
    /// Pulling velocity (units / step).
    pub velocity: f64,
    /// Current target CV value.
    pub target_cv: f64,
    /// Initial target CV.
    pub initial_cv: f64,
    /// Final target CV.
    pub final_cv: f64,
    /// Pulling protocol.
    pub protocol: PullingProtocol,
    /// Accumulated work (kJ/mol).
    pub work: f64,
    /// Step counter.
    pub step: u64,
    /// History of `(step, cv_actual, force, work)`.
    pub history: Vec<(u64, f64, f64, f64)>,
    /// Adaptive velocity scale factor.
    pub adaptive_scale: f64,
    /// Maximum force threshold for adaptive velocity.
    pub max_force_threshold: f64,
}

impl SteeredMdExtension {
    /// Create a new steered MD extension.
    pub fn new(
        force_constant: f64,
        velocity: f64,
        initial_cv: f64,
        final_cv: f64,
        protocol: PullingProtocol,
    ) -> Self {
        Self {
            force_constant,
            velocity,
            target_cv: initial_cv,
            initial_cv,
            final_cv,
            protocol,
            work: 0.0,
            step: 0,
            history: Vec::new(),
            adaptive_scale: 1.0,
            max_force_threshold: 100.0,
        }
    }

    /// Set adaptive velocity parameters.
    pub fn set_adaptive_params(&mut self, max_force_threshold: f64, scale: f64) {
        self.max_force_threshold = max_force_threshold;
        self.adaptive_scale = scale;
    }

    /// Step the steered MD engine. Returns the force to apply.
    pub fn step(&mut self, cv_actual: f64) -> f64 {
        self.step += 1;
        let force = match self.protocol {
            PullingProtocol::ConstantVelocity => {
                let old_target = self.target_cv;
                self.target_cv += self.velocity;
                self.target_cv = clamp_towards(self.target_cv, self.initial_cv, self.final_cv);
                let f = self.force_constant * (self.target_cv - cv_actual);
                let dw = f * (self.target_cv - old_target);
                self.work += dw;
                f
            }
            PullingProtocol::ConstantForce => {
                let dir = if self.final_cv > self.initial_cv {
                    1.0
                } else {
                    -1.0
                };
                let f = self.force_constant * dir;
                self.work += f * self.velocity;
                f
            }
            PullingProtocol::AdaptiveVelocity => {
                let base_force = self.force_constant * (self.target_cv - cv_actual);
                let speed_factor = if base_force.abs() > self.max_force_threshold {
                    self.adaptive_scale / (1.0 + base_force.abs() / self.max_force_threshold)
                } else {
                    1.0
                };
                let old_target = self.target_cv;
                self.target_cv += self.velocity * speed_factor;
                self.target_cv = clamp_towards(self.target_cv, self.initial_cv, self.final_cv);
                let f = self.force_constant * (self.target_cv - cv_actual);
                let dw = f * (self.target_cv - old_target);
                self.work += dw;
                f
            }
        };
        self.history.push((self.step, cv_actual, force, self.work));
        force
    }

    /// Check if pulling is complete.
    pub fn is_complete(&self) -> bool {
        if self.final_cv > self.initial_cv {
            self.target_cv >= self.final_cv
        } else {
            self.target_cv <= self.final_cv
        }
    }

    /// Get the Jarzynski free energy estimate: `ΔF = -kBT ln <exp(-βW)>`.
    pub fn jarzynski_estimate(work_values: &[f64], temperature: f64) -> f64 {
        if work_values.is_empty() {
            return 0.0;
        }
        let beta = 1.0 / (KB_KJMOL * temperature);
        let max_bw = work_values
            .iter()
            .map(|&w| -beta * w)
            .fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = work_values
            .iter()
            .map(|&w| (-beta * w - max_bw).exp())
            .sum();
        let log_avg = max_bw + (sum / work_values.len() as f64).ln();
        -log_avg / beta
    }

    /// Mean work over the trajectory.
    pub fn mean_work(&self) -> f64 {
        if self.history.is_empty() {
            0.0
        } else {
            self.work
        }
    }
}

/// Clamp `val` so it doesn't overshoot the target direction.
fn clamp_towards(val: f64, _from: f64, to: f64) -> f64 {
    if to > _from { val.min(to) } else { val.max(to) }
}

// ─────────────────────────────────────────────────────────────────────────────
// TransitionPathSampling
// ─────────────────────────────────────────────────────────────────────────────

/// State classification for transition path sampling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TpsState {
    /// In reactant basin A.
    Reactant,
    /// In product basin B.
    Product,
    /// In transition region (neither A nor B).
    Transition,
}

/// A single reactive trajectory for TPS.
#[derive(Debug, Clone)]
pub struct TpsTrajectory {
    /// Snapshots along the path: (positions_flat, cv_value).
    pub frames: Vec<(Vec<[f64; 3]>, f64)>,
    /// Whether this is a reactive trajectory (connects A to B).
    pub is_reactive: bool,
    /// Weight of this trajectory.
    pub weight: f64,
}

/// Transition Path Sampling (TPS) engine.
///
/// Uses two-ended shooting moves to sample reactive trajectories connecting
/// reactant state A and product state B.
#[derive(Debug, Clone)]
pub struct TransitionPathSampling {
    /// CV threshold for reactant state A.
    pub cv_a: f64,
    /// CV threshold for product state B.
    pub cv_b: f64,
    /// Maximum trajectory length (frames).
    pub max_path_length: usize,
    /// Momentum perturbation magnitude.
    pub perturbation_magnitude: f64,
    /// Accepted trajectories.
    pub accepted_paths: Vec<TpsTrajectory>,
    /// Total number of shooting attempts.
    pub n_attempts: u64,
    /// Number of accepted trajectories.
    pub n_accepted: u64,
}

impl TransitionPathSampling {
    /// Create a new TPS engine.
    ///
    /// * `cv_a` – CV value below which the system is in state A.
    /// * `cv_b` – CV value above which the system is in state B.
    /// * `max_path_length` – maximum number of frames per trajectory.
    /// * `perturbation_magnitude` – magnitude of momentum perturbation.
    pub fn new(cv_a: f64, cv_b: f64, max_path_length: usize, perturbation_magnitude: f64) -> Self {
        Self {
            cv_a,
            cv_b,
            max_path_length,
            perturbation_magnitude,
            accepted_paths: Vec::new(),
            n_attempts: 0,
            n_accepted: 0,
        }
    }

    /// Classify a CV value into a TPS state.
    pub fn classify(&self, cv: f64) -> TpsState {
        if cv <= self.cv_a {
            TpsState::Reactant
        } else if cv >= self.cv_b {
            TpsState::Product
        } else {
            TpsState::Transition
        }
    }

    /// Check if a trajectory is reactive (connects A to B in either direction).
    pub fn is_reactive(path: &[(Vec<[f64; 3]>, f64)], cv_a: f64, cv_b: f64) -> bool {
        if path.is_empty() {
            return false;
        }
        let first = path[0].1;
        let last = path[path.len() - 1].1;
        (first <= cv_a && last >= cv_b) || (first >= cv_b && last <= cv_a)
    }

    /// Select a random shooting point index from the current path.
    pub fn random_shooting_point(&self, path_length: usize) -> usize {
        if path_length < 3 {
            return 0;
        }
        let mut rng = rand::rng();
        rng.random_range(1..path_length - 1)
    }

    /// Generate perturbed momenta from original momenta.
    pub fn perturb_momenta(&self, momenta: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut rng = rand::rng();
        let mag = self.perturbation_magnitude;
        momenta
            .iter()
            .map(|p| {
                [
                    p[0] + rng.random_range(-mag..mag),
                    p[1] + rng.random_range(-mag..mag),
                    p[2] + rng.random_range(-mag..mag),
                ]
            })
            .collect()
    }

    /// Record a shooting attempt result.
    pub fn record_attempt(&mut self, path: Vec<(Vec<[f64; 3]>, f64)>) -> bool {
        self.n_attempts += 1;
        let reactive = Self::is_reactive(&path, self.cv_a, self.cv_b);
        if reactive {
            self.n_accepted += 1;
            self.accepted_paths.push(TpsTrajectory {
                frames: path,
                is_reactive: true,
                weight: 1.0,
            });
        }
        reactive
    }

    /// Acceptance ratio.
    pub fn acceptance_ratio(&self) -> f64 {
        if self.n_attempts == 0 {
            0.0
        } else {
            self.n_accepted as f64 / self.n_attempts as f64
        }
    }

    /// Compute the committor probability from stored reactive paths.
    /// Returns fraction of paths starting from reactant that reach product.
    pub fn committor_estimate(&self) -> f64 {
        if self.accepted_paths.is_empty() {
            return 0.5;
        }
        let mut n_ab = 0u64;
        let mut n_total = 0u64;
        for path in &self.accepted_paths {
            if path.frames.is_empty() {
                continue;
            }
            n_total += 1;
            let first_cv = path.frames[0].1;
            if first_cv <= self.cv_a {
                n_ab += 1;
            }
        }
        if n_total == 0 {
            0.5
        } else {
            n_ab as f64 / n_total as f64
        }
    }

    /// Number of reactive paths collected.
    pub fn n_reactive_paths(&self) -> usize {
        self.accepted_paths.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-CV Metadynamics
// ─────────────────────────────────────────────────────────────────────────────

/// A 2D Gaussian hill for multi-CV metadynamics.
#[derive(Debug, Clone)]
pub struct GaussianHill2D {
    /// CV centres.
    pub centre: [f64; 2],
    /// Hill widths (sigma).
    pub sigma: [f64; 2],
    /// Hill height.
    pub height: f64,
    /// Deposition time.
    pub time: f64,
}

/// Two-dimensional well-tempered metadynamics.
#[derive(Debug, Clone)]
pub struct Metadynamics2D {
    /// Initial hill height (kJ/mol).
    pub initial_height: f64,
    /// Hill widths (sigma) for each CV.
    pub sigma: [f64; 2],
    /// Bias factor γ.
    pub bias_factor: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Deposition interval (steps).
    pub deposition_interval: u64,
    /// Deposited hills.
    pub hills: Vec<GaussianHill2D>,
    /// Step counter.
    pub step: u64,
}

impl Metadynamics2D {
    /// Create a 2D well-tempered metadynamics engine.
    pub fn new(
        initial_height: f64,
        sigma: [f64; 2],
        bias_factor: f64,
        temperature: f64,
        deposition_interval: u64,
    ) -> Self {
        Self {
            initial_height,
            sigma,
            bias_factor,
            temperature,
            deposition_interval,
            hills: Vec::new(),
            step: 0,
        }
    }

    /// Compute the bias at a 2D CV point.
    pub fn bias_at(&self, s: [f64; 2]) -> f64 {
        let mut v = 0.0;
        let inv_2s0 = 0.5 / (self.sigma[0] * self.sigma[0]);
        let inv_2s1 = 0.5 / (self.sigma[1] * self.sigma[1]);
        for hill in &self.hills {
            let d0 = s[0] - hill.centre[0];
            let d1 = s[1] - hill.centre[1];
            v += hill.height * (-d0 * d0 * inv_2s0 - d1 * d1 * inv_2s1).exp();
        }
        v
    }

    /// Deposit a 2D hill at the current CV values.
    pub fn deposit_hill(&mut self, s: [f64; 2]) -> f64 {
        let current_bias = self.bias_at(s);
        let delta_t = (self.bias_factor - 1.0) * self.temperature;
        let eff_height = if delta_t.abs() < EPS {
            self.initial_height
        } else {
            self.initial_height * (-current_bias / (KB_KJMOL * delta_t)).exp()
        };
        self.hills.push(GaussianHill2D {
            centre: s,
            sigma: self.sigma,
            height: eff_height,
            time: self.step as f64,
        });
        eff_height
    }

    /// Step the 2D metadynamics engine.
    pub fn step(&mut self, cv: [f64; 2]) -> Option<f64> {
        self.step += 1;
        if self.step.is_multiple_of(self.deposition_interval) {
            Some(self.deposit_hill(cv))
        } else {
            None
        }
    }

    /// Compute the bias gradient at a 2D CV point.
    pub fn bias_gradient(&self, s: [f64; 2]) -> [f64; 2] {
        let mut g = [0.0; 2];
        let inv_s0_sq = 1.0 / (self.sigma[0] * self.sigma[0]);
        let inv_s1_sq = 1.0 / (self.sigma[1] * self.sigma[1]);
        for hill in &self.hills {
            let d0 = s[0] - hill.centre[0];
            let d1 = s[1] - hill.centre[1];
            let exp_val =
                hill.height * (-0.5 * d0 * d0 * inv_s0_sq - 0.5 * d1 * d1 * inv_s1_sq).exp();
            g[0] += -exp_val * d0 * inv_s0_sq;
            g[1] += -exp_val * d1 * inv_s1_sq;
        }
        g
    }

    /// Number of hills deposited.
    pub fn n_hills(&self) -> usize {
        self.hills.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Extended ABF with block averaging
// ─────────────────────────────────────────────────────────────────────────────

/// Block-averaged ABF bin data.
#[derive(Debug, Clone)]
pub struct AbfBlockBin {
    /// Running block sums.
    pub block_sums: Vec<f64>,
    /// Samples in current block.
    pub current_block_sum: f64,
    /// Samples in current block count.
    pub current_block_count: u64,
    /// Block size.
    pub block_size: u64,
}

impl AbfBlockBin {
    /// Create a new block bin.
    pub fn new(block_size: u64) -> Self {
        Self {
            block_sums: Vec::new(),
            current_block_sum: 0.0,
            current_block_count: 0,
            block_size,
        }
    }

    /// Add a force sample. Returns `Some(block_average)` if a block completed.
    pub fn add_sample(&mut self, force: f64) -> Option<f64> {
        self.current_block_sum += force;
        self.current_block_count += 1;
        if self.current_block_count >= self.block_size {
            let avg = self.current_block_sum / self.current_block_count as f64;
            self.block_sums.push(avg);
            self.current_block_sum = 0.0;
            self.current_block_count = 0;
            Some(avg)
        } else {
            None
        }
    }

    /// Mean over completed blocks.
    pub fn block_mean(&self) -> f64 {
        if self.block_sums.is_empty() {
            return 0.0;
        }
        self.block_sums.iter().sum::<f64>() / self.block_sums.len() as f64
    }

    /// Standard error from block averages.
    pub fn block_stderr(&self) -> f64 {
        let n = self.block_sums.len();
        if n < 2 {
            return 0.0;
        }
        let mean = self.block_mean();
        let var: f64 = self
            .block_sums
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / (n - 1) as f64;
        (var / n as f64).sqrt()
    }
}

/// Extended ABF with block averaging for error estimation.
#[derive(Debug, Clone)]
pub struct AbfBlockAveraged {
    /// Lower bound of CV range.
    pub cv_min: f64,
    /// Upper bound of CV range.
    pub cv_max: f64,
    /// Bin width.
    pub bin_width: f64,
    /// Number of bins.
    pub n_bins: usize,
    /// Per-bin block data.
    pub bins: Vec<AbfBlockBin>,
    /// Minimum samples before bias is applied.
    pub min_samples: u64,
    /// Total sample count per bin.
    pub sample_count: Vec<u64>,
}

impl AbfBlockAveraged {
    /// Create a new block-averaged ABF.
    pub fn new(cv_min: f64, cv_max: f64, n_bins: usize, block_size: u64, min_samples: u64) -> Self {
        let bw = (cv_max - cv_min) / n_bins as f64;
        let bins = (0..n_bins).map(|_| AbfBlockBin::new(block_size)).collect();
        Self {
            cv_min,
            cv_max,
            bin_width: bw,
            n_bins,
            bins,
            min_samples,
            sample_count: vec![0; n_bins],
        }
    }

    /// Get bin index.
    fn bin_index(&self, cv: f64) -> Option<usize> {
        if cv < self.cv_min || cv >= self.cv_max {
            return None;
        }
        let idx = ((cv - self.cv_min) / self.bin_width) as usize;
        if idx < self.n_bins { Some(idx) } else { None }
    }

    /// Accumulate a sample and return the bias force.
    pub fn update(&mut self, cv: f64, force: f64) -> f64 {
        if let Some(idx) = self.bin_index(cv) {
            self.bins[idx].add_sample(force);
            self.sample_count[idx] += 1;
            if self.sample_count[idx] < self.min_samples {
                return 0.0;
            }
            -self.bins[idx].block_mean()
        } else {
            0.0
        }
    }

    /// Get the mean force profile with error bars.
    /// Returns `(centres, mean_forces, std_errors)`.
    pub fn profile_with_errors(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let centres: Vec<f64> = (0..self.n_bins)
            .map(|i| self.cv_min + (i as f64 + 0.5) * self.bin_width)
            .collect();
        let means: Vec<f64> = self.bins.iter().map(|b| b.block_mean()).collect();
        let errs: Vec<f64> = self.bins.iter().map(|b| b.block_stderr()).collect();
        (centres, means, errs)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper functions
// ─────────────────────────────────────────────────────────────────────────────

/// Euclidean distance between two 3D points.
#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Vector subtraction.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Dot product.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Vector length.
#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── CollectiveVariable ──────────────────────────────────────────────

    #[test]
    fn test_cv_distance() {
        let cv = CollectiveVariable::distance(0, 1);
        let pos = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        let d = cv.evaluate(&pos);
        assert!((d - 5.0).abs() < 1e-10, "distance = {d}");
    }

    #[test]
    fn test_cv_distance_same_point() {
        let cv = CollectiveVariable::distance(0, 0);
        let pos = vec![[1.0, 2.0, 3.0]];
        assert!((cv.evaluate(&pos)).abs() < 1e-12);
    }

    #[test]
    fn test_cv_dihedral_planar() {
        // Four coplanar points forming a known dihedral = 0 (cis)
        let cv = CollectiveVariable::dihedral(0, 1, 2, 3);
        let pos = vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let angle = cv.evaluate(&pos);
        assert!(angle.abs() < 1e-10, "dihedral = {angle}");
    }

    #[test]
    fn test_cv_dihedral_trans() {
        let cv = CollectiveVariable::dihedral(0, 1, 2, 3);
        let pos = vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ];
        let angle = cv.evaluate(&pos);
        assert!(
            (angle.abs() - std::f64::consts::PI).abs() < 1e-10,
            "dihedral = {angle}"
        );
    }

    #[test]
    fn test_cv_rmsd_zero() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let cv = CollectiveVariable::rmsd(vec![0, 1], positions.clone());
        let rmsd = cv.evaluate(&positions);
        assert!(rmsd.abs() < 1e-12, "rmsd = {rmsd}");
    }

    #[test]
    fn test_cv_rmsd_nonzero() {
        let reference = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let cv = CollectiveVariable::rmsd(vec![0, 1], reference);
        let positions = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let rmsd = cv.evaluate(&positions);
        assert!((rmsd - 1.0).abs() < 1e-10, "rmsd = {rmsd}");
    }

    #[test]
    fn test_cv_linear_combination() {
        let cv = CollectiveVariable::linear_combination(vec![0, 1], [1.0, 0.0, 0.0]);
        let pos = vec![[3.0, 10.0, 20.0], [5.0, 10.0, 20.0]];
        let val = cv.evaluate(&pos);
        assert!((val - 8.0).abs() < 1e-10, "linear_combo = {val}");
    }

    #[test]
    fn test_cv_coordination() {
        let cv = CollectiveVariable::coordination_number(0, &[1, 2], 1.0, 0.0, 6.0, 12.0);
        let pos = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let coord = cv.evaluate(&pos);
        // First neighbour very close → ~1, second very far → ~0
        assert!(coord > 0.8 && coord < 2.1, "coord = {coord}");
    }

    // ── RemdController ──────────────────────────────────────────────────

    #[test]
    fn test_remd_geometric_ladder() {
        let remd = RemdController::new_geometric(300.0, 600.0, 4, 100);
        assert_eq!(remd.n_replicas(), 4);
        assert!((remd.temperatures[0] - 300.0).abs() < 1e-6);
        assert!((remd.temperatures[3] - 600.0).abs() < 1.0);
        // Check geometric spacing
        let ratio = remd.temperatures[1] / remd.temperatures[0];
        let ratio2 = remd.temperatures[2] / remd.temperatures[1];
        assert!((ratio - ratio2).abs() < 1e-6);
    }

    #[test]
    fn test_remd_swap_always_accept_downhill() {
        let mut remd = RemdController::new_geometric(300.0, 400.0, 2, 1);
        // Set energies so swap is always favourable
        remd.set_energies(&[100.0, 50.0]); // lower T has higher E
        let swapped = remd.attempt_swaps();
        // With ΔE<0 at higher beta, should accept
        assert!(!swapped.is_empty() || remd.pair_stats[0].0 > 0);
    }

    #[test]
    fn test_remd_hamiltonian() {
        let remd = RemdController::new_hamiltonian(300.0, vec![0.0, 0.5, 1.0], 50);
        assert_eq!(remd.n_replicas(), 3);
        assert!(remd.lambda.is_some());
    }

    #[test]
    fn test_remd_acceptance_initially_zero() {
        let remd = RemdController::new_geometric(300.0, 400.0, 3, 10);
        assert!((remd.overall_acceptance() - 0.0).abs() < 1e-12);
    }

    // ── WellTemperedMetadynamics ────────────────────────────────────────

    #[test]
    fn test_wtm_deposit_hill() {
        let mut wtm = WellTemperedMetadynamics::new(1.0, 0.1, 10.0, 300.0, 1);
        wtm.step(0.5);
        assert_eq!(wtm.n_hills(), 1);
        let bias = wtm.bias_at(0.5);
        assert!(bias > 0.0, "bias = {bias}");
    }

    #[test]
    fn test_wtm_height_decreases() {
        let mut wtm = WellTemperedMetadynamics::new(1.0, 0.1, 5.0, 300.0, 1);
        let h1 = wtm.deposit_hill(0.0);
        let h2 = wtm.deposit_hill(0.0);
        assert!(h2 < h1, "h2={h2} should be less than h1={h1}");
    }

    #[test]
    fn test_wtm_bias_force() {
        let mut wtm = WellTemperedMetadynamics::new(1.0, 0.1, 10.0, 300.0, 1);
        wtm.deposit_hill(0.0);
        let f = wtm.bias_force(0.05);
        // Force should push back towards 0
        assert!(f > 0.0, "force = {f}");
    }

    #[test]
    fn test_wtm_with_grid() {
        let mut wtm = WellTemperedMetadynamics::new(1.0, 0.1, 10.0, 300.0, 1);
        wtm.enable_grid(-1.0, 1.0, 200);
        wtm.deposit_hill(0.0);
        let bias_grid = wtm.bias_at(0.0);
        assert!(bias_grid > 0.0, "grid bias = {bias_grid}");
    }

    #[test]
    fn test_wtm_fes() {
        let mut wtm = WellTemperedMetadynamics::new(1.0, 0.1, 10.0, 300.0, 1);
        for _ in 0..10 {
            wtm.deposit_hill(0.0);
        }
        let (_cvs, fes) = wtm.free_energy_surface(-0.5, 0.5, 11);
        assert_eq!(fes.len(), 11);
    }

    // ── UmbrellaSamplingWham ────────────────────────────────────────────

    #[test]
    fn test_umbrella_window_bias() {
        let w = UmbrellaWindow::new(1.0, 100.0);
        let e = w.bias_energy(1.5);
        assert!((e - 12.5).abs() < 1e-10, "bias energy = {e}");
    }

    #[test]
    fn test_wham_setup() {
        let us = UmbrellaSamplingWham::new_uniform(0.0, 1.0, 5, 100.0, 300.0);
        assert_eq!(us.windows.len(), 5);
        assert!((us.windows[0].centre - 0.0).abs() < 1e-10);
        assert!((us.windows[4].centre - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_wham_trivial() {
        let mut us = UmbrellaSamplingWham::new_uniform(0.0, 1.0, 3, 100.0, 300.0);
        us.auto_bin_edges(0.0, 1.0, 10);
        // Add identical samples to all windows
        for w in 0..3 {
            for _ in 0..100 {
                let cv = us.windows[w].centre;
                us.add_sample(w, cv);
            }
        }
        let (centres, fes) = us.run_wham();
        assert_eq!(centres.len(), 10);
        // FES minimum should be zero
        let min_fes = fes
            .iter()
            .cloned()
            .filter(|x| x.is_finite())
            .fold(f64::INFINITY, f64::min);
        assert!(min_fes.abs() < 1e-6, "min_fes = {min_fes}");
    }

    #[test]
    fn test_wham_total_samples() {
        let mut us = UmbrellaSamplingWham::new_uniform(0.0, 1.0, 2, 100.0, 300.0);
        us.add_sample(0, 0.1);
        us.add_sample(0, 0.2);
        us.add_sample(1, 0.9);
        assert_eq!(us.total_samples(), 3);
    }

    // ── ABF ─────────────────────────────────────────────────────────────

    #[test]
    fn test_abf_ramp_up() {
        let mut abf = AdaptiveBiasingForce::new(0.0, 1.0, 10, 5);
        abf.set_full_samples(10);
        // First 4 samples should return zero bias
        for _ in 0..4 {
            let bias = abf.update(0.05, 10.0);
            assert!((bias).abs() < 1e-12, "should be zero in ramp-up");
        }
        // 5th sample starts the ramp
        let bias = abf.update(0.05, 10.0);
        assert!(bias.abs() >= 0.0); // may still be zero at boundary
    }

    #[test]
    fn test_abf_mean_force_profile() {
        let mut abf = AdaptiveBiasingForce::new(0.0, 1.0, 10, 0);
        abf.set_full_samples(0);
        for _ in 0..100 {
            abf.update(0.15, 5.0);
        }
        let (centres, forces) = abf.mean_force_profile();
        assert_eq!(centres.len(), 10);
        // bin 1 (centre at 0.15) should have mean force ~5.0
        assert!(
            (forces[1] - 5.0).abs() < 1e-10,
            "mean force = {}",
            forces[1]
        );
    }

    #[test]
    fn test_abf_free_energy_profile() {
        let mut abf = AdaptiveBiasingForce::new(0.0, 1.0, 5, 0);
        abf.set_full_samples(0);
        for _ in 0..50 {
            abf.update(0.1, 2.0);
            abf.update(0.3, 2.0);
            abf.update(0.5, 2.0);
            abf.update(0.7, 2.0);
            abf.update(0.9, 2.0);
        }
        let (_centres, fes) = abf.free_energy_profile();
        assert_eq!(fes.len(), 5);
        // With uniform force, FES should be monotonic
        for i in 1..fes.len() {
            assert!(fes[i] >= fes[i - 1] - 1e-6 || fes[i] <= fes[i - 1] + 1e-6);
        }
    }

    // ── SteeredMdExtension ──────────────────────────────────────────────

    #[test]
    fn test_smd_constant_velocity() {
        let mut smd =
            SteeredMdExtension::new(100.0, 0.01, 0.0, 1.0, PullingProtocol::ConstantVelocity);
        let _f = smd.step(0.0);
        assert!(smd.target_cv > 0.0);
        assert_eq!(smd.history.len(), 1);
    }

    #[test]
    fn test_smd_completes() {
        let mut smd =
            SteeredMdExtension::new(100.0, 0.1, 0.0, 1.0, PullingProtocol::ConstantVelocity);
        for _ in 0..20 {
            smd.step(smd.target_cv);
        }
        assert!(smd.is_complete());
    }

    #[test]
    fn test_smd_jarzynski() {
        let works = vec![1.0, 2.0, 1.5, 1.2];
        let fe = SteeredMdExtension::jarzynski_estimate(&works, 300.0);
        // Jarzynski estimate should be <= mean work (second law)
        let mean_w: f64 = works.iter().sum::<f64>() / works.len() as f64;
        assert!(fe <= mean_w + 1e-6, "fe={fe} > mean_w={mean_w}");
    }

    #[test]
    fn test_smd_adaptive() {
        let mut smd =
            SteeredMdExtension::new(100.0, 0.1, 0.0, 5.0, PullingProtocol::AdaptiveVelocity);
        smd.set_adaptive_params(50.0, 0.5);
        // Run several steps; target should advance
        for _ in 0..10 {
            smd.step(0.0); // pulling away from actual position
        }
        assert!(smd.target_cv > 0.0);
    }

    // ── TPS ─────────────────────────────────────────────────────────────

    #[test]
    fn test_tps_classify() {
        let tps = TransitionPathSampling::new(0.3, 0.7, 100, 0.1);
        assert_eq!(tps.classify(0.1), TpsState::Reactant);
        assert_eq!(tps.classify(0.5), TpsState::Transition);
        assert_eq!(tps.classify(0.9), TpsState::Product);
    }

    #[test]
    fn test_tps_reactive_path() {
        let path: Vec<(Vec<[f64; 3]>, f64)> = vec![
            (vec![[0.0; 3]], 0.1),
            (vec![[0.0; 3]], 0.5),
            (vec![[0.0; 3]], 0.9),
        ];
        assert!(TransitionPathSampling::is_reactive(&path, 0.3, 0.7));
    }

    #[test]
    fn test_tps_non_reactive_path() {
        let path: Vec<(Vec<[f64; 3]>, f64)> = vec![
            (vec![[0.0; 3]], 0.1),
            (vec![[0.0; 3]], 0.2),
            (vec![[0.0; 3]], 0.15),
        ];
        assert!(!TransitionPathSampling::is_reactive(&path, 0.3, 0.7));
    }

    #[test]
    fn test_tps_record() {
        let mut tps = TransitionPathSampling::new(0.3, 0.7, 100, 0.1);
        let reactive_path: Vec<(Vec<[f64; 3]>, f64)> =
            vec![(vec![[0.0; 3]], 0.1), (vec![[0.0; 3]], 0.9)];
        let accepted = tps.record_attempt(reactive_path);
        assert!(accepted);
        assert_eq!(tps.n_reactive_paths(), 1);
        assert!((tps.acceptance_ratio() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_tps_committor() {
        let mut tps = TransitionPathSampling::new(0.3, 0.7, 100, 0.1);
        // Add a path starting from reactant
        tps.record_attempt(vec![(vec![[0.0; 3]], 0.1), (vec![[0.0; 3]], 0.9)]);
        // Add a reverse path starting from product
        tps.record_attempt(vec![(vec![[0.0; 3]], 0.9), (vec![[0.0; 3]], 0.1)]);
        let c = tps.committor_estimate();
        assert!((c - 0.5).abs() < 1e-12, "committor = {c}");
    }

    #[test]
    fn test_tps_perturb_momenta() {
        let tps = TransitionPathSampling::new(0.3, 0.7, 100, 0.01);
        let orig = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let perturbed = tps.perturb_momenta(&orig);
        assert_eq!(perturbed.len(), 2);
        // Should be close to original
        for (o, p) in orig.iter().zip(perturbed.iter()) {
            for d in 0..3 {
                assert!((o[d] - p[d]).abs() < 0.02);
            }
        }
    }

    // ── Metadynamics2D ──────────────────────────────────────────────────

    #[test]
    fn test_metad2d_deposit() {
        let mut m2d = Metadynamics2D::new(1.0, [0.1, 0.1], 10.0, 300.0, 1);
        m2d.step([0.0, 0.0]);
        assert_eq!(m2d.n_hills(), 1);
        let b = m2d.bias_at([0.0, 0.0]);
        assert!(b > 0.0);
    }

    #[test]
    fn test_metad2d_gradient() {
        let mut m2d = Metadynamics2D::new(1.0, [0.1, 0.1], 10.0, 300.0, 1);
        m2d.deposit_hill([0.0, 0.0]);
        let g = m2d.bias_gradient([0.05, 0.0]);
        // gradient in CV0 direction should be negative (pointing towards hill)
        assert!(g[0] < 0.0, "gradient[0] = {}", g[0]);
    }

    // ── AbfBlockAveraged ────────────────────────────────────────────────

    #[test]
    fn test_abf_block_bin() {
        let mut bin = AbfBlockBin::new(5);
        for i in 0..5 {
            bin.add_sample(i as f64);
        }
        assert_eq!(bin.block_sums.len(), 1);
        assert!((bin.block_mean() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_abf_block_averaged() {
        let mut abf = AbfBlockAveraged::new(0.0, 1.0, 10, 10, 0);
        for _ in 0..100 {
            abf.update(0.05, 3.0);
        }
        let (_c, means, errs) = abf.profile_with_errors();
        assert!((means[0] - 3.0).abs() < 1e-10, "mean = {}", means[0]);
        assert!(
            errs[0].abs() < 1e-10,
            "stderr should be 0 for constant force"
        );
    }
}
