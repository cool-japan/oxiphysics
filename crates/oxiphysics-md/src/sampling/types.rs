//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::KB_SAMP;
#[allow(unused_imports)]
use super::functions::*;

/// Nudged elastic band (NEB) for MEP optimisation.
///
/// Adds elastic springs between images to maintain equal spacing,
/// while projecting forces to stay on the MEP.
#[allow(dead_code)]
pub struct NudgedElasticBand {
    /// Images along the band.
    pub images: Vec<PathImage>,
    /// Spring constant for elastic forces between images.
    pub spring_constant: f64,
    /// Step size for gradient descent moves.
    pub step_size: f64,
    /// Convergence threshold (max perpendicular force component).
    pub convergence_tol: f64,
    /// Number of optimisation steps performed.
    pub iterations: usize,
    /// Whether to use climbing image NEB (CI-NEB).
    pub climbing_image: bool,
}
impl NudgedElasticBand {
    /// Create a NEB with linearly interpolated initial images.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: Vec<f64>,
        end: Vec<f64>,
        n_images: usize,
        spring_constant: f64,
        step_size: f64,
        convergence_tol: f64,
        climbing_image: bool,
    ) -> Self {
        assert!(n_images >= 2);
        let dim = start.len();
        assert_eq!(dim, end.len());
        let images: Vec<PathImage> = (0..n_images)
            .map(|i| {
                let t = i as f64 / (n_images - 1) as f64;
                let coords: Vec<f64> = (0..dim)
                    .map(|d| start[d] + t * (end[d] - start[d]))
                    .collect();
                PathImage::new(coords, 0.0)
            })
            .collect();
        Self {
            images,
            spring_constant,
            step_size,
            convergence_tol,
            iterations: 0,
            climbing_image,
        }
    }
    /// Perform one NEB optimisation step.
    ///
    /// `energy_force_fn` computes `(forces, energy)` for given coordinates.
    pub fn step<F>(&mut self, mut energy_force_fn: F) -> f64
    where
        F: FnMut(&[f64]) -> (Vec<f64>, f64),
    {
        let n = self.images.len();
        let dim = self.images[0].coords.len();
        let mut forces: Vec<Vec<f64>> = vec![vec![0.0; dim]; n];
        for i in 0..n {
            let (f, e) = energy_force_fn(&self.images[i].coords);
            forces[i] = f;
            self.images[i].energy = e;
        }
        let ci_idx = if self.climbing_image {
            (1..n - 1)
                .max_by(|&a, &b| {
                    self.images[a]
                        .energy
                        .partial_cmp(&self.images[b].energy)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap_or(0)
        } else {
            usize::MAX
        };
        let mut max_force = 0.0f64;
        for i in 1..n - 1 {
            let e_prev = self.images[i - 1].energy;
            let e_cur = self.images[i].energy;
            let e_next = self.images[i + 1].energy;
            let mut tangent = vec![0.0f64; dim];
            if e_next > e_cur && e_cur > e_prev {
                for d in 0..dim {
                    tangent[d] = self.images[i + 1].coords[d] - self.images[i].coords[d];
                }
            } else if e_next < e_cur && e_cur < e_prev {
                for d in 0..dim {
                    tangent[d] = self.images[i].coords[d] - self.images[i - 1].coords[d];
                }
            } else {
                let de_max = (e_next - e_cur).abs().max((e_prev - e_cur).abs());
                let de_min = (e_next - e_cur).abs().min((e_prev - e_cur).abs());
                for d in 0..dim {
                    let t_plus = self.images[i + 1].coords[d] - self.images[i].coords[d];
                    let t_minus = self.images[i].coords[d] - self.images[i - 1].coords[d];
                    if e_next > e_prev {
                        tangent[d] = t_plus * de_max + t_minus * de_min;
                    } else {
                        tangent[d] = t_plus * de_min + t_minus * de_max;
                    }
                }
            }
            let tan_len: f64 = tangent.iter().map(|&t| t * t).sum::<f64>().sqrt();
            if tan_len > 1e-15 {
                for t in &mut tangent {
                    *t /= tan_len;
                }
            }
            let dist_next: f64 = {
                let mut d2 = 0.0;
                for d in 0..dim {
                    let dd = self.images[i + 1].coords[d] - self.images[i].coords[d];
                    d2 += dd * dd;
                }
                d2.sqrt()
            };
            let dist_prev: f64 = {
                let mut d2 = 0.0;
                for d in 0..dim {
                    let dd = self.images[i].coords[d] - self.images[i - 1].coords[d];
                    d2 += dd * dd;
                }
                d2.sqrt()
            };
            let spring_force_mag = self.spring_constant * (dist_next - dist_prev);
            let total_force: Vec<f64> = if i == ci_idx {
                let f_dot_t: f64 = forces[i]
                    .iter()
                    .zip(tangent.iter())
                    .map(|(f, t)| f * t)
                    .sum();
                (0..dim)
                    .map(|d| forces[i][d] - 2.0 * f_dot_t * tangent[d])
                    .collect()
            } else {
                let f_dot_t: f64 = forces[i]
                    .iter()
                    .zip(tangent.iter())
                    .map(|(f, t)| f * t)
                    .sum();
                (0..dim)
                    .map(|d| {
                        let f_perp = forces[i][d] - f_dot_t * tangent[d];
                        let f_spring = spring_force_mag * tangent[d];
                        f_perp + f_spring
                    })
                    .collect()
            };
            let f_mag: f64 = total_force.iter().map(|&f| f * f).sum::<f64>().sqrt();
            if f_mag > max_force {
                max_force = f_mag;
            }
            self.images[i].perp_force = total_force;
        }
        for i in 1..n - 1 {
            let pf = self.images[i].perp_force.clone();
            for d in 0..dim {
                self.images[i].coords[d] += self.step_size * pf[d];
            }
        }
        self.iterations += 1;
        max_force
    }
    /// Check for convergence (max force < tolerance).
    pub fn is_converged(&self) -> bool {
        let n = self.images.len();
        (1..n - 1).all(|i| {
            let f_mag: f64 = self.images[i]
                .perp_force
                .iter()
                .map(|&f| f * f)
                .sum::<f64>()
                .sqrt();
            f_mag < self.convergence_tol
        })
    }
    /// Number of images.
    pub fn n_images(&self) -> usize {
        self.images.len()
    }
    /// Activation energy (max energy along path minus start energy).
    pub fn activation_energy(&self) -> f64 {
        let e_start = self.images[0].energy;
        self.images
            .iter()
            .map(|im| im.energy - e_start)
            .fold(f64::NEG_INFINITY, f64::max)
    }
}
/// 1D umbrella sampling window.
/// Adds a harmonic bias potential: U_bias = k/2 * (xi - xi0)^2
/// where xi is the collective variable (e.g. distance, dihedral).
#[allow(dead_code)]
pub struct UmbrellaSampling {
    /// Window center (collective variable reference value).
    pub xi0: f64,
    /// Force constant (kJ/mol/nm² or MD units).
    pub k_spring: f64,
    /// Collected xi values.
    pub samples: Vec<f64>,
}
impl UmbrellaSampling {
    /// Create a new umbrella sampling window.
    pub fn new(xi0: f64, k_spring: f64) -> Self {
        Self {
            xi0,
            k_spring,
            samples: Vec::new(),
        }
    }
    /// Compute bias force on collective variable xi.
    /// F_bias = -k * (xi - xi0)
    pub fn bias_force(&self, xi: f64) -> f64 {
        -self.k_spring * (xi - self.xi0)
    }
    /// Add a sample point.
    pub fn collect(&mut self, xi: f64) {
        self.samples.push(xi);
    }
    /// Compute free energy estimate from histogram (Boltzmann inversion).
    /// F(xi) = -kT * ln(P(xi)) + const
    /// Uses bins centered at bin_edges\[i\] .. bin_edges\[i+1\].
    /// Returns (bin_centers, free_energies).
    /// Returns (vec!\[\], vec!\[\]) if no samples collected (avoid log(0)).
    #[allow(non_snake_case)]
    pub fn free_energy_histogram(&self, n_bins: usize, kT: f64) -> (Vec<f64>, Vec<f64>) {
        if self.samples.is_empty() || n_bins == 0 {
            return (vec![], vec![]);
        }
        let xmin = self.samples.iter().cloned().fold(f64::INFINITY, f64::min);
        let xmax = self
            .samples
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        if (xmax - xmin).abs() < 1e-15 {
            let center = xmin;
            return (vec![center], vec![0.0]);
        }
        let bin_width = (xmax - xmin) / n_bins as f64;
        let mut counts = vec![0usize; n_bins];
        for &xi in &self.samples {
            let idx = ((xi - xmin) / bin_width) as usize;
            let idx = idx.min(n_bins - 1);
            counts[idx] += 1;
        }
        let n_total = self.samples.len() as f64;
        let mut bin_centers = Vec::with_capacity(n_bins);
        let mut free_energies = Vec::with_capacity(n_bins);
        let max_count = *counts.iter().max().unwrap_or(&1).max(&1);
        for (i, &cnt) in counts.iter().enumerate() {
            let center = xmin + (i as f64 + 0.5) * bin_width;
            bin_centers.push(center);
            if cnt == 0 {
                free_energies.push(f64::INFINITY);
            } else {
                let f_val = -kT * (cnt as f64 / max_count as f64).ln();
                let _ = n_total;
                free_energies.push(f_val);
            }
        }
        (bin_centers, free_energies)
    }
    /// Mean and variance of collected samples.
    pub fn statistics(&self) -> (f64, f64) {
        let n = self.samples.len();
        if n == 0 {
            return (0.0, 0.0);
        }
        let mean = self.samples.iter().sum::<f64>() / n as f64;
        let variance = self
            .samples
            .iter()
            .map(|&x| (x - mean) * (x - mean))
            .sum::<f64>()
            / n as f64;
        (mean, variance)
    }
}
/// Simulated annealing schedule.
///
/// Provides temperature as a function of step number for use in
/// Monte Carlo or MD simulations.
#[allow(dead_code)]
pub struct SimulatedAnnealing {
    /// Initial temperature.
    pub t_initial: f64,
    /// Final temperature.
    pub t_final: f64,
    /// Total number of steps.
    pub n_steps: usize,
    /// Current step.
    pub current_step: usize,
    /// Schedule type.
    pub schedule: AnnealingSchedule,
}
impl SimulatedAnnealing {
    /// Create a new simulated annealing schedule.
    pub fn new(t_initial: f64, t_final: f64, n_steps: usize, schedule: AnnealingSchedule) -> Self {
        Self {
            t_initial,
            t_final,
            n_steps,
            current_step: 0,
            schedule,
        }
    }
    /// Current temperature at the current step.
    pub fn temperature(&self) -> f64 {
        let frac = if self.n_steps == 0 {
            1.0
        } else {
            (self.current_step as f64 / self.n_steps as f64).min(1.0)
        };
        match self.schedule {
            AnnealingSchedule::Linear => self.t_initial + (self.t_final - self.t_initial) * frac,
            AnnealingSchedule::Exponential => {
                if self.t_initial <= 0.0 || self.t_final <= 0.0 {
                    return self.t_initial * (1.0 - frac) + self.t_final * frac;
                }
                self.t_initial * (self.t_final / self.t_initial).powf(frac)
            }
        }
    }
    /// Advance one step and return the new temperature.
    pub fn step(&mut self) -> f64 {
        if self.current_step < self.n_steps {
            self.current_step += 1;
        }
        self.temperature()
    }
    /// Metropolis acceptance criterion at current temperature.
    /// Returns `true` if `delta_e` should be accepted with random number `u01`.
    pub fn metropolis_accept(&self, delta_e: f64, u01: f64) -> bool {
        if delta_e <= 0.0 {
            return true;
        }
        let t = self.temperature();
        if t <= 0.0 {
            return false;
        }
        u01 < (-delta_e / (KB_SAMP * t)).exp()
    }
    /// Reset the annealing schedule.
    pub fn reset(&mut self) {
        self.current_step = 0;
    }
}
/// Path collective variables (s, z) as defined by Branduardi et al.
///
/// s measures progress along the path (0 to 1).
/// z measures deviation from the path (perpendicular distance).
///
/// Reference: Branduardi, Gervasio, Parrinello, J. Chem. Phys. 126, 054103 (2007).
#[allow(dead_code)]
pub struct PathCV {
    /// Reference path images (coordinates only).
    pub reference_images: Vec<Vec<f64>>,
    /// Smoothing parameter λ (larger = sharper path definition).
    pub lambda: f64,
}
impl PathCV {
    /// Create a new path CV from reference images.
    pub fn new(reference_images: Vec<Vec<f64>>, lambda: f64) -> Self {
        Self {
            reference_images,
            lambda,
        }
    }
    /// Compute path coordinate s and deviation z for a configuration `x`.
    ///
    /// Returns `(s, z)` where s ∈ \[0, 1\] and z ≥ 0.
    pub fn compute(&self, x: &[f64]) -> (f64, f64) {
        let n = self.reference_images.len();
        if n == 0 {
            return (0.0, 0.0);
        }
        let dist_sq: Vec<f64> = self
            .reference_images
            .iter()
            .map(|ref_im| {
                ref_im
                    .iter()
                    .zip(x.iter())
                    .map(|(&r, &xi)| (xi - r) * (xi - r))
                    .sum::<f64>()
            })
            .collect();
        let weights: Vec<f64> = dist_sq.iter().map(|&d| (-self.lambda * d).exp()).collect();
        let w_sum: f64 = weights.iter().sum();
        if w_sum < 1e-300 {
            return (0.0, 0.0);
        }
        let s_num: f64 = weights.iter().enumerate().map(|(i, &w)| i as f64 * w).sum();
        let s = s_num / (w_sum * (n - 1) as f64);
        let z = -(w_sum.ln()) / self.lambda;
        (s, z)
    }
    /// Number of reference images.
    pub fn n_images(&self) -> usize {
        self.reference_images.len()
    }
}
/// Wang-Landau flat-histogram sampling.
///
/// Iteratively estimates the density of states g(E) by performing a random
/// walk in energy space with an adaptive modification factor.
///
/// Reference: Wang & Landau, Phys. Rev. Lett. 86, 2050 (2001).
#[allow(dead_code)]
pub struct WangLandau {
    /// Energy bin edges.
    pub bin_edges: Vec<f64>,
    /// Logarithm of the density of states ln(g(E)) for each bin.
    pub ln_g: Vec<f64>,
    /// Histogram of visits to each bin.
    pub histogram: Vec<usize>,
    /// Current modification factor (starts at 1.0, halved when histogram is flat).
    pub ln_f: f64,
    /// Flatness criterion (fraction of mean, typically 0.8).
    pub flatness: f64,
    /// Minimum modification factor to stop iteration.
    pub ln_f_min: f64,
    /// Total number of MC moves attempted.
    pub total_moves: usize,
}
impl WangLandau {
    /// Create a new Wang-Landau sampler.
    ///
    /// # Arguments
    /// * `e_min` – minimum energy
    /// * `e_max` – maximum energy
    /// * `n_bins` – number of energy bins
    /// * `flatness` – flatness criterion (0.8 is typical)
    pub fn new(e_min: f64, e_max: f64, n_bins: usize, flatness: f64) -> Self {
        let mut bin_edges = Vec::with_capacity(n_bins + 1);
        let de = (e_max - e_min) / n_bins as f64;
        for i in 0..=n_bins {
            bin_edges.push(e_min + i as f64 * de);
        }
        Self {
            bin_edges,
            ln_g: vec![0.0; n_bins],
            histogram: vec![0; n_bins],
            ln_f: 1.0,
            flatness,
            ln_f_min: 1e-8,
            total_moves: 0,
        }
    }
    /// Find the bin index for a given energy. Returns `None` if out of range.
    pub fn bin_index(&self, energy: f64) -> Option<usize> {
        let n_bins = self.ln_g.len();
        if n_bins == 0 {
            return None;
        }
        let e_min = self.bin_edges[0];
        let e_max = self.bin_edges[n_bins];
        if energy < e_min || energy >= e_max {
            return None;
        }
        let de = (e_max - e_min) / n_bins as f64;
        let idx = ((energy - e_min) / de) as usize;
        Some(idx.min(n_bins - 1))
    }
    /// Wang-Landau acceptance criterion for a move from energy `e_old` to `e_new`.
    ///
    /// Accepts with probability min(1, g(e_old)/g(e_new)) = min(1, exp(ln_g_old - ln_g_new)).
    /// `u01` is a uniform random number in \[0, 1).
    pub fn accept_move(&self, e_old: f64, e_new: f64, u01: f64) -> bool {
        let idx_old = match self.bin_index(e_old) {
            Some(i) => i,
            None => return true,
        };
        let idx_new = match self.bin_index(e_new) {
            Some(i) => i,
            None => return false,
        };
        let delta = self.ln_g[idx_old] - self.ln_g[idx_new];
        delta >= 0.0 || u01 < delta.exp()
    }
    /// Update ln(g) and histogram after visiting energy `energy`.
    pub fn visit(&mut self, energy: f64) {
        if let Some(idx) = self.bin_index(energy) {
            self.ln_g[idx] += self.ln_f;
            self.histogram[idx] += 1;
            self.total_moves += 1;
        }
    }
    /// Check if the histogram is flat (all bins ≥ flatness * mean).
    pub fn is_flat(&self) -> bool {
        let n_bins = self.histogram.len();
        if n_bins == 0 {
            return true;
        }
        let total: usize = self.histogram.iter().sum();
        if total == 0 {
            return false;
        }
        let mean = total as f64 / n_bins as f64;
        let threshold = self.flatness * mean;
        self.histogram.iter().all(|&h| h as f64 >= threshold)
    }
    /// Reset the histogram and halve ln_f. Returns `true` if ln_f is still
    /// above `ln_f_min`.
    pub fn update_modification_factor(&mut self) -> bool {
        self.ln_f *= 0.5;
        self.histogram.iter_mut().for_each(|h| *h = 0);
        self.ln_f > self.ln_f_min
    }
    /// Number of energy bins.
    pub fn n_bins(&self) -> usize {
        self.ln_g.len()
    }
    /// Bin centers.
    pub fn bin_centers(&self) -> Vec<f64> {
        let n = self.n_bins();
        (0..n)
            .map(|i| 0.5 * (self.bin_edges[i] + self.bin_edges[i + 1]))
            .collect()
    }
}
/// A transition path: a sequence of states (represented as indices or
/// configurations) connecting two metastable basins.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TransitionPath {
    /// Sequence of states (e.g., collective variable values).
    pub states: Vec<f64>,
    /// Weight (log-weight) associated with this path.
    pub log_weight: f64,
}
impl TransitionPath {
    /// Create a new transition path.
    pub fn new(states: Vec<f64>) -> Self {
        Self {
            states,
            log_weight: 0.0,
        }
    }
    /// Path length (number of frames).
    pub fn len(&self) -> usize {
        self.states.len()
    }
    /// Check if the path is empty.
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }
    /// Check if this path connects basin A (state ≤ `a_boundary`) to
    /// basin B (state ≥ `b_boundary`).
    pub fn connects_basins(&self, a_boundary: f64, b_boundary: f64) -> bool {
        if self.states.is_empty() {
            return false;
        }
        let starts_in_a = self.states[0] <= a_boundary;
        let ends_in_b = self.states[self.states.len() - 1] >= b_boundary;
        let starts_in_b = self.states[0] >= b_boundary;
        let ends_in_a = self.states[self.states.len() - 1] <= a_boundary;
        (starts_in_a && ends_in_b) || (starts_in_b && ends_in_a)
    }
}
/// Transition matrix for Wang-Landau or multicanonical entropy estimation.
///
/// Records accepted MC transitions between energy bins and estimates the
/// density of states via detailed-balance.
#[allow(dead_code)]
pub struct TransitionMatrix {
    /// Number of energy bins.
    pub n_bins: usize,
    /// Transition counts T\[i\]\[j\]: number of moves from bin i to bin j.
    pub counts: Vec<Vec<f64>>,
    /// Current entropy estimate S\[i\] = ln g(E_i).
    pub ln_g: Vec<f64>,
}
impl TransitionMatrix {
    /// Create a new `TransitionMatrix` with `n_bins` energy bins.
    pub fn new(n_bins: usize) -> Self {
        Self {
            n_bins,
            counts: vec![vec![0.0; n_bins]; n_bins],
            ln_g: vec![0.0; n_bins],
        }
    }
    /// Record a transition attempt from bin `from` to bin `to`.
    ///
    /// For rejected moves, `from == to`.
    pub fn record_transition(&mut self, from: usize, to: usize) {
        if from < self.n_bins && to < self.n_bins {
            self.counts[from][to] += 1.0;
        }
    }
    /// Compute the entropy difference ΔS = S(j) - S(i) between two bins
    /// using the detailed-balance relation:
    ///
    /// ```text
    /// T(i→j) / T(j→i) = g(j) / g(i) = exp(S_j - S_i)
    /// ```
    ///
    /// Returns `None` if the reverse transition count is zero.
    pub fn compute_entropy_difference(&self, i: usize, j: usize) -> Option<f64> {
        if i >= self.n_bins || j >= self.n_bins {
            return None;
        }
        let t_ij = self.counts[i][j];
        let t_ji = self.counts[j][i];
        if t_ji < 1e-30 || t_ij < 1e-30 {
            return None;
        }
        Some((t_ij / t_ji).ln())
    }
    /// Update `ln_g` estimates using all available detailed-balance pairs.
    ///
    /// For each adjacent pair (i, i+1) with non-zero counts in both directions,
    /// the entropy difference is estimated and accumulated relative to ln_g\[0\] = 0.
    pub fn update_entropy_estimates(&mut self) {
        self.ln_g[0] = 0.0;
        for i in 0..self.n_bins.saturating_sub(1) {
            let j = i + 1;
            if let Some(ds) = self.compute_entropy_difference(i, j) {
                self.ln_g[j] = self.ln_g[i] + ds;
            }
        }
    }
    /// Total number of recorded transitions (off-diagonal).
    pub fn total_transitions(&self) -> f64 {
        let mut total = 0.0;
        for i in 0..self.n_bins {
            for j in 0..self.n_bins {
                if i != j {
                    total += self.counts[i][j];
                }
            }
        }
        total
    }
}
/// Committor probability estimator.
///
/// The committor p_B(x) is the probability that a trajectory starting
/// from configuration x will reach basin B before basin A.
/// p_B = 0 means "certain to reach A first"; p_B = 1 means "certain to reach B first".
#[allow(dead_code)]
pub struct CommittorAnalysis {
    /// Basin A boundary (collective variable ≤ a_boundary → in A).
    pub a_boundary: f64,
    /// Basin B boundary (collective variable ≥ b_boundary → in B).
    pub b_boundary: f64,
    /// Stored (xi, committor) pairs.
    pub data: Vec<(f64, f64)>,
}
impl CommittorAnalysis {
    /// Create a new committor analysis.
    pub fn new(a_boundary: f64, b_boundary: f64) -> Self {
        Self {
            a_boundary,
            b_boundary,
            data: Vec::new(),
        }
    }
    /// Add a committor estimate for configuration `xi`.
    ///
    /// `p_b` should be in \[0, 1\]: fraction of trial trajectories that reached B.
    pub fn add_point(&mut self, xi: f64, p_b: f64) {
        self.data.push((xi, p_b.clamp(0.0, 1.0)));
    }
    /// Compute committor from trial trajectories.
    ///
    /// `paths` is a list of (starting_xi, ended_in_b) pairs.
    /// Groups by starting xi value (binned) and averages committor.
    pub fn from_trajectories(&mut self, paths: &[(f64, bool)], n_bins: usize) {
        if paths.is_empty() || n_bins == 0 {
            return;
        }
        let xmin = paths
            .iter()
            .map(|&(xi, _)| xi)
            .fold(f64::INFINITY, f64::min);
        let xmax = paths
            .iter()
            .map(|&(xi, _)| xi)
            .fold(f64::NEG_INFINITY, f64::max);
        if (xmax - xmin).abs() < 1e-15 {
            return;
        }
        let bw = (xmax - xmin) / n_bins as f64;
        let mut counts = vec![0usize; n_bins];
        let mut b_counts = vec![0usize; n_bins];
        for &(xi, ended_b) in paths {
            let idx = ((xi - xmin) / bw) as usize;
            let idx = idx.min(n_bins - 1);
            counts[idx] += 1;
            if ended_b {
                b_counts[idx] += 1;
            }
        }
        for i in 0..n_bins {
            if counts[i] > 0 {
                let xi = xmin + (i as f64 + 0.5) * bw;
                let p_b = b_counts[i] as f64 / counts[i] as f64;
                self.add_point(xi, p_b);
            }
        }
    }
    /// Find the transition state xi (where committor ≈ 0.5).
    pub fn transition_state_xi(&self) -> Option<f64> {
        if self.data.is_empty() {
            return None;
        }
        self.data
            .iter()
            .min_by(|(_, p1), (_, p2)| {
                (p1 - 0.5)
                    .abs()
                    .partial_cmp(&(p2 - 0.5).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|&(xi, _)| xi)
    }
    /// Mean committor value.
    pub fn mean_committor(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().map(|&(_, p)| p).sum::<f64>() / self.data.len() as f64
    }
    /// Number of data points.
    pub fn n_points(&self) -> usize {
        self.data.len()
    }
}
/// Multicanonical (MUCA) simulation state.
///
/// Uses iteratively updated weights W(E) = 1/g(E) to achieve flat energy
/// sampling across the entire energy range.
#[allow(dead_code)]
pub struct MulticanonicalSimulation {
    /// Energy bin edges (n_bins + 1 values).
    pub bin_edges: Vec<f64>,
    /// Log density of states ln g(E) for each bin.
    pub ln_dos: Vec<f64>,
    /// Histogram H(E): visit counts for each bin.
    pub histogram: Vec<u64>,
}
impl MulticanonicalSimulation {
    /// Create a new `MulticanonicalSimulation` over \[`e_min`, `e_max`\] with `n_bins` bins.
    pub fn new(e_min: f64, e_max: f64, n_bins: usize) -> Self {
        let de = if n_bins > 0 {
            (e_max - e_min) / n_bins as f64
        } else {
            1.0
        };
        let bin_edges: Vec<f64> = (0..=n_bins).map(|i| e_min + i as f64 * de).collect();
        Self {
            bin_edges,
            ln_dos: vec![0.0; n_bins],
            histogram: vec![0; n_bins],
        }
    }
    /// Find the bin index for energy `e`. Returns `None` if out of range.
    pub fn bin_index(&self, e: f64) -> Option<usize> {
        let n = self.ln_dos.len();
        if n == 0 {
            return None;
        }
        let e_min = self.bin_edges[0];
        let e_max = self.bin_edges[self.bin_edges.len() - 1];
        if e < e_min || e >= e_max {
            return None;
        }
        let de = (e_max - e_min) / n as f64;
        let idx = ((e - e_min) / de) as usize;
        Some(idx.min(n - 1))
    }
    /// Record a visit to energy `e` and update the histogram.
    pub fn visit(&mut self, e: f64) {
        if let Some(idx) = self.bin_index(e) {
            self.histogram[idx] += 1;
        }
    }
    /// Compute/update the density of states from the flat-histogram visit counts.
    ///
    /// The density of states estimate is updated as:
    ///
    /// ```text
    /// ln g(E) += ln H(E) - ln `H`
    /// ```
    ///
    /// This shifts the DoS so that the most visited state has ln g = 0 offset.
    /// After the update the histogram is reset for the next iteration.
    pub fn compute_density_of_states(&mut self) {
        let n = self.ln_dos.len();
        if n == 0 {
            return;
        }
        let total: u64 = self.histogram.iter().sum();
        if total == 0 {
            return;
        }
        let mean = total as f64 / n as f64;
        for (lng, &h) in self.ln_dos.iter_mut().zip(self.histogram.iter()) {
            if h > 0 {
                *lng += (h as f64 / mean).ln();
            }
        }
        let offset = self.ln_dos[0];
        for v in self.ln_dos.iter_mut() {
            *v -= offset;
        }
        self.histogram.iter_mut().for_each(|h| *h = 0);
    }
    /// Multicanonical weight for energy `e`: W(E) = exp(-ln_dos(E)).
    ///
    /// A move is accepted with probability min(1, W(E_new)/W(E_old)).
    pub fn weight(&self, e: f64) -> f64 {
        match self.bin_index(e) {
            Some(idx) => (-self.ln_dos[idx]).exp(),
            None => 0.0,
        }
    }
    /// Number of energy bins.
    pub fn n_bins(&self) -> usize {
        self.ln_dos.len()
    }
}
/// Cooling schedule types for simulated annealing.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub enum AnnealingSchedule {
    /// Linear cooling: T(t) = T_i + (T_f - T_i) * t / n_steps
    Linear,
    /// Exponential cooling: T(t) = T_i * (T_f / T_i)^(t / n_steps)
    Exponential,
}
/// Replica exchange (parallel tempering) manager.
///
/// Manages N replicas at different temperatures and proposes swap moves
/// between adjacent replicas using the Metropolis criterion:
///
/// P_swap = min(1, exp(Δβ · ΔE))
///
/// where Δβ = β_j - β_i and ΔE = E_j - E_i for adjacent replicas i, j.
#[allow(dead_code)]
pub struct ReplicaExchange {
    /// Temperatures for each replica (K).
    pub temperatures: Vec<f64>,
    /// Current energies for each replica.
    pub energies: Vec<f64>,
    /// Total accepted swaps.
    pub accepted: usize,
    /// Total attempted swaps.
    pub attempted: usize,
}
impl ReplicaExchange {
    /// Create a new replica exchange manager.
    ///
    /// `temperatures` must be sorted in ascending order.
    pub fn new(temperatures: Vec<f64>) -> Self {
        let n = temperatures.len();
        Self {
            temperatures,
            energies: vec![0.0; n],
            accepted: 0,
            attempted: 0,
        }
    }
    /// Set the energy of replica `idx`.
    pub fn set_energy(&mut self, idx: usize, energy: f64) {
        self.energies[idx] = energy;
    }
    /// Attempt a swap between replicas `i` and `j`.
    ///
    /// Returns `true` if the swap was accepted (using a deterministic
    /// Metropolis criterion with a provided uniform random number `u01`
    /// in \[0, 1)).
    pub fn try_swap(&mut self, i: usize, j: usize, u01: f64) -> bool {
        self.attempted += 1;
        let beta_i = 1.0 / (KB_SAMP * self.temperatures[i]);
        let beta_j = 1.0 / (KB_SAMP * self.temperatures[j]);
        let delta = (beta_j - beta_i) * (self.energies[i] - self.energies[j]);
        let accepted = delta <= 0.0 || u01 < delta.exp();
        if accepted {
            self.accepted += 1;
            self.energies.swap(i, j);
        }
        accepted
    }
    /// Acceptance rate.
    pub fn acceptance_rate(&self) -> f64 {
        if self.attempted == 0 {
            0.0
        } else {
            self.accepted as f64 / self.attempted as f64
        }
    }
    /// Number of replicas.
    pub fn n_replicas(&self) -> usize {
        self.temperatures.len()
    }
}
/// A single image (bead) on the string / NEB path.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PathImage {
    /// Collective variable values for this image.
    pub coords: Vec<f64>,
    /// Potential energy at this image.
    pub energy: f64,
    /// Forces projected out of the path direction.
    pub perp_force: Vec<f64>,
}
impl PathImage {
    /// Create a new path image.
    pub fn new(coords: Vec<f64>, energy: f64) -> Self {
        let n = coords.len();
        Self {
            coords,
            energy,
            perp_force: vec![0.0; n],
        }
    }
}
/// Multi-window WHAM (Weighted Histogram Analysis Method) — simplified.
/// Given N umbrella windows, compute the unbiased free energy profile.
#[allow(dead_code, non_snake_case)]
pub struct WhamAnalysis {
    /// Collection of umbrella sampling windows.
    pub windows: Vec<UmbrellaSampling>,
    /// Thermal energy kT (same units as k_spring * xi²).
    pub kT: f64,
}
impl WhamAnalysis {
    /// Create a new WHAM analysis with given thermal energy kT.
    #[allow(non_snake_case)]
    pub fn new(kT: f64) -> Self {
        Self {
            windows: Vec::new(),
            kT,
        }
    }
    /// Add an umbrella sampling window.
    pub fn add_window(&mut self, window: UmbrellaSampling) {
        self.windows.push(window);
    }
    /// Simple WHAM iteration to find unbiased free energies.
    /// Returns (xi_values, F_values) for the combined PMF.
    pub fn compute_pmf(&self, n_bins: usize, max_iter: usize) -> (Vec<f64>, Vec<f64>) {
        if self.windows.is_empty() || n_bins == 0 {
            return (vec![], vec![]);
        }
        let all_samples: Vec<f64> = self
            .windows
            .iter()
            .flat_map(|w| w.samples.iter().cloned())
            .collect();
        if all_samples.is_empty() {
            return (vec![], vec![]);
        }
        let xmin = all_samples.iter().cloned().fold(f64::INFINITY, f64::min);
        let xmax = all_samples
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        if (xmax - xmin).abs() < 1e-15 {
            return (vec![xmin], vec![0.0]);
        }
        let bin_width = (xmax - xmin) / n_bins as f64;
        let bin_centers: Vec<f64> = (0..n_bins)
            .map(|i| xmin + (i as f64 + 0.5) * bin_width)
            .collect();
        let n_windows = self.windows.len();
        let mut hist = vec![vec![0usize; n_bins]; n_windows];
        let mut n_i = vec![0usize; n_windows];
        for (wi, window) in self.windows.iter().enumerate() {
            for &xi in &window.samples {
                let idx = ((xi - xmin) / bin_width) as usize;
                let idx = idx.min(n_bins - 1);
                hist[wi][idx] += 1;
            }
            n_i[wi] = window.samples.len();
        }
        let mut f_i = vec![0.0f64; n_windows];
        let bias: Vec<Vec<f64>> = self
            .windows
            .iter()
            .map(|w| {
                bin_centers
                    .iter()
                    .map(|&xi| 0.5 * w.k_spring * (xi - w.xi0) * (xi - w.xi0))
                    .collect()
            })
            .collect();
        let mut rho = vec![0.0f64; n_bins];
        for _iter in 0..max_iter {
            for b in 0..n_bins {
                let numerator: f64 = (0..n_windows).map(|i| hist[i][b] as f64).sum();
                let denominator: f64 = (0..n_windows)
                    .map(|i| n_i[i] as f64 * (f_i[i] / self.kT - bias[i][b] / self.kT).exp())
                    .sum();
                rho[b] = if denominator > 0.0 {
                    numerator / denominator
                } else {
                    0.0
                };
            }
            let f_i_new: Vec<f64> = (0..n_windows)
                .map(|i| {
                    let s: f64 = (0..n_bins)
                        .map(|b| rho[b] * (-bias[i][b] / self.kT).exp())
                        .sum();
                    if s > 0.0 { -self.kT * s.ln() } else { 0.0 }
                })
                .collect();
            let delta: f64 = f_i_new
                .iter()
                .zip(f_i.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            f_i = f_i_new;
            if delta < 1e-8 {
                break;
            }
        }
        let max_rho = rho.iter().cloned().fold(0.0_f64, f64::max);
        let free_energies: Vec<f64> = rho
            .iter()
            .map(|&r| {
                if r <= 0.0 || max_rho <= 0.0 {
                    f64::INFINITY
                } else {
                    -self.kT * (r / max_rho).ln()
                }
            })
            .collect();
        (bin_centers, free_energies)
    }
}
/// String method for minimum energy path (MEP) finding.
///
/// Maintains a chain of images between two endpoints and relaxes
/// them onto the MEP while keeping equal arc-length spacing.
#[allow(dead_code)]
pub struct StringMethod {
    /// Images along the path.
    pub images: Vec<PathImage>,
    /// Step size for image relaxation.
    pub step_size: f64,
    /// Number of relaxation iterations performed.
    pub iterations: usize,
}
impl StringMethod {
    /// Create a string method with linearly interpolated initial images.
    ///
    /// `start` and `end` are the endpoint coordinates.
    /// `n_images` is the total number of images (including endpoints).
    pub fn new(start: Vec<f64>, end: Vec<f64>, n_images: usize, step_size: f64) -> Self {
        assert!(n_images >= 2, "need at least 2 images");
        let dim = start.len();
        assert_eq!(dim, end.len(), "start/end must have same dimension");
        let images: Vec<PathImage> = (0..n_images)
            .map(|i| {
                let t = i as f64 / (n_images - 1) as f64;
                let coords: Vec<f64> = (0..dim)
                    .map(|d| start[d] + t * (end[d] - start[d]))
                    .collect();
                PathImage::new(coords, 0.0)
            })
            .collect();
        Self {
            images,
            step_size,
            iterations: 0,
        }
    }
    /// Re-parameterise the string to equal arc-length spacing.
    ///
    /// Computes cumulative arc length and redistributes images by
    /// linear interpolation along the path.
    pub fn reparameterise(&mut self) {
        let n = self.images.len();
        if n < 3 {
            return;
        }
        let dim = self.images[0].coords.len();
        let mut arc = vec![0.0f64; n];
        for i in 1..n {
            let mut ds = 0.0;
            for d in 0..dim {
                let dd = self.images[i].coords[d] - self.images[i - 1].coords[d];
                ds += dd * dd;
            }
            arc[i] = arc[i - 1] + ds.sqrt();
        }
        let total_len = arc[n - 1];
        if total_len < 1e-15 {
            return;
        }
        let old_coords: Vec<Vec<f64>> = self.images.iter().map(|im| im.coords.clone()).collect();
        for i in 1..n - 1 {
            let target = total_len * i as f64 / (n - 1) as f64;
            let seg = arc
                .windows(2)
                .position(|w| w[0] <= target && target <= w[1])
                .unwrap_or(n - 2);
            let dlen = arc[seg + 1] - arc[seg];
            let t = if dlen < 1e-15 {
                0.0
            } else {
                (target - arc[seg]) / dlen
            };
            for d in 0..dim {
                self.images[i].coords[d] =
                    old_coords[seg][d] + t * (old_coords[seg + 1][d] - old_coords[seg][d]);
            }
        }
    }
    /// Apply one gradient descent step to interior images.
    ///
    /// `force_fn` computes the (negative gradient of energy, energy) for given coords.
    pub fn step<F>(&mut self, mut force_fn: F)
    where
        F: FnMut(&[f64]) -> (Vec<f64>, f64),
    {
        let n = self.images.len();
        let dim = self.images[0].coords.len();
        for i in 1..n - 1 {
            let (forces, energy) = force_fn(&self.images[i].coords);
            self.images[i].energy = energy;
            let mut tangent = vec![0.0f64; dim];
            for d in 0..dim {
                tangent[d] = self.images[i + 1].coords[d] - self.images[i - 1].coords[d];
            }
            let tan_len: f64 = tangent.iter().map(|&t| t * t).sum::<f64>().sqrt();
            if tan_len > 1e-15 {
                for t in &mut tangent {
                    *t /= tan_len;
                }
            }
            let f_dot_t: f64 = forces.iter().zip(tangent.iter()).map(|(f, t)| f * t).sum();
            for d in 0..dim {
                self.images[i].perp_force[d] = forces[d] - f_dot_t * tangent[d];
            }
        }
        for i in 1..n - 1 {
            let pf = self.images[i].perp_force.clone();
            for d in 0..dim {
                self.images[i].coords[d] += self.step_size * pf[d];
            }
        }
        self.iterations += 1;
    }
    /// Path length (total arc length of current images).
    pub fn path_length(&self) -> f64 {
        let dim = self.images[0].coords.len();
        self.images
            .windows(2)
            .map(|w| {
                let mut ds = 0.0;
                for d in 0..dim {
                    let dd = w[1].coords[d] - w[0].coords[d];
                    ds += dd * dd;
                }
                ds.sqrt()
            })
            .sum()
    }
    /// Find the image with maximum energy (approximate saddle point).
    pub fn transition_state_index(&self) -> usize {
        self.images
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.energy
                    .partial_cmp(&b.energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}
/// Nested sampling state for computing the Bayesian evidence (marginal likelihood).
///
/// Nested sampling (Skilling 2006) works by maintaining a set of "live points"
/// with likelihood > some threshold and iteratively replacing the lowest-likelihood
/// point with a new sample drawn from the prior with higher likelihood.  The
/// evidence integral Z = ∫ L(θ) π(θ) dθ is accumulated as a Riemann-like sum.
#[allow(dead_code)]
pub struct NestedSampling {
    /// Number of live points.
    pub n_live: usize,
    /// Log-likelihood values of all dead points (in order of removal).
    pub dead_log_likelihoods: Vec<f64>,
    /// Log-prior-volume fraction for each dead point.
    pub log_prior_volumes: Vec<f64>,
    /// Accumulated log-evidence.
    pub log_evidence: f64,
}
impl NestedSampling {
    /// Create a new nested sampler with `n_live` live points.
    pub fn new(n_live: usize) -> Self {
        Self {
            n_live,
            dead_log_likelihoods: Vec::new(),
            log_prior_volumes: Vec::new(),
            log_evidence: f64::NEG_INFINITY,
        }
    }
    /// Record a dead point with given log-likelihood.
    ///
    /// Each time the lowest-likelihood live point is removed, the prior-volume
    /// fraction shrinks by factor exp(-1/N_live).  This method records one step.
    ///
    /// # Arguments
    /// * `log_like` – log-likelihood of the removed point.
    pub fn record_dead_point(&mut self, log_like: f64) {
        let iteration = self.dead_log_likelihoods.len();
        let log_x = -(iteration as f64 + 1.0) / self.n_live as f64;
        self.dead_log_likelihoods.push(log_like);
        self.log_prior_volumes.push(log_x);
    }
    /// Compute the Bayesian evidence estimate (log Z).
    ///
    /// Uses the trapezoidal rule over the prior-volume axis:
    ///
    /// ```text
    /// log Z ≈ log Σ_i L_i * w_i
    /// ```
    ///
    /// where w_i = ½ (X_{i-1} - X_{i+1}) (widths) and X_i = exp(log_X_i).
    /// Returns `f64::NEG_INFINITY` if no dead points have been recorded.
    pub fn compute_evidence(&self) -> f64 {
        let n = self.dead_log_likelihoods.len();
        if n == 0 {
            return f64::NEG_INFINITY;
        }
        let mut x_vals: Vec<f64> = std::iter::once(1.0_f64)
            .chain(self.log_prior_volumes.iter().map(|&lx| lx.exp()))
            .chain(std::iter::once(0.0_f64))
            .collect();
        let mut log_terms: Vec<f64> = Vec::with_capacity(n);
        for i in 0..n {
            let w_i = 0.5 * (x_vals[i] - x_vals[i + 2]);
            if w_i > 0.0 {
                let log_term = self.dead_log_likelihoods[i] + w_i.ln();
                log_terms.push(log_term);
            }
        }
        if log_terms.is_empty() {
            return f64::NEG_INFINITY;
        }
        let max_term = log_terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = log_terms.iter().map(|&v| (v - max_term).exp()).sum();
        let _ = x_vals.pop();
        max_term + sum.ln()
    }
    /// Number of dead points accumulated.
    pub fn n_dead(&self) -> usize {
        self.dead_log_likelihoods.len()
    }
}
/// Transition path sampling (TPS) engine.
///
/// Stores an ensemble of reactive paths and provides shooting-move
/// proposal logic.
#[allow(dead_code)]
pub struct TransitionPathSampling {
    /// Collection of accepted paths.
    pub paths: Vec<TransitionPath>,
    /// Basin A boundary (collective variable value).
    pub a_boundary: f64,
    /// Basin B boundary (collective variable value).
    pub b_boundary: f64,
    /// Number of accepted shooting moves.
    pub accepted_shoots: usize,
    /// Number of attempted shooting moves.
    pub attempted_shoots: usize,
}
impl TransitionPathSampling {
    /// Create a new TPS engine.
    pub fn new(a_boundary: f64, b_boundary: f64) -> Self {
        Self {
            paths: Vec::new(),
            a_boundary,
            b_boundary,
            accepted_shoots: 0,
            attempted_shoots: 0,
        }
    }
    /// Add a path to the ensemble if it connects the two basins.
    /// Returns `true` if the path was accepted.
    pub fn add_path(&mut self, path: TransitionPath) -> bool {
        if path.connects_basins(self.a_boundary, self.b_boundary) {
            self.paths.push(path);
            true
        } else {
            false
        }
    }
    /// Propose a shooting move: pick a random time slice from the current
    /// path and perturb it. The caller provides:
    /// - `path_idx`: index of the current path
    /// - `shoot_idx`: time slice to perturb
    /// - `new_path`: the newly generated path after perturbation
    ///
    /// Returns `true` if the new path was accepted (connects basins).
    pub fn shooting_move(&mut self, path_idx: usize, new_path: TransitionPath) -> bool {
        self.attempted_shoots += 1;
        if new_path.connects_basins(self.a_boundary, self.b_boundary) {
            self.accepted_shoots += 1;
            self.paths[path_idx] = new_path;
            true
        } else {
            false
        }
    }
    /// Acceptance rate for shooting moves.
    pub fn acceptance_rate(&self) -> f64 {
        if self.attempted_shoots == 0 {
            0.0
        } else {
            self.accepted_shoots as f64 / self.attempted_shoots as f64
        }
    }
}
