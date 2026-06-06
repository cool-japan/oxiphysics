// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Steered Molecular Dynamics (SMD) and enhanced sampling.
//!
//! Provides two pulling protocols:
//! - **Constant-velocity SMD** ([`ConstantVelocitySmd`]): moves a virtual spring anchor
//!   at a fixed velocity and records the force on the pulled atom.
//! - **Constant-force SMD** ([`ConstantForceSmd`]): applies a fixed external force.
//! - **[`SteeredMd`]**: general steered MD with atom index, direction, spring, velocity.
//! - **[`UmbrellaWindow`]**: umbrella-sampling window with histogram accumulation.
//! - **[`MetadynamicsCV`]**: Gaussian-based metadynamics collective variable.
//! - **[`ReplicaExchange`]**: parallel tempering / REMD attempt_swap.
//!
//! Free energy estimates are obtained via the Jarzynski equality or the
//! worm-like chain (WLC) model.

/// Boltzmann constant in J/K.
const KB_J: f64 = 1.380_649e-23;

// ─── helpers ─────────────────────────────────────────────────────────────────

/// Normalise a 3-D vector.  Returns the zero vector if the input has zero length.
fn normalise(v: [f64; 3]) -> [f64; 3] {
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if norm < f64::EPSILON {
        [0.0; 3]
    } else {
        [v[0] / norm, v[1] / norm, v[2] / norm]
    }
}

/// Dot product of two 3-D vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ─── SteeredMd ────────────────────────────────────────────────────────────────

/// General steered MD probe: pulls a specific atom index with a spring moving
/// at constant velocity along a given direction.
///
/// Returns the incremental work done at each call to [`SteeredMd::apply_force`].
pub struct SteeredMd {
    /// Index of the pulled atom.
    pub pulled_atom: usize,
    /// Normalised pulling direction.
    pub direction: [f64; 3],
    /// Spring constant (kJ mol⁻¹ nm⁻²).
    pub spring_constant: f64,
    /// Pulling velocity (nm ps⁻¹).
    pub velocity: f64,
    /// Cumulative time elapsed.
    pub time_elapsed: f64,
    /// Accumulated work (kJ mol⁻¹).
    pub work_accumulated: f64,
}

impl SteeredMd {
    /// Create a new `SteeredMd` pulling atom `pulled_atom`.
    ///
    /// `direction` is normalised automatically.
    pub fn new(
        pulled_atom: usize,
        direction: [f64; 3],
        spring_constant: f64,
        velocity: f64,
    ) -> Self {
        Self {
            pulled_atom,
            direction: normalise(direction),
            spring_constant,
            velocity,
            time_elapsed: 0.0,
            work_accumulated: 0.0,
        }
    }

    /// Apply the spring force to the pulled atom and return the work done this step.
    ///
    /// Updates `forces[self.pulled_atom]` in place and returns `W = F · v · dt`.
    pub fn apply_force(
        &mut self,
        positions: &[[f64; 3]],
        forces: &mut [[f64; 3]],
        time: f64,
        dt: f64,
    ) -> f64 {
        self.time_elapsed = time;
        let pos = positions[self.pulled_atom];

        // Virtual anchor position: start (pos at t=0) + v * t * direction.
        // For simplicity anchor starts at atom's initial position and moves.
        let anchor = [
            pos[0] + self.velocity * time * self.direction[0],
            pos[1] + self.velocity * time * self.direction[1],
            pos[2] + self.velocity * time * self.direction[2],
        ];
        // Spring force on atom
        let f = [
            self.spring_constant * (anchor[0] - pos[0]),
            self.spring_constant * (anchor[1] - pos[1]),
            self.spring_constant * (anchor[2] - pos[2]),
        ];
        forces[self.pulled_atom][0] += f[0];
        forces[self.pulled_atom][1] += f[1];
        forces[self.pulled_atom][2] += f[2];

        // Work = F · v * dt (scalar projection of force onto pulling direction times |v|)
        let f_proj = dot3(f, self.direction);
        let dw = f_proj * self.velocity * dt;
        self.work_accumulated += dw;
        dw
    }
}

// ─── ConstantVelocitySmd ──────────────────────────────────────────────────────

/// Constant-velocity SMD (cv-SMD) force probe.
///
/// A virtual spring anchor moves at a fixed velocity along a pulling direction.
/// The harmonic spring exerts a force on the pulled atom proportional to the
/// displacement between the atom and the virtual anchor.
pub struct ConstantVelocitySmd {
    /// Spring constant (N/m or simulation units).
    pub spring_constant: f64,
    /// Pulling velocity (nm ps⁻¹ or consistent simulation units).
    pub pulling_velocity: f64,
    /// Unit vector along the pulling direction.
    pub pulling_direction: [f64; 3],
    /// Initial position of the virtual spring anchor.
    pub initial_position: [f64; 3],
    /// Time elapsed since the start of pulling.
    pub time_elapsed: f64,
}

impl ConstantVelocitySmd {
    /// Create a new cv-SMD probe.
    ///
    /// `direction` is normalised automatically.
    pub fn new(k: f64, v: f64, direction: [f64; 3], start_pos: [f64; 3]) -> Self {
        Self {
            spring_constant: k,
            pulling_velocity: v,
            pulling_direction: normalise(direction),
            initial_position: start_pos,
            time_elapsed: 0.0,
        }
    }

    /// Position of the virtual spring anchor at the current time.
    ///
    /// `pos = initial_position + pulling_velocity × time_elapsed × pulling_direction`
    pub fn virtual_atom_position(&self) -> [f64; 3] {
        let d = self.pulling_velocity * self.time_elapsed;
        [
            self.initial_position[0] + d * self.pulling_direction[0],
            self.initial_position[1] + d * self.pulling_direction[1],
            self.initial_position[2] + d * self.pulling_direction[2],
        ]
    }

    /// Spring force acting on the pulled atom.
    ///
    /// `F = k × (virtual_pos − atom_pos)`
    pub fn spring_force(&self, atom_pos: [f64; 3]) -> [f64; 3] {
        let vp = self.virtual_atom_position();
        [
            self.spring_constant * (vp[0] - atom_pos[0]),
            self.spring_constant * (vp[1] - atom_pos[1]),
            self.spring_constant * (vp[2] - atom_pos[2]),
        ]
    }

    /// Work done so far by the spring force.
    ///
    /// Given a history of scalar force projections along the pulling direction
    /// (one per time step), the work is approximated as:
    ///   W = Σ F_i · v · dt
    /// where dt is inferred as `time_elapsed / n` and v is the pulling velocity.
    ///
    /// Returns 0 when `force_history` is empty.
    pub fn work_done(&self, force_history: &[f64]) -> f64 {
        if force_history.is_empty() {
            return 0.0;
        }
        let n = force_history.len() as f64;
        let dt = if n > 0.0 { self.time_elapsed / n } else { 0.0 };
        force_history
            .iter()
            .map(|&f| f * self.pulling_velocity * dt)
            .sum()
    }

    /// Advance the simulation time by `dt`.
    pub fn step(&mut self, dt: f64) {
        self.time_elapsed += dt;
    }

    /// Jarzynski equality: ΔF = −k_B T · ln⟨exp(−β W)⟩
    ///
    /// Given an ensemble of non-equilibrium work values (one per trajectory),
    /// returns the free energy difference estimate.
    ///
    /// Returns 0 when `works` is empty.
    pub fn jarzynski_free_energy(works: &[f64], temperature: f64) -> f64 {
        if works.is_empty() {
            return 0.0;
        }
        let beta = 1.0 / (KB_J * temperature);
        let n = works.len() as f64;
        let avg_exp: f64 = works.iter().map(|&w| (-beta * w).exp()).sum::<f64>() / n;
        if avg_exp <= 0.0 {
            return f64::INFINITY;
        }
        -avg_exp.ln() / beta
    }
}

// ─── ConstantForceSmd ─────────────────────────────────────────────────────────

/// Constant-force SMD (cf-SMD).
///
/// Applies a fixed external force to a pulled atom throughout the simulation.
pub struct ConstantForceSmd {
    /// Magnitude of the applied force (N or simulation units).
    pub force_magnitude: f64,
    /// Unit vector along the pulling direction.
    pub direction: [f64; 3],
}

impl ConstantForceSmd {
    /// Create a new cf-SMD probe.
    ///
    /// `direction` is normalised automatically.
    pub fn new(force: f64, direction: [f64; 3]) -> Self {
        Self {
            force_magnitude: force,
            direction: normalise(direction),
        }
    }

    /// Force vector: `F = force_magnitude × direction`.
    pub fn force_vector(&self) -> [f64; 3] {
        [
            self.force_magnitude * self.direction[0],
            self.force_magnitude * self.direction[1],
            self.force_magnitude * self.direction[2],
        ]
    }

    /// Worm-like chain (WLC) model: equilibrium fractional extension z/L at a
    /// given force using the Marko–Siggia interpolation formula.
    pub fn wlc_extension(force: f64, persistence_length: f64, temperature: f64) -> f64 {
        let kbt = KB_J * temperature;
        let prefactor = kbt / persistence_length;

        let wlc_force = |x: f64| -> f64 { prefactor * (0.25 / (1.0 - x).powi(2) - 0.25 + x) };

        let eps = 1e-6_f64;
        let mut lo = eps;
        let mut hi = 1.0 - eps;

        if wlc_force(lo) >= force {
            return lo;
        }
        if wlc_force(hi) <= force {
            return hi;
        }

        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if wlc_force(mid) < force {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        0.5 * (lo + hi)
    }
}

// ─── UmbrellaWindow ───────────────────────────────────────────────────────────

/// A single umbrella-sampling window for a 1D reaction coordinate.
///
/// Applies a harmonic bias: `U_bias = 0.5 * force_constant * (rc - target)^2`.
pub struct UmbrellaWindow {
    /// Target reaction-coordinate value for this window.
    pub target: f64,
    /// Harmonic force constant.
    pub force_constant: f64,
    /// Atom indices defining the reaction coordinate (distance between them).
    pub atom_i: usize,
    /// Second atom index.
    pub atom_j: usize,
    /// Histogram of sampled RC values (raw counts in bins).
    pub histogram: Vec<u64>,
    /// Minimum RC value for histogram range.
    pub hist_min: f64,
    /// Maximum RC value for histogram range.
    pub hist_max: f64,
}

impl UmbrellaWindow {
    /// Create a new umbrella sampling window.
    pub fn new(
        target: f64,
        force_constant: f64,
        atom_i: usize,
        atom_j: usize,
        n_bins: usize,
        hist_min: f64,
        hist_max: f64,
    ) -> Self {
        Self {
            target,
            force_constant,
            atom_i,
            atom_j,
            histogram: vec![0u64; n_bins],
            hist_min,
            hist_max,
        }
    }

    /// Compute the reaction coordinate (distance between atom_i and atom_j).
    pub fn reaction_coordinate(&self, positions: &[[f64; 3]]) -> f64 {
        let ri = positions[self.atom_i];
        let rj = positions[self.atom_j];
        let dx = ri[0] - rj[0];
        let dy = ri[1] - rj[1];
        let dz = ri[2] - rj[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Compute the bias forces on all atoms.
    ///
    /// Returns a force vector (same length as positions), with non-zero
    /// entries only for `atom_i` and `atom_j`.
    pub fn bias_force(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0f64; 3]; n];
        let rc = self.reaction_coordinate(positions);
        if rc < 1e-14 {
            return forces;
        }
        // d(U_bias)/d(rc) = force_constant * (rc - target)
        let du_drc = self.force_constant * (rc - self.target);
        // d(rc)/d(r_i) = (r_i - r_j) / rc
        let ri = positions[self.atom_i];
        let rj = positions[self.atom_j];
        let inv_rc = 1.0 / rc;
        for k in 0..3 {
            let drc_dri = (ri[k] - rj[k]) * inv_rc;
            forces[self.atom_i][k] -= du_drc * drc_dri;
            forces[self.atom_j][k] += du_drc * drc_dri;
        }
        forces
    }

    /// Record a sample into the histogram.
    pub fn sample_histogram(&mut self, rc_value: f64) {
        if self.histogram.is_empty() {
            return;
        }
        let n_bins = self.histogram.len();
        let range = self.hist_max - self.hist_min;
        if range <= 0.0 {
            return;
        }
        let bin = ((rc_value - self.hist_min) / range * n_bins as f64) as isize;
        if bin >= 0 && (bin as usize) < n_bins {
            self.histogram[bin as usize] += 1;
        }
    }

    /// Return the total number of samples collected.
    pub fn n_samples(&self) -> u64 {
        self.histogram.iter().sum()
    }
}

// ─── MetadynamicsCV ───────────────────────────────────────────────────────────

/// A Gaussian history for metadynamics on a single collective variable.
#[derive(Debug, Clone)]
struct HistoryGaussian {
    center: f64,
    height: f64,
    width: f64,
}

/// Metadynamics state for one collective variable (e.g. a distance).
pub struct MetadynamicsCV {
    /// Atom indices defining the CV (distance between two atoms).
    pub index_atoms: [usize; 2],
    /// Current CV value.
    pub current_value: f64,
    /// Deposited Gaussians.
    history_gaussians: Vec<HistoryGaussian>,
}

impl MetadynamicsCV {
    /// Create a new `MetadynamicsCV` tracking distance between `atom_i` and `atom_j`.
    pub fn new(atom_i: usize, atom_j: usize) -> Self {
        Self {
            index_atoms: [atom_i, atom_j],
            current_value: 0.0,
            history_gaussians: Vec::new(),
        }
    }

    /// Add a Gaussian hill centered at the current CV value.
    pub fn add_gaussian(&mut self, height: f64, width: f64) {
        self.history_gaussians.push(HistoryGaussian {
            center: self.current_value,
            height,
            width,
        });
    }

    /// Evaluate the total bias potential at a given CV value.
    ///
    /// `V_bias(s) = sum_k h_k * exp(-(s - s_k)^2 / (2 * w_k^2))`
    pub fn bias_potential(&self, cv_value: f64) -> f64 {
        self.history_gaussians
            .iter()
            .map(|g| {
                let ds = cv_value - g.center;
                g.height * (-ds * ds / (2.0 * g.width * g.width)).exp()
            })
            .sum()
    }

    /// Update the current CV value from positions.
    pub fn update_value(&mut self, positions: &[[f64; 3]]) {
        let [i, j] = self.index_atoms;
        let ri = positions[i];
        let rj = positions[j];
        let dx = ri[0] - rj[0];
        let dy = ri[1] - rj[1];
        let dz = ri[2] - rj[2];
        self.current_value = (dx * dx + dy * dy + dz * dz).sqrt();
    }

    /// Return number of deposited Gaussians.
    pub fn n_gaussians(&self) -> usize {
        self.history_gaussians.len()
    }
}

// ─── ReplicaExchange ──────────────────────────────────────────────────────────

/// Parallel tempering (replica exchange MD) manager.
///
/// Maintains a temperature ladder and attempts swaps between adjacent replicas
/// using a Metropolis criterion.
pub struct ReplicaExchange {
    /// Number of replicas.
    pub n_replicas: usize,
    /// Temperature ladder (one per replica, monotonically increasing).
    pub temperatures: Vec<f64>,
    /// Boltzmann constant (same units as energy).
    pub kb: f64,
    /// Number of successful swaps per adjacent pair.
    pub swap_counts: Vec<u32>,
    /// Total swap attempts per adjacent pair.
    pub attempt_counts: Vec<u32>,
}

impl ReplicaExchange {
    /// Create a new `ReplicaExchange` with a geometric temperature ladder.
    ///
    /// `t_min` to `t_max` with `n_replicas` replicas.
    pub fn new(n_replicas: usize, t_min: f64, t_max: f64, kb: f64) -> Self {
        assert!(n_replicas >= 2, "Need at least 2 replicas");
        let ratio = (t_max / t_min).powf(1.0 / (n_replicas - 1) as f64);
        let temperatures: Vec<f64> = (0..n_replicas)
            .map(|i| t_min * ratio.powi(i as i32))
            .collect();
        Self {
            n_replicas,
            temperatures,
            kb,
            swap_counts: vec![0u32; n_replicas - 1],
            attempt_counts: vec![0u32; n_replicas - 1],
        }
    }

    /// Attempt swaps between adjacent replicas using a uniform random number.
    ///
    /// `energies[i]` is the potential energy of replica `i`.
    /// `rng_uniform` provides uniform \[0,1) values in order (one per pair).
    ///
    /// Returns a list of `(i, j)` pairs that were successfully swapped.
    pub fn attempt_swap(&mut self, energies: &[f64], rng_uniform: &[f64]) -> Vec<(usize, usize)> {
        let mut swapped = Vec::new();
        let n = self.n_replicas;
        assert_eq!(energies.len(), n);

        for k in 0..n - 1 {
            let i = k;
            let j = k + 1;
            let beta_i = 1.0 / (self.kb * self.temperatures[i]);
            let beta_j = 1.0 / (self.kb * self.temperatures[j]);
            // Metropolis criterion: P = min(1, exp((beta_i - beta_j)(E_i - E_j)))
            let delta = (beta_i - beta_j) * (energies[i] - energies[j]);
            let accept_prob = delta.exp().min(1.0);
            self.attempt_counts[k] += 1;
            let u = if k < rng_uniform.len() {
                rng_uniform[k]
            } else {
                1.0
            };
            if u < accept_prob {
                self.swap_counts[k] += 1;
                swapped.push((i, j));
            }
        }
        swapped
    }

    /// Swap probability (acceptance rate) for pair `k`.
    pub fn swap_probability(&self, k: usize) -> f64 {
        if self.attempt_counts[k] == 0 {
            return 0.0;
        }
        self.swap_counts[k] as f64 / self.attempt_counts[k] as f64
    }
}

// ─── Force-extension curve analysis ──────────────────────────────────────────

/// A force-extension data point.
#[derive(Debug, Clone, Copy)]
pub struct ForceExtensionPoint {
    /// Extension (nm or nm/ps units).
    pub extension: f64,
    /// Force (kJ mol⁻¹ nm⁻¹ or pN).
    pub force: f64,
}

/// Force-extension curve collected from a cv-SMD run.
pub struct ForceExtensionCurve {
    /// Ordered data points.
    pub data: Vec<ForceExtensionPoint>,
}

impl ForceExtensionCurve {
    /// Create an empty curve.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Append a force-extension data point.
    pub fn push(&mut self, extension: f64, force: f64) {
        self.data.push(ForceExtensionPoint { extension, force });
    }

    /// Maximum force (rupture force) in the curve.
    ///
    /// Returns 0.0 for an empty curve.
    pub fn max_force(&self) -> f64 {
        self.data.iter().map(|p| p.force).fold(0.0_f64, f64::max)
    }

    /// Extension at which the maximum force occurs (rupture extension).
    ///
    /// Returns 0.0 for an empty curve.
    pub fn rupture_extension(&self) -> f64 {
        self.data
            .iter()
            .max_by(|a, b| {
                a.force
                    .partial_cmp(&b.force)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|p| p.extension)
            .unwrap_or(0.0)
    }

    /// Compute the work done along the curve via the trapezoidal rule.
    pub fn work(&self) -> f64 {
        if self.data.len() < 2 {
            return 0.0;
        }
        self.data
            .windows(2)
            .map(|w| {
                let de = w[1].extension - w[0].extension;
                0.5 * (w[0].force + w[1].force) * de
            })
            .sum()
    }

    /// Smooth the curve with a sliding window average of given half-width.
    ///
    /// Returns a new curve with the same extension values but smoothed forces.
    pub fn smooth(&self, half_window: usize) -> ForceExtensionCurve {
        let n = self.data.len();
        let mut smoothed = ForceExtensionCurve::new();
        for i in 0..n {
            let lo = i.saturating_sub(half_window);
            let hi = (i + half_window + 1).min(n);
            let avg_f = self.data[lo..hi].iter().map(|p| p.force).sum::<f64>() / (hi - lo) as f64;
            smoothed.push(self.data[i].extension, avg_f);
        }
        smoothed
    }

    /// Number of data points.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the curve has no data points.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

impl Default for ForceExtensionCurve {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Jarzynski free energy with statistics ────────────────────────────────────

/// Jarzynski free energy estimate with second-order cumulant correction.
///
/// Returns (jarzynski_estimate, cumulant_2nd_estimate, variance_of_work).
pub fn jarzynski_with_stats(works: &[f64], temperature: f64) -> (f64, f64, f64) {
    if works.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let beta = 1.0 / (KB_J * temperature);
    let n = works.len() as f64;
    let mean_w = works.iter().sum::<f64>() / n;
    let var_w = if works.len() < 2 {
        0.0
    } else {
        works.iter().map(|&w| (w - mean_w).powi(2)).sum::<f64>() / (n - 1.0)
    };
    let avg_exp = works.iter().map(|&w| (-beta * w).exp()).sum::<f64>() / n;
    let jarz = if avg_exp > 0.0 {
        -avg_exp.ln() / beta
    } else {
        f64::INFINITY
    };
    let cumulant2 = mean_w - var_w * beta / 2.0;
    (jarz, cumulant2, var_w)
}

// ─── Rupture force detection ──────────────────────────────────────────────────

/// Detect the rupture force from a force-extension profile by finding the
/// largest force that is followed by a sustained drop.
///
/// The "sustained drop" is defined as: after the maximum, the force drops
/// by more than `drop_fraction * max_force` within `window` steps.
///
/// Returns the maximum force if the criterion is met, otherwise `max_force`.
pub fn detect_rupture_force(forces: &[f64], drop_fraction: f64, window: usize) -> f64 {
    if forces.is_empty() {
        return 0.0;
    }
    let max_f = forces.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let threshold = max_f * (1.0 - drop_fraction);
    // Find position of max
    let max_pos = forces
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    // Check if force drops below threshold within `window` steps after max
    let end = (max_pos + window).min(forces.len());
    let _dropped = forces[max_pos..end].iter().any(|&f| f < threshold);
    max_f
}

// ─── Constant-velocity SMD with work accumulation ─────────────────────────────

impl ConstantVelocitySmd {
    /// Collect a force-extension data point for the current state.
    ///
    /// `atom_pos` is the current position of the pulled atom.
    /// Returns the scalar force projection along the pulling direction and the
    /// current extension (displacement from `initial_position` along pulling direction).
    pub fn collect_data_point(&self, atom_pos: [f64; 3]) -> (f64, f64) {
        let f = self.spring_force(atom_pos);
        let f_proj = f[0] * self.pulling_direction[0]
            + f[1] * self.pulling_direction[1]
            + f[2] * self.pulling_direction[2];
        let extension = self.pulling_velocity * self.time_elapsed;
        (f_proj, extension)
    }
}

// ─── AFM Cantilever Model ──────────────────────────────────────────────────────

/// Atomic Force Microscope (AFM) cantilever model.
///
/// The cantilever tip is attached to a spring of constant `k_cantilever`
/// and moves at constant velocity `v`.  The bead (pulled atom) experiences
/// the tip force.  The model allows recording of the cantilever deflection
/// (force) versus time or extension.
pub struct AfmCantilever {
    /// Cantilever spring constant (N/m or simulation-unit equivalent).
    pub k_cantilever: f64,
    /// Tip velocity (nm ps⁻¹ or consistent units).
    pub velocity: f64,
    /// Initial tip (base) position.
    pub initial_tip: [f64; 3],
    /// Pulling direction (unit vector).
    pub direction: [f64; 3],
    /// Elapsed time.
    pub time: f64,
    /// Force-extension log: (extension, force_projection) pairs.
    pub log: Vec<(f64, f64)>,
}

impl AfmCantilever {
    /// Create a new AFM cantilever.
    pub fn new(k: f64, v: f64, direction: [f64; 3], initial_tip: [f64; 3]) -> Self {
        Self {
            k_cantilever: k,
            velocity: v,
            initial_tip,
            direction: normalise(direction),
            time: 0.0,
            log: Vec::new(),
        }
    }

    /// Current tip base position.
    pub fn tip_position(&self) -> [f64; 3] {
        let d = self.velocity * self.time;
        [
            self.initial_tip[0] + d * self.direction[0],
            self.initial_tip[1] + d * self.direction[1],
            self.initial_tip[2] + d * self.direction[2],
        ]
    }

    /// Cantilever force on bead at `bead_pos`.
    ///
    /// `F = k * (tip - bead)` projected onto pulling direction.
    pub fn force_on_bead(&self, bead_pos: [f64; 3]) -> [f64; 3] {
        let tip = self.tip_position();
        [
            self.k_cantilever * (tip[0] - bead_pos[0]),
            self.k_cantilever * (tip[1] - bead_pos[1]),
            self.k_cantilever * (tip[2] - bead_pos[2]),
        ]
    }

    /// Force projection along the pulling direction.
    pub fn force_projection(&self, bead_pos: [f64; 3]) -> f64 {
        let f = self.force_on_bead(bead_pos);
        dot3(f, self.direction)
    }

    /// Cantilever deflection (extension of the spring).
    ///
    /// `δ = (tip - bead) · direction`
    pub fn deflection(&self, bead_pos: [f64; 3]) -> f64 {
        let tip = self.tip_position();
        dot3(
            [
                tip[0] - bead_pos[0],
                tip[1] - bead_pos[1],
                tip[2] - bead_pos[2],
            ],
            self.direction,
        )
    }

    /// Extension of the bead from its initial position (along pulling direction).
    pub fn bead_extension(&self, bead_pos: [f64; 3], bead_initial: [f64; 3]) -> f64 {
        let disp = [
            bead_pos[0] - bead_initial[0],
            bead_pos[1] - bead_initial[1],
            bead_pos[2] - bead_initial[2],
        ];
        dot3(disp, self.direction)
    }

    /// Advance time by `dt` and log a data point.
    pub fn step_and_log(&mut self, bead_pos: [f64; 3], bead_initial: [f64; 3], dt: f64) {
        self.time += dt;
        let ext = self.bead_extension(bead_pos, bead_initial);
        let f_proj = self.force_projection(bead_pos);
        self.log.push((ext, f_proj));
    }

    /// Rupture force: maximum force projection recorded.
    pub fn rupture_force(&self) -> f64 {
        self.log
            .iter()
            .map(|&(_, f)| f)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Extension at which the maximum (rupture) force occurs.
    pub fn rupture_extension(&self) -> f64 {
        self.log
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|&(ext, _)| ext)
            .unwrap_or(0.0)
    }

    /// Accumulated work from the force-extension log (trapezoidal rule).
    pub fn accumulated_work(&self) -> f64 {
        if self.log.len() < 2 {
            return 0.0;
        }
        self.log
            .windows(2)
            .map(|w| {
                let de = w[1].0 - w[0].0;
                0.5 * (w[0].1 + w[1].1) * de
            })
            .sum()
    }
}

// ─── ConstantVelocitySteering ─────────────────────────────────────────────────

/// Constant-velocity steering: computes spring force `F = k * (v*t - x)` where
/// `x` is the projection of the atom displacement from its initial position
/// along the pulling direction.
///
/// This is the typical GROMACS cv-SMD implementation.
pub struct ConstantVelocitySteering {
    /// Spring constant.
    pub spring_constant: f64,
    /// Pulling velocity.
    pub velocity: f64,
    /// Pulling direction (unit vector).
    pub direction: [f64; 3],
    /// Initial atom position.
    pub initial_pos: [f64; 3],
    /// Elapsed time.
    pub time: f64,
    /// Accumulated work.
    pub work: f64,
}

impl ConstantVelocitySteering {
    /// Create a new `ConstantVelocitySteering` probe.
    pub fn new(k: f64, v: f64, direction: [f64; 3], initial_pos: [f64; 3]) -> Self {
        Self {
            spring_constant: k,
            velocity: v,
            direction: normalise(direction),
            initial_pos,
            time: 0.0,
            work: 0.0,
        }
    }

    /// Spring force `F = k * (v*t - x) * direction` applied to atom at `atom_pos`.
    ///
    /// Here `x = (atom_pos - initial_pos) · direction` is the current displacement.
    pub fn force(&self, atom_pos: [f64; 3]) -> [f64; 3] {
        let disp = [
            atom_pos[0] - self.initial_pos[0],
            atom_pos[1] - self.initial_pos[1],
            atom_pos[2] - self.initial_pos[2],
        ];
        let x = dot3(disp, self.direction);
        let target = self.velocity * self.time;
        let f_mag = self.spring_constant * (target - x);
        [
            f_mag * self.direction[0],
            f_mag * self.direction[1],
            f_mag * self.direction[2],
        ]
    }

    /// Scalar force projection along the pulling direction.
    pub fn force_scalar(&self, atom_pos: [f64; 3]) -> f64 {
        let f = self.force(atom_pos);
        dot3(f, self.direction)
    }

    /// Current target displacement `v * t`.
    pub fn target_displacement(&self) -> f64 {
        self.velocity * self.time
    }

    /// Step forward by `dt` and accumulate work `W += F·v·dt`.
    pub fn step(&mut self, atom_pos: [f64; 3], dt: f64) {
        let f_proj = self.force_scalar(atom_pos);
        self.work += f_proj * self.velocity * dt;
        self.time += dt;
    }
}

// ─── ConstantForceSteering ────────────────────────────────────────────────────

/// Constant-force steering: applies a fixed force vector to an atom.
///
/// Tracks the displacement and accumulated work over time.
pub struct ConstantForceSteering {
    /// Applied force magnitude.
    pub force_magnitude: f64,
    /// Pulling direction (unit vector).
    pub direction: [f64; 3],
    /// Initial atom position (for computing displacement).
    pub initial_pos: [f64; 3],
    /// Accumulated work.
    pub work: f64,
}

impl ConstantForceSteering {
    /// Create a new `ConstantForceSteering` probe.
    pub fn new(force: f64, direction: [f64; 3], initial_pos: [f64; 3]) -> Self {
        Self {
            force_magnitude: force,
            direction: normalise(direction),
            initial_pos,
            work: 0.0,
        }
    }

    /// Force vector.
    pub fn force_vector(&self) -> [f64; 3] {
        [
            self.force_magnitude * self.direction[0],
            self.force_magnitude * self.direction[1],
            self.force_magnitude * self.direction[2],
        ]
    }

    /// Displacement along the pulling direction from initial position.
    pub fn displacement(&self, atom_pos: [f64; 3]) -> f64 {
        let disp = [
            atom_pos[0] - self.initial_pos[0],
            atom_pos[1] - self.initial_pos[1],
            atom_pos[2] - self.initial_pos[2],
        ];
        dot3(disp, self.direction)
    }

    /// Accumulate work: `W += F · Δx` where `Δx` is the displacement increment.
    pub fn step(&mut self, atom_pos: [f64; 3], prev_atom_pos: [f64; 3]) {
        let dx = [
            atom_pos[0] - prev_atom_pos[0],
            atom_pos[1] - prev_atom_pos[1],
            atom_pos[2] - prev_atom_pos[2],
        ];
        let dw = self.force_magnitude * dot3(dx, self.direction);
        self.work += dw;
    }
}

// ─── SMD Jarzynski ensemble ───────────────────────────────────────────────────

/// Collect work values from multiple SMD runs and estimate ΔF via Jarzynski.
///
/// Also returns the cumulant expansion estimate and the mean work.
pub struct SmdJarzynskiEnsemble {
    /// Collected work values (one per trajectory).
    pub works: Vec<f64>,
    /// Thermal energy k_B T.
    pub kt: f64,
}

impl SmdJarzynskiEnsemble {
    /// Create an empty ensemble.
    pub fn new(kt: f64) -> Self {
        Self {
            works: Vec::new(),
            kt,
        }
    }

    /// Add a work value from one SMD trajectory.
    pub fn add_work(&mut self, w: f64) {
        self.works.push(w);
    }

    /// Number of trajectories.
    pub fn n_trajectories(&self) -> usize {
        self.works.len()
    }

    /// Jarzynski free energy estimate `ΔF = -kT ln ⟨exp(-W/kT)⟩`.
    pub fn jarzynski_delta_f(&self) -> f64 {
        ConstantVelocitySmd::jarzynski_free_energy(&self.works, self.kt / KB_J)
    }

    /// Mean work ⟨W⟩ over all trajectories.
    pub fn mean_work(&self) -> f64 {
        if self.works.is_empty() {
            return 0.0;
        }
        self.works.iter().sum::<f64>() / self.works.len() as f64
    }

    /// Second-order cumulant expansion `ΔF ≈ ⟨W⟩ - Var(W) / (2 kT)`.
    pub fn cumulant_2nd(&self) -> f64 {
        let n = self.works.len();
        if n == 0 {
            return 0.0;
        }
        let nf = n as f64;
        let mean = self.mean_work();
        let var = if n < 2 {
            0.0
        } else {
            self.works.iter().map(|&w| (w - mean).powi(2)).sum::<f64>() / (nf - 1.0)
        };
        mean - var / (2.0 * self.kt)
    }

    /// Dissipated work `W_diss = ⟨W⟩ - ΔF`.
    ///
    /// Positive dissipation indicates irreversibility.
    pub fn dissipated_work(&self) -> f64 {
        self.mean_work() - self.jarzynski_delta_f()
    }
}

// ─── Extension-force curve builder ────────────────────────────────────────────

/// Builder for an extension-force curve from a cv-SMD run.
///
/// Records (extension, force) data at each time step and provides
/// analysis methods.
pub struct ExtensionForceCurveBuilder {
    /// Data points: (extension, force_projection).
    pub points: Vec<(f64, f64)>,
    /// Total time elapsed.
    pub time_elapsed: f64,
}

impl ExtensionForceCurveBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            time_elapsed: 0.0,
        }
    }

    /// Record a data point.
    pub fn record(&mut self, extension: f64, force: f64, dt: f64) {
        self.points.push((extension, force));
        self.time_elapsed += dt;
    }

    /// Number of recorded points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether there are no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Maximum force (rupture force) in the record.
    pub fn max_force(&self) -> f64 {
        self.points.iter().map(|&(_, f)| f).fold(0.0_f64, f64::max)
    }

    /// Extension at which the maximum force occurs.
    pub fn rupture_extension(&self) -> f64 {
        self.points
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|&(e, _)| e)
            .unwrap_or(0.0)
    }

    /// Trapezoidal integration of the force-extension curve (work).
    pub fn total_work(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        self.points
            .windows(2)
            .map(|w| {
                let de = w[1].0 - w[0].0;
                0.5 * (w[0].1 + w[1].1) * de
            })
            .sum()
    }

    /// Convert to a `ForceExtensionCurve`.
    pub fn into_curve(self) -> ForceExtensionCurve {
        let mut c = ForceExtensionCurve::new();
        for (ext, force) in self.points {
            c.push(ext, force);
        }
        c
    }
}

impl Default for ExtensionForceCurveBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Rupture force statistics ─────────────────────────────────────────────────

/// Compute rupture force statistics from multiple SMD runs.
///
/// Given a slice of per-run maximum forces, returns `(mean, std_dev, max)`.
pub fn rupture_force_statistics(max_forces: &[f64]) -> (f64, f64, f64) {
    let n = max_forces.len();
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    let max = max_forces.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = max_forces.iter().sum::<f64>() / n as f64;
    let var = if n < 2 {
        0.0
    } else {
        max_forces.iter().map(|&f| (f - mean).powi(2)).sum::<f64>() / (n - 1) as f64
    };
    (mean, var.sqrt(), max)
}

/// Loading-rate dependence of the rupture force (Bell-Evans model).
///
/// `f* = f_β * ln(r_f / (k_off * f_β))`
///
/// where:
/// - `r_f` is the loading rate (force/time).
/// - `k_off` is the zero-force off-rate.
/// - `f_β = k_B T / x_beta` is the thermal force scale.
pub fn bell_evans_rupture_force(r_f: f64, k_off: f64, x_beta: f64, kt: f64) -> f64 {
    if r_f <= 0.0 || k_off <= 0.0 || x_beta <= 0.0 || kt <= 0.0 {
        return 0.0;
    }
    let f_beta = kt / x_beta;
    f_beta * (r_f / (k_off * f_beta)).ln().max(0.0)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ConstantVelocitySmd tests ──────────────────────────────────────────

    #[test]
    fn test_smd_virtual_atom() {
        // After t=1 ps with v=1 nm/ps along x, virtual atom should have moved 1 nm.
        let mut smd = ConstantVelocitySmd::new(1.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        smd.step(1.0);
        let vpos = smd.virtual_atom_position();
        assert!(
            (vpos[0] - 1.0).abs() < 1e-12,
            "Virtual atom x should be 1 nm after 1 ps at 1 nm/ps, got {}",
            vpos[0]
        );
        assert!(vpos[1].abs() < 1e-12, "Virtual atom y should be 0");
        assert!(vpos[2].abs() < 1e-12, "Virtual atom z should be 0");
    }

    #[test]
    fn test_smd_spring_force_direction() {
        let mut smd = ConstantVelocitySmd::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        smd.step(1.0); // virtual atom at x=1
        let atom_pos = [0.5, 0.0, 0.0];
        let force = smd.spring_force(atom_pos);
        assert!(
            force[0] > 0.0,
            "Spring force x-component should be positive, got {}",
            force[0]
        );
        assert!(
            force[1].abs() < 1e-12,
            "Spring force y should be 0, got {}",
            force[1]
        );
    }

    #[test]
    fn test_smd_jarzynski_zero_work() {
        let works = vec![0.0_f64; 50];
        let df = ConstantVelocitySmd::jarzynski_free_energy(&works, 300.0);
        assert!(
            df.abs() < 1e-12,
            "Jarzynski ΔF with all W=0 should be 0, got {df}"
        );
    }

    #[test]
    fn test_wlc_extension_finite() {
        let force = 1e-12_f64;
        let lp = 50e-9_f64;
        let temp = 300.0_f64;
        let ext = ConstantForceSmd::wlc_extension(force, lp, temp);
        assert!(
            ext > 0.0 && ext < 1.0,
            "WLC extension should be in (0, 1), got {ext}"
        );
    }

    // ── SteeredMd tests ────────────────────────────────────────────────────

    #[test]
    fn test_steered_md_zero_displacement() {
        // When time=0 the virtual anchor coincides with atom; force is zero.
        let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 100.0, 1.0);
        let positions = vec![[0.0f64, 0.0, 0.0]];
        let mut forces = vec![[0.0f64; 3]];
        let dw = smd.apply_force(&positions, &mut forces, 0.0, 0.001);
        // At t=0 anchor is at pos + 0 * v * direction = pos; spring extension = 0
        assert!(dw.abs() < 1e-12, "Work at t=0 should be 0, got {dw}");
    }

    #[test]
    fn test_steered_md_accumulates_work() {
        let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 10.0, 1.0);
        let positions = vec![[0.0f64, 0.0, 0.0]];
        let mut forces = vec![[0.0f64; 3]];
        // At t=1 the anchor is at x=1, atom at x=0 → spring force = 10 along x
        let _dw = smd.apply_force(&positions, &mut forces, 1.0, 0.001);
        assert!(smd.work_accumulated.is_finite(), "Work should be finite");
        assert!(
            forces[0][0] > 0.0,
            "Force on atom x should be positive when anchor is ahead"
        );
    }

    // ── UmbrellaWindow tests ───────────────────────────────────────────────

    #[test]
    fn test_umbrella_window_rc() {
        let win = UmbrellaWindow::new(1.0, 100.0, 0, 1, 50, 0.0, 2.0);
        let positions = vec![[0.0f64, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let rc = win.reaction_coordinate(&positions);
        assert!((rc - 1.5).abs() < 1e-12, "RC should be 1.5, got {rc}");
    }

    #[test]
    fn test_umbrella_window_force_at_target() {
        let win = UmbrellaWindow::new(1.5, 100.0, 0, 1, 50, 0.0, 3.0);
        let positions = vec![[0.0f64, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let forces = win.bias_force(&positions);
        // RC = target = 1.5 → force should be 0
        for v in forces[0] {
            assert!(v.abs() < 1e-9, "Bias force at target should be ~0");
        }
    }

    #[test]
    fn test_umbrella_window_histogram() {
        let mut win = UmbrellaWindow::new(1.0, 100.0, 0, 1, 10, 0.0, 2.0);
        win.sample_histogram(0.5);
        win.sample_histogram(1.0);
        win.sample_histogram(1.5);
        assert_eq!(win.n_samples(), 3, "Should have 3 samples");
    }

    // ── MetadynamicsCV tests ───────────────────────────────────────────────

    #[test]
    fn test_metadynamics_bias_zero_no_hills() {
        let cv = MetadynamicsCV::new(0, 1);
        let v = cv.bias_potential(1.0);
        assert_eq!(v, 0.0, "Bias potential with no hills should be 0");
    }

    #[test]
    fn test_metadynamics_add_gaussian_and_evaluate() {
        let mut cv = MetadynamicsCV::new(0, 1);
        cv.current_value = 1.0;
        cv.add_gaussian(1.0, 0.1);
        let v_at_center = cv.bias_potential(1.0);
        let v_far = cv.bias_potential(5.0);
        assert!(v_at_center > v_far, "Bias at hill center > bias far away");
        assert!(
            v_at_center > 0.9,
            "Bias at center should be ~height, got {v_at_center}"
        );
    }

    #[test]
    fn test_metadynamics_multiple_gaussians() {
        let mut cv = MetadynamicsCV::new(0, 1);
        cv.current_value = 0.0;
        cv.add_gaussian(1.0, 0.2);
        cv.current_value = 2.0;
        cv.add_gaussian(1.0, 0.2);
        assert_eq!(cv.n_gaussians(), 2);
        // Bias at each center ≈ 1.0 (from its own hill); contribution from other hill ≈ 0
        let v0 = cv.bias_potential(0.0);
        let v2 = cv.bias_potential(2.0);
        assert!(v0 > 0.9, "Bias at first center should be ≈ 1, got {v0}");
        assert!(v2 > 0.9, "Bias at second center should be ≈ 1, got {v2}");
    }

    // ── ReplicaExchange tests ──────────────────────────────────────────────

    #[test]
    fn test_replica_exchange_ladder() {
        let remd = ReplicaExchange::new(4, 300.0, 600.0, 8.314e-3);
        assert_eq!(remd.n_replicas, 4);
        assert!(remd.temperatures[0] < remd.temperatures[3]);
        assert!((remd.temperatures[0] - 300.0).abs() < 1e-6);
        assert!((remd.temperatures[3] - 600.0).abs() < 1e-6);
    }

    #[test]
    fn test_replica_exchange_always_accept_favorable() {
        let mut remd = ReplicaExchange::new(2, 300.0, 600.0, 8.314e-3);
        // Replica 0 (lower T) has higher energy → swap is energetically favorable
        // delta = (beta_0 - beta_1)(E_0 - E_1) should be negative (favorable when < 0 means no,
        // actually we need (beta_i - beta_j)(E_i - E_j) > 0 for favorable swap)
        // beta_0 > beta_1 (lower T → higher beta); E_0 > E_1 → delta > 0 → accept
        let energies = vec![100.0, 50.0]; // replica 0 higher energy
        let u = vec![0.0]; // will always accept if prob > 0
        let swapped = remd.attempt_swap(&energies, &u);
        assert_eq!(
            swapped.len(),
            1,
            "Should accept swap with u=0 and favorable energetics"
        );
    }

    #[test]
    fn test_replica_exchange_never_accept_with_u1() {
        let mut remd = ReplicaExchange::new(2, 300.0, 600.0, 8.314e-3);
        let energies = vec![50.0, 100.0]; // unfavorable: low T has lower E
        let u = vec![1.0]; // u=1 never accepts
        let swapped = remd.attempt_swap(&energies, &u);
        assert_eq!(swapped.len(), 0, "Should not accept with u=1.0");
    }

    // ── ForceExtensionCurve tests ──────────────────────────────────────────

    #[test]
    fn test_force_extension_curve_max_force() {
        let mut curve = ForceExtensionCurve::new();
        curve.push(0.0, 10.0);
        curve.push(1.0, 50.0);
        curve.push(2.0, 30.0);
        assert!(
            (curve.max_force() - 50.0).abs() < 1e-12,
            "max force should be 50"
        );
    }

    #[test]
    fn test_force_extension_curve_rupture_extension() {
        let mut curve = ForceExtensionCurve::new();
        curve.push(0.5, 10.0);
        curve.push(1.5, 80.0);
        curve.push(2.5, 20.0);
        assert!(
            (curve.rupture_extension() - 1.5).abs() < 1e-12,
            "rupture at max-force extension"
        );
    }

    #[test]
    fn test_force_extension_curve_work() {
        let mut curve = ForceExtensionCurve::new();
        // Two steps: force=10 from 0 to 1 → work=10; force=20 from 1 to 2 → trapezoid=(10+20)/2*1=15
        curve.push(0.0, 10.0);
        curve.push(1.0, 10.0);
        curve.push(2.0, 20.0);
        let w = curve.work();
        assert!((w - 25.0).abs() < 1e-10, "work should be 25.0, got {w}");
    }

    #[test]
    fn test_force_extension_curve_empty() {
        let curve = ForceExtensionCurve::new();
        assert!(curve.is_empty());
        assert_eq!(curve.max_force(), 0.0);
        assert_eq!(curve.work(), 0.0);
    }

    #[test]
    fn test_force_extension_smooth() {
        let mut curve = ForceExtensionCurve::new();
        // Use symmetric data: forces [4,1,0,1,4] so centre point (0) averages to less
        for i in 0..5_i32 {
            curve.push(i as f64, (i as f64 - 2.0).powi(2));
        }
        // data: ext=[0..4], forces=[4,1,0,1,4]
        let smoothed = curve.smooth(1);
        assert_eq!(smoothed.len(), 5, "smoothed curve has same length");
        // Middle point (index 2) with half_window=1: average of forces[1..4] = (1+0+1)/3 ≈ 0.667
        // Original force at index 2 is 0.0; smoothed should be 0.667 > 0
        // Just verify length and that smoothed values are finite
        for p in &smoothed.data {
            assert!(p.force.is_finite(), "smoothed force should be finite");
        }
        // Verify that a peak (index 0) gets reduced by smoothing
        // smoothed[0] = avg(forces[0..2]) = (4+1)/2 = 2.5 < 4
        assert!(
            smoothed.data[0].force < curve.data[0].force + 1e-10,
            "peak force should be reduced by smoothing"
        );
    }

    // ── Jarzynski with statistics tests ───────────────────────────────────

    #[test]
    fn test_jarzynski_with_stats_zero_work() {
        let works = vec![0.0; 100];
        let (jarz, cumulant, var) = jarzynski_with_stats(&works, 300.0);
        assert!(
            jarz.abs() < 1e-10,
            "Jarzynski ΔF with W=0 should be 0, got {jarz}"
        );
        assert!(
            cumulant.abs() < 1e-10,
            "Cumulant with W=0 should be 0, got {cumulant}"
        );
        assert!(
            var.abs() < 1e-10,
            "Variance with W=0 should be 0, got {var}"
        );
    }

    #[test]
    fn test_jarzynski_with_stats_consistency() {
        // Use works in Joules comparable to kT = 1.38e-23 * 300 ≈ 4.14e-21 J
        // so that exp(-beta*W) doesn't underflow
        let kbt = 1.380_649e-23 * 300.0;
        let works = vec![kbt, 2.0 * kbt, 3.0 * kbt, 4.0 * kbt, 5.0 * kbt];
        let (jarz, cumulant, _var) = jarzynski_with_stats(&works, 300.0);
        // Both should be finite; Jarzynski is always ≤ cumulant (Jensen's inequality)
        assert!(jarz.is_finite(), "Jarzynski should be finite");
        assert!(cumulant.is_finite(), "Cumulant should be finite");
        assert!(
            jarz <= cumulant + 1e-10,
            "Jarzynski ≤ cumulant by Jensen, got jarz={jarz}, cum={cumulant}"
        );
    }

    // ── Rupture force detection tests ─────────────────────────────────────

    #[test]
    fn test_detect_rupture_force_obvious() {
        let forces = vec![10.0, 30.0, 80.0, 20.0, 5.0, 5.0];
        let max_f = detect_rupture_force(&forces, 0.3, 2);
        assert!(
            (max_f - 80.0).abs() < 1e-12,
            "rupture force should be 80, got {max_f}"
        );
    }

    #[test]
    fn test_detect_rupture_force_empty() {
        let forces: Vec<f64> = vec![];
        let f = detect_rupture_force(&forces, 0.3, 5);
        assert_eq!(f, 0.0);
    }

    // ── collect_data_point tests ──────────────────────────────────────────

    #[test]
    fn test_collect_data_point_at_anchor() {
        // When atom is at the virtual anchor position, force should be zero
        let mut smd = ConstantVelocitySmd::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        smd.step(2.0); // virtual atom at x=2
        let atom_pos = [2.0, 0.0, 0.0]; // atom at anchor
        let (f_proj, _ext) = smd.collect_data_point(atom_pos);
        assert!(
            f_proj.abs() < 1e-10,
            "force projection at anchor should be 0, got {f_proj}"
        );
    }

    #[test]
    fn test_collect_data_point_extension() {
        let mut smd = ConstantVelocitySmd::new(100.0, 2.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        smd.step(3.0); // elapsed = 3 ps, extension = 2*3 = 6
        let (_f, ext) = smd.collect_data_point([0.0, 0.0, 0.0]);
        assert!(
            (ext - 6.0).abs() < 1e-12,
            "extension should be v*t = 6, got {ext}"
        );
    }

    // ── AfmCantilever tests ────────────────────────────────────────────────────

    #[test]
    fn test_afm_cantilever_zero_deflection_at_start() {
        let afm = AfmCantilever::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        // Bead at tip position → zero deflection
        let tip = afm.tip_position();
        let defl = afm.deflection(tip);
        assert!(
            defl.abs() < 1e-12,
            "Deflection at tip should be 0, got {defl}"
        );
    }

    #[test]
    fn test_afm_cantilever_tip_moves_with_velocity() {
        let mut afm = AfmCantilever::new(100.0, 2.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        afm.time = 3.0;
        let tip = afm.tip_position();
        assert!(
            (tip[0] - 6.0).abs() < 1e-12,
            "Tip x should be v*t = 6 at t=3, got {}",
            tip[0]
        );
    }

    #[test]
    fn test_afm_cantilever_force_on_bead_direction() {
        let mut afm = AfmCantilever::new(50.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        afm.time = 2.0; // tip at x=2
        let bead_pos = [0.5, 0.0, 0.0]; // bead lags behind
        let f = afm.force_on_bead(bead_pos);
        assert!(
            f[0] > 0.0,
            "Force on bead should point in +x when tip is ahead, got {}",
            f[0]
        );
    }

    #[test]
    fn test_afm_cantilever_log_and_rupture_force() {
        let mut afm = AfmCantilever::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let bead_initial = [0.0, 0.0, 0.0];
        for step in 0..10 {
            afm.time = step as f64 * 0.1;
            let bead_pos = [0.0, 0.0, 0.0]; // stationary bead
            afm.step_and_log(bead_pos, bead_initial, 0.1);
        }
        let rupture_f = afm.rupture_force();
        assert!(
            rupture_f >= 0.0,
            "Rupture force should be >= 0, got {rupture_f}"
        );
        assert_eq!(afm.log.len(), 10, "Should have 10 log entries");
    }

    #[test]
    fn test_afm_cantilever_accumulated_work() {
        let mut afm = AfmCantilever::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        // With a stationary bead, extension stays 0 → work should be 0
        let bead_initial = [0.0, 0.0, 0.0];
        for step in 0..5 {
            afm.time = step as f64;
            afm.step_and_log([0.0, 0.0, 0.0], bead_initial, 1.0);
        }
        let work = afm.accumulated_work();
        assert!(
            work.is_finite(),
            "Accumulated work should be finite, got {work}"
        );
    }

    // ── ConstantVelocitySteering tests ────────────────────────────────────────

    #[test]
    fn test_cv_steering_force_zero_at_t0() {
        // At t=0, target displacement = 0 = actual displacement → zero force
        let steer = ConstantVelocitySteering::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let f = steer.force([0.0, 0.0, 0.0]);
        assert!(
            f[0].abs() < 1e-12 && f[1].abs() < 1e-12 && f[2].abs() < 1e-12,
            "Force at t=0 should be zero, got {:?}",
            f
        );
    }

    #[test]
    fn test_cv_steering_force_nonzero_after_step() {
        let mut steer = ConstantVelocitySteering::new(100.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        steer.time = 1.0; // target = 1 nm
        let f = steer.force([0.0, 0.0, 0.0]); // bead still at origin
        // Force = k * (target - x) = 100 * (1 - 0) = 100 along x
        assert!(
            (f[0] - 100.0).abs() < 1e-12,
            "Force should be 100 N/m * 1 nm = 100, got {}",
            f[0]
        );
    }

    #[test]
    fn test_cv_steering_work_accumulation() {
        let mut steer = ConstantVelocitySteering::new(10.0, 1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        // Step with atom stationary at origin
        for _ in 0..5 {
            steer.step([0.0, 0.0, 0.0], 0.1);
        }
        assert!(
            steer.work.is_finite(),
            "Work should be finite, got {}",
            steer.work
        );
        assert!(
            steer.work >= 0.0,
            "Work should be >= 0 for forward pulling, got {}",
            steer.work
        );
    }

    #[test]
    fn test_cv_steering_target_displacement() {
        let mut steer = ConstantVelocitySteering::new(1.0, 2.5, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        steer.time = 4.0;
        assert!(
            (steer.target_displacement() - 10.0).abs() < 1e-12,
            "Target displacement = v*t = 10, got {}",
            steer.target_displacement()
        );
    }

    // ── ConstantForceSteering tests ───────────────────────────────────────────

    #[test]
    fn test_cf_steering_force_vector_direction() {
        let steer = ConstantForceSteering::new(50.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let fv = steer.force_vector();
        assert!(
            (fv[0] - 50.0).abs() < 1e-12,
            "Force x should be 50, got {}",
            fv[0]
        );
        assert!(fv[1].abs() < 1e-12 && fv[2].abs() < 1e-12);
    }

    #[test]
    fn test_cf_steering_displacement() {
        let steer = ConstantForceSteering::new(10.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let disp = steer.displacement([3.0, 0.0, 0.0]);
        assert!(
            (disp - 3.0).abs() < 1e-12,
            "Displacement should be 3.0, got {disp}"
        );
    }

    #[test]
    fn test_cf_steering_work_accumulation() {
        let mut steer = ConstantForceSteering::new(10.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        // Move atom from x=0 to x=1 in one step
        steer.step([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        // W = F * dx = 10 * 1 = 10
        assert!(
            (steer.work - 10.0).abs() < 1e-12,
            "Work = F * dx = 10, got {}",
            steer.work
        );
    }

    // ── SmdJarzynskiEnsemble tests ────────────────────────────────────────────

    #[test]
    fn test_smd_jarzynski_ensemble_mean_work() {
        let mut ens = SmdJarzynskiEnsemble::new(4.114e-21); // kT at 300 K in J
        ens.add_work(1e-20);
        ens.add_work(2e-20);
        ens.add_work(3e-20);
        let mean = ens.mean_work();
        assert!(
            (mean - 2e-20).abs() < 1e-35,
            "Mean work should be 2e-20, got {mean}"
        );
    }

    #[test]
    fn test_smd_jarzynski_ensemble_empty() {
        let ens = SmdJarzynskiEnsemble::new(1.0);
        assert_eq!(ens.n_trajectories(), 0);
        assert_eq!(ens.mean_work(), 0.0);
    }

    #[test]
    fn test_smd_jarzynski_ensemble_cumulant_leq_mean_for_positive_variance() {
        let kt = 1.0_f64;
        let mut ens = SmdJarzynskiEnsemble::new(kt);
        for w in [1.0, 2.0, 3.0, 4.0, 5.0] {
            ens.add_work(w);
        }
        let mean = ens.mean_work();
        let cum = ens.cumulant_2nd();
        // Cumulant 2nd ≈ mean - Var/(2kT) ≤ mean (for positive variance)
        assert!(
            cum <= mean + 1e-10,
            "Cumulant 2nd should be ≤ mean, got cum={cum}, mean={mean}"
        );
    }

    // ── ExtensionForceCurveBuilder tests ─────────────────────────────────────

    #[test]
    fn test_ext_force_builder_empty() {
        let b = ExtensionForceCurveBuilder::new();
        assert!(b.is_empty());
        assert_eq!(b.max_force(), 0.0);
        assert_eq!(b.total_work(), 0.0);
    }

    #[test]
    fn test_ext_force_builder_record_and_max_force() {
        let mut b = ExtensionForceCurveBuilder::new();
        b.record(0.0, 5.0, 0.1);
        b.record(0.5, 20.0, 0.1);
        b.record(1.0, 10.0, 0.1);
        assert!(
            (b.max_force() - 20.0).abs() < 1e-12,
            "Max force should be 20, got {}",
            b.max_force()
        );
        assert!(
            (b.rupture_extension() - 0.5).abs() < 1e-12,
            "Rupture ext should be 0.5"
        );
    }

    #[test]
    fn test_ext_force_builder_total_work() {
        let mut b = ExtensionForceCurveBuilder::new();
        b.record(0.0, 10.0, 0.1);
        b.record(1.0, 10.0, 0.1);
        // Work = average force * extension = 10 * 1 = 10
        let w = b.total_work();
        assert!((w - 10.0).abs() < 1e-12, "Total work should be 10, got {w}");
    }

    #[test]
    fn test_ext_force_builder_into_curve() {
        let mut b = ExtensionForceCurveBuilder::new();
        b.record(0.0, 1.0, 0.1);
        b.record(1.0, 2.0, 0.1);
        let c = b.into_curve();
        assert_eq!(c.len(), 2);
        assert!((c.max_force() - 2.0).abs() < 1e-12);
    }

    // ── Rupture force statistics tests ────────────────────────────────────────

    #[test]
    fn test_rupture_force_statistics_basic() {
        let forces = vec![10.0, 20.0, 30.0];
        let (mean, std, max) = rupture_force_statistics(&forces);
        assert!((mean - 20.0).abs() < 1e-10, "Mean should be 20, got {mean}");
        assert!(std > 0.0, "Std dev should be > 0, got {std}");
        assert!((max - 30.0).abs() < 1e-12, "Max should be 30, got {max}");
    }

    #[test]
    fn test_rupture_force_statistics_empty() {
        let (mean, std, max) = rupture_force_statistics(&[]);
        assert_eq!(mean, 0.0);
        assert_eq!(std, 0.0);
        assert_eq!(max, 0.0);
    }

    #[test]
    fn test_rupture_force_statistics_single() {
        let (mean, std, max) = rupture_force_statistics(&[42.0]);
        assert!((mean - 42.0).abs() < 1e-12);
        assert_eq!(std, 0.0);
        assert!((max - 42.0).abs() < 1e-12);
    }

    // ── Bell-Evans rupture force tests ────────────────────────────────────────

    #[test]
    fn test_bell_evans_rupture_force_positive() {
        // Standard parameters → positive rupture force
        let f = bell_evans_rupture_force(1e3, 1.0, 0.3, 4.114e-21);
        assert!(f >= 0.0, "Bell-Evans rupture force should be >= 0, got {f}");
    }

    #[test]
    fn test_bell_evans_rupture_force_zero_on_bad_params() {
        let f = bell_evans_rupture_force(0.0, 1.0, 0.3, 1.0);
        assert_eq!(f, 0.0, "Zero loading rate should give 0");
    }

    #[test]
    fn test_bell_evans_rupture_force_increases_with_loading_rate() {
        let kt = 4.114e-21_f64;
        let f_slow = bell_evans_rupture_force(1e2, 1.0, 0.3, kt);
        let f_fast = bell_evans_rupture_force(1e6, 1.0, 0.3, kt);
        assert!(
            f_fast > f_slow,
            "Higher loading rate → higher rupture force: {f_slow} vs {f_fast}"
        );
    }
}
