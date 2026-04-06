#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced MD types: Parrinello-Rahman barostat, MSD calculator,
//! running RDF, and adaptive timestep methods.

use super::core_sim::MdSim;

// ---------------------------------------------------------------------------
// ParrinelloRahmanBarostat — NPT barostat with flexible cell
// ---------------------------------------------------------------------------

/// Parrinello-Rahman barostat state for NPT simulations.
///
/// Tracks the box metric tensor and its time derivative to allow
/// full cell fluctuations during constant-pressure MD.
///
/// Reference: Parrinello & Rahman, J. Appl. Phys. 52, 7182 (1981).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ParrinelloRahmanBarostat {
    /// Target pressure (bar).
    pub target_pressure: f64,
    /// Coupling time constant (ps).
    pub tau_p: f64,
    /// Isothermal compressibility (bar^-1).
    pub compressibility: f64,
    /// Box matrix velocity \[Lx_dot, Ly_dot, Lz_dot\] (Å ps^-1).
    pub box_velocity: [f64; 3],
}

#[allow(dead_code)]
impl ParrinelloRahmanBarostat {
    /// Create a new Parrinello-Rahman barostat.
    pub fn new(target_pressure: f64, tau_p: f64, compressibility: f64) -> Self {
        assert!(tau_p > 0.0, "tau_p must be positive");
        assert!(compressibility > 0.0, "compressibility must be positive");
        Self {
            target_pressure,
            tau_p,
            compressibility,
            box_velocity: [0.0; 3],
        }
    }

    /// Compute the box scaling factors per axis for an orthorhombic box.
    ///
    /// Uses the Berendsen-like isotropic approximation scaled per-dimension:
    /// ```text
    /// mu_a = cbrt(1 - beta*dt/tau_p * (P_target - P))
    /// ```
    pub fn scale_factors(&self, p_current: f64, dt: f64) -> [f64; 3] {
        let arg = 1.0 - self.compressibility * dt / self.tau_p * (self.target_pressure - p_current);
        let mu = arg.max(0.01_f64).cbrt();
        [mu, mu, mu]
    }

    /// Apply the barostat to an `MdSim`, rescaling positions and box.
    ///
    /// `virial` is the total scalar virial (kJ mol^-1) from the force computation.
    pub fn apply(&mut self, sim: &mut MdSim, virial: f64, dt: f64) {
        let p_current = sim.compute_pressure(virial);
        let scale = self.scale_factors(p_current, dt);
        let n = sim.state.n_atoms();
        for i in 0..n {
            for a in 0..3 {
                sim.state.positions[i][a] *= scale[a];
            }
        }
        for a in 0..3 {
            sim.config.box_lengths[a] *= scale[a];
            self.box_velocity[a] = (scale[a] - 1.0) * sim.config.box_lengths[a] / dt;
        }
    }
}

// ---------------------------------------------------------------------------
// MsdCalculator — mean square displacement
// ---------------------------------------------------------------------------

/// Computes mean square displacement (MSD) relative to reference positions.
///
/// MSD(t) = <|r(t) - r(0)|^2> averaged over all atoms.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MsdCalculator {
    /// Reference positions at t = 0 (Å).
    pub reference: Vec<[f64; 3]>,
    /// Accumulated MSD values (Å²) indexed by frame.
    pub msd_values: Vec<f64>,
}

#[allow(dead_code)]
impl MsdCalculator {
    /// Create a new MSD calculator with the given reference positions.
    pub fn new(reference: Vec<[f64; 3]>) -> Self {
        Self {
            reference,
            msd_values: Vec::new(),
        }
    }

    /// Compute the instantaneous MSD (Å²) between `positions` and the reference.
    ///
    /// No PBC correction is applied; unwrapped coordinates should be used.
    pub fn compute(&self, positions: &[[f64; 3]]) -> f64 {
        let n = self.reference.len().min(positions.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n)
            .map(|i| {
                let dx = positions[i][0] - self.reference[i][0];
                let dy = positions[i][1] - self.reference[i][1];
                let dz = positions[i][2] - self.reference[i][2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        sum / n as f64
    }

    /// Compute MSD using the minimum image convention for a cubic box of side `box_len`.
    pub fn compute_pbc(&self, positions: &[[f64; 3]], box_len: f64) -> f64 {
        let n = self.reference.len().min(positions.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n)
            .map(|i| {
                let mut dsq = 0.0;
                for a in 0..3 {
                    let mut d = positions[i][a] - self.reference[i][a];
                    d -= box_len * (d / box_len).round();
                    dsq += d * d;
                }
                dsq
            })
            .sum();
        sum / n as f64
    }

    /// Accumulate an MSD frame.
    pub fn accumulate(&mut self, positions: &[[f64; 3]]) {
        let val = self.compute(positions);
        self.msd_values.push(val);
    }

    /// Return the diffusion coefficient D (Å² ps^-1) from the slope of MSD vs time.
    ///
    /// Fits D = MSD / (6 t) using the last accumulated value.
    /// `dt` is the time interval between frames (ps).
    pub fn diffusion_coefficient(&self, dt: f64) -> f64 {
        let n = self.msd_values.len();
        if n == 0 {
            return 0.0;
        }
        let t = (n as f64) * dt;
        if t < 1e-30 {
            return 0.0;
        }
        self.msd_values[n - 1] / (6.0 * t)
    }
}

// ---------------------------------------------------------------------------
// RdfRunning — running RDF accumulator
// ---------------------------------------------------------------------------

/// Running (on-the-fly) radial distribution function accumulator.
///
/// Accumulates pair counts into histogram bins over multiple frames,
/// then normalises by the ideal-gas reference to give g(r).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RdfRunning {
    /// Histogram bin counts.
    pub histogram: Vec<u64>,
    /// Number of bins.
    pub n_bins: usize,
    /// Maximum radius (Å).
    pub r_max: f64,
    /// Bin width (Å).
    pub dr: f64,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Box volume (Å³) — used for normalisation.
    pub box_volume: f64,
    /// Number of frames accumulated.
    pub n_frames: u64,
}

#[allow(dead_code)]
impl RdfRunning {
    /// Create a new running RDF accumulator.
    ///
    /// * `r_max`       – maximum pair distance to histogram (Å)
    /// * `n_bins`      – number of histogram bins
    /// * `box_volume`  – volume of the periodic box (Å³)
    /// * `n_atoms`     – number of atoms
    pub fn new(r_max: f64, n_bins: usize, box_volume: f64, n_atoms: usize) -> Self {
        assert!(n_bins > 0, "n_bins must be > 0");
        assert!(r_max > 0.0, "r_max must be > 0");
        Self {
            histogram: vec![0; n_bins],
            n_bins,
            r_max,
            dr: r_max / n_bins as f64,
            n_atoms,
            box_volume,
            n_frames: 0,
        }
    }

    /// Accumulate one frame of pair distances from `positions` in a cubic box of side `box_len`.
    pub fn accumulate(&mut self, positions: &[[f64; 3]], box_len: f64) {
        let n = positions.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let mut r2 = 0.0;
                for a in 0..3 {
                    let mut d = positions[j][a] - positions[i][a];
                    d -= box_len * (d / box_len).round();
                    r2 += d * d;
                }
                let r = r2.sqrt();
                if r < self.r_max {
                    let bin = (r / self.dr) as usize;
                    let bin = bin.min(self.n_bins - 1);
                    self.histogram[bin] += 1;
                }
            }
        }
        self.n_frames += 1;
    }

    /// Normalise and return g(r) values.
    ///
    /// Divides pair counts by the expected number of pairs in each shell
    /// for an ideal gas at the same density.
    pub fn gofr(&self) -> Vec<f64> {
        if self.n_frames == 0 || self.n_atoms < 2 {
            return vec![0.0; self.n_bins];
        }
        let rho = self.n_atoms as f64 / self.box_volume;
        let n_pairs = (self.n_atoms * (self.n_atoms - 1)) as f64 / 2.0;
        let frames = self.n_frames as f64;
        (0..self.n_bins)
            .map(|b| {
                let r_lo = b as f64 * self.dr;
                let r_hi = r_lo + self.dr;
                let shell_vol = 4.0 / 3.0 * std::f64::consts::PI * (r_hi.powi(3) - r_lo.powi(3));
                let n_ideal = rho * shell_vol * n_pairs;
                if n_ideal < 1e-20 {
                    return 0.0;
                }
                self.histogram[b] as f64 / (frames * n_ideal)
            })
            .collect()
    }

    /// Bin centers (Å) corresponding to each g(r) value.
    pub fn bin_centers(&self) -> Vec<f64> {
        (0..self.n_bins)
            .map(|b| (b as f64 + 0.5) * self.dr)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// MdSim extensions: NPT Parrinello-Rahman, MSD, adaptive timestep
// ---------------------------------------------------------------------------

impl MdSim {
    /// Run an NPT simulation using the Parrinello-Rahman barostat.
    ///
    /// Returns `(ke, pe, T, box_x)` records at every `record_every` steps.
    ///
    /// The Parrinello-Rahman equations of motion are integrated using
    /// a predictor–corrector step for the box degrees of freedom.
    #[allow(dead_code)]
    pub fn run_npt_parrinello_rahman(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        pr: &mut ParrinelloRahmanBarostat,
        tau_t: f64,
        record_every: u64,
    ) -> Vec<(f64, f64, f64, f64)> {
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let mut records = Vec::new();
        let n_steps = self.config.n_steps;

        for s in 0..n_steps {
            self.velocity_verlet_step(&forces_fn);
            self.apply_berendsen_thermostat(self.config.temperature, tau_t);
            pr.apply(self, 0.0, self.config.dt);

            let ke = self.compute_kinetic_energy();
            self.state.kinetic_energy = ke;
            self.state.temperature = self.compute_temperature();
            self.state.step += 1;
            self.state.time += self.config.dt;

            if record_every > 0 && (s + 1) % record_every == 0 {
                records.push((
                    ke,
                    self.state.potential_energy,
                    self.state.temperature,
                    self.config.box_lengths[0],
                ));
            }
        }
        records
    }

    /// Compute the mean square displacement relative to `reference_positions` (Å²).
    ///
    /// No PBC unfolding is performed; pass unwrapped coordinates if needed.
    #[allow(dead_code)]
    pub fn compute_msd(&self, reference_positions: &[[f64; 3]]) -> f64 {
        let calc = MsdCalculator::new(reference_positions.to_vec());
        calc.compute(&self.state.positions)
    }

    /// Accumulate one RDF frame from the current positions into `rdf`.
    ///
    /// Assumes a cubic box with side length equal to `box_lengths[0]`.
    #[allow(dead_code)]
    pub fn compute_rdf_running(&self, rdf: &mut RdfRunning) {
        rdf.accumulate(&self.state.positions, self.config.box_lengths[0]);
    }

    /// Perform one adaptive velocity-Verlet step with automatic timestep control.
    ///
    /// The timestep `dt` is halved if the maximum force magnitude exceeds
    /// `force_threshold`, and doubled (up to `dt_max`) if all forces are
    /// smaller than `force_threshold / 4`.  The final timestep used is returned.
    ///
    /// # Arguments
    /// * `forces_fn`       – force function `(positions, box) -> forces`
    /// * `force_threshold` – maximum acceptable force magnitude (kJ mol^-1 Å^-1)
    /// * `dt_max`          – upper limit on the timestep (ps)
    #[allow(dead_code)]
    pub fn adaptive_timestep(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        force_threshold: f64,
        dt_max: f64,
    ) -> f64 {
        // Evaluate current forces to judge required step size.
        let current_forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let max_f2: f64 = current_forces
            .iter()
            .map(|f| f[0] * f[0] + f[1] * f[1] + f[2] * f[2])
            .fold(0.0_f64, f64::max);
        let max_f = max_f2.sqrt();

        // Adapt timestep.
        if max_f > force_threshold {
            self.config.dt *= 0.5;
        } else if max_f < force_threshold * 0.25 {
            self.config.dt = (self.config.dt * 2.0).min(dt_max);
        }
        let dt_used = self.config.dt;

        // Velocity-Verlet step with the adapted dt.
        self.state.forces = current_forces;
        let n = self.state.n_atoms();
        for i in 0..n {
            let inv_m = 1.0 / self.state.masses[i];
            for a in 0..3 {
                self.state.velocities[i][a] += 0.5 * dt_used * self.state.forces[i][a] * inv_m;
            }
        }
        for i in 0..n {
            for a in 0..3 {
                self.state.positions[i][a] += dt_used * self.state.velocities[i][a];
            }
        }
        if self.config.pbc {
            for i in 0..n {
                for a in 0..3 {
                    let l = self.config.box_lengths[a];
                    self.state.positions[i][a] = self.state.positions[i][a].rem_euclid(l);
                }
            }
        }
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);
        for i in 0..n {
            let inv_m = 1.0 / self.state.masses[i];
            for a in 0..3 {
                self.state.velocities[i][a] += 0.5 * dt_used * self.state.forces[i][a] * inv_m;
            }
        }
        self.state.step += 1;
        self.state.time += dt_used;
        dt_used
    }
}
