// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced MD types: Parrinello-Rahman barostat, MSD calculator,
//! running RDF, and adaptive timestep methods.

use std::f64::consts::PI;

use super::core_sim::MdSim;
use super::tracking::{KJ_MOL_PER_ANG3_TO_BAR, normalise_rdf};

/// Smallest ratio `L_new / L_old` that one barostat step may give a box edge.
///
/// Only a blown-up configuration can reach it; the bound keeps the box positive.
const MIN_EDGE_SCALE: f64 = 0.5;

// ---------------------------------------------------------------------------
// ParrinelloRahmanBarostat — extended-Lagrangian NPT barostat
// ---------------------------------------------------------------------------

/// Parrinello-Rahman barostat for an orthorhombic box with isotropic pressure
/// coupling.
///
/// The box edges `L_a` (a = x, y, z) are dynamical variables with their own
/// inertia `W`.  With scaled coordinates `x_ia = L_a s_ia` the extended
/// Lagrangian is
///
/// ```text
/// ℒ = Σ_i ½ m_i L_a² ṡ_ia² − U(x) + ½ W Σ_a L̇_a² − P₀ V
/// ```
///
/// and its equations of motion are
///
/// ```text
/// W L̈_a   = (V / L_a) (P − P₀)                 (box)
/// ẋ_ia    = v_ia + (L̇_a / L_a) x_ia            (affine drift of the coordinates)
/// v̇_ia    = F_ia / m_i − (L̇_a / L_a) v_ia      (peculiar velocities)
/// ```
///
/// where `P` is the instantaneous pressure `(2·KE + W_vir) / (3V)` obtained from
/// the pair virial supplied by the force routine.  The box velocity
/// `L̇_a` ([`ParrinelloRahmanBarostat::box_velocity`]) is **persistent state**: it
/// is accumulated over successive steps, so the box oscillates around the
/// equilibrium volume (damped by the thermostatted particles) instead of
/// relaxing monotonically as under Berendsen scaling.  Dilating the box scales
/// the peculiar velocities with `1/L_a`, so an expanding box cools the system.
///
/// The barostat mass follows the GROMACS convention
/// `W = 3 τ_p² L_max / (4π² β)`, see [`ParrinelloRahmanBarostat::mass`]: for a
/// system whose isothermal compressibility is `β` the volume oscillates with the
/// period `τ_p`.
///
/// Only the scalar pressure is available from the scalar virial, so all three
/// edges are driven by the same pressure (GROMACS `pcoupltype = isotropic`); a
/// cubic box stays cubic.
///
/// Reference: Parrinello & Rahman, J. Appl. Phys. 52, 7182 (1981).
#[derive(Debug, Clone)]
pub struct ParrinelloRahmanBarostat {
    /// Target pressure (bar).
    pub target_pressure: f64,
    /// Period of the box oscillation for a system of compressibility
    /// [`ParrinelloRahmanBarostat::compressibility`] (ps).
    pub tau_p: f64,
    /// Isothermal compressibility of the system (bar⁻¹); sets the barostat mass.
    pub compressibility: f64,
    /// Box-edge velocities \[Lx_dot, Ly_dot, Lz_dot\] (Å ps⁻¹).
    pub box_velocity: [f64; 3],
}

impl ParrinelloRahmanBarostat {
    /// Create a new Parrinello-Rahman barostat with the box at rest.
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

    /// Inertia `W` of the box variables (kJ mol⁻¹ ps² Å⁻², which is 100 amu).
    ///
    /// ```text
    /// W = 3 τ_p² L_max / (4π² β)
    /// ```
    ///
    /// with the compressibility `β` converted to (kJ mol⁻¹ Å⁻³)⁻¹ and `L_max` the
    /// longest box edge.  Linearising `W L̈ = (V / L)(P − P₀)` for a system of
    /// bulk modulus `1/β`, a volume change `δV / V = 3 δL / L` produces
    /// `δP = −3 δL / (β L)`, hence the angular frequency `ω² = 3 L / (β W)` for a
    /// cubic box and the period `2π/ω = τ_p`.
    ///
    /// The box variables are driven by energies (`V · ΔP` in kJ mol⁻¹) and
    /// move in Å and ps, so `W` carries the energy-based unit
    /// kJ mol⁻¹ ps² Å⁻² directly and needs no amu conversion; the particles use
    /// the factor `1 kJ mol⁻¹ = 100 amu Å² ps⁻²` of
    /// [`ACCEL_PER_FORCE_OVER_MASS`](super::ACCEL_PER_FORCE_OVER_MASS) instead.
    pub fn mass(&self, box_lengths: &[f64; 3]) -> f64 {
        let longest_edge = box_lengths.iter().copied().fold(0.0_f64, f64::max);
        let compressibility = self.compressibility * KJ_MOL_PER_ANG3_TO_BAR;
        3.0 * self.tau_p * self.tau_p * longest_edge / (4.0 * PI * PI * compressibility)
    }

    /// Acceleration `L̈_a = V (P − P₀) / (W L_a)` of each box edge (Å ps⁻²) at the
    /// instantaneous pressure `p_current` (bar).
    ///
    /// Positive when the system is above the target pressure.  Returns zeros if the
    /// barostat mass or the box volume is not positive and finite (a degenerate
    /// box has no well-defined pressure force).
    pub fn box_acceleration(&self, p_current: f64, box_lengths: &[f64; 3]) -> [f64; 3] {
        let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
        let mass = self.mass(box_lengths);
        let well_defined = volume.is_finite() && volume > 0.0 && mass.is_finite() && mass > 0.0;
        if !well_defined {
            return [0.0; 3];
        }
        // Pressure excess in kJ mol⁻¹ Å⁻³; V · ΔP is the energy that drives the box.
        let drive = volume * (p_current - self.target_pressure) / KJ_MOL_PER_ANG3_TO_BAR;
        [
            drive / (mass * box_lengths[0]),
            drive / (mass * box_lengths[1]),
            drive / (mass * box_lengths[2]),
        ]
    }

    /// Advance the box edges by one time step `dt` (ps) at the instantaneous
    /// pressure `p_current` (bar), and return the per-axis edge scale factors
    /// `L_new / L_old`.
    ///
    /// The box velocities are kicked with the box acceleration and each edge is
    /// advanced with the new velocity (symplectic Euler).  The caller dilates the
    /// particle coordinates affinely with the returned factors and divides the
    /// peculiar velocities by them (`v ∝ 1/L`); [`ParrinelloRahmanBarostat::advance`]
    /// does so for an [`MdSim`].  An edge that is not positive and finite is left
    /// unchanged (factor 1).
    pub fn advance_edges(
        &mut self,
        box_lengths: &mut [f64; 3],
        p_current: f64,
        dt: f64,
    ) -> [f64; 3] {
        let acceleration = self.box_acceleration(p_current, box_lengths);
        let mut scales = [1.0_f64; 3];
        for axis in 0..3 {
            self.box_velocity[axis] += dt * acceleration[axis];
            let edge = box_lengths[axis];
            if !(edge.is_finite() && edge > 0.0) {
                continue;
            }
            let new_edge = (edge + dt * self.box_velocity[axis]).max(MIN_EDGE_SCALE * edge);
            scales[axis] = new_edge / edge;
            box_lengths[axis] = new_edge;
        }
        scales
    }

    /// Advance the box by one time step `dt` (ps) at the instantaneous pressure
    /// `p_current` (bar).
    ///
    /// The box velocities are kicked with the box acceleration, each edge is
    /// advanced with the new velocity (symplectic Euler, see
    /// [`ParrinelloRahmanBarostat::advance_edges`]), the particle coordinates are
    /// dilated affinely with the edge, and the peculiar velocities are divided by
    /// the same factor (`v ∝ 1/L`).
    pub fn advance(&mut self, sim: &mut MdSim, p_current: f64, dt: f64) {
        let scales = self.advance_edges(&mut sim.config.box_lengths, p_current, dt);
        for pos in sim.state.positions.iter_mut() {
            for (x, scale) in pos.iter_mut().zip(scales.iter()) {
                *x *= scale;
            }
        }
        for vel in sim.state.velocities.iter_mut() {
            for (v, scale) in vel.iter_mut().zip(scales.iter()) {
                *v /= scale;
            }
        }
    }

    /// Apply the barostat to an `MdSim` for one step of length `dt` (ps).
    ///
    /// `virial` is the scalar pair virial `W = Σ_{i<j} r_ij · F_ij` (kJ mol⁻¹) of the
    /// current configuration; the pressure follows from it, the current kinetic
    /// energy and the current box volume (see [`MdSim::compute_pressure`]).
    pub fn apply(&mut self, sim: &mut MdSim, virial: f64, dt: f64) {
        let p_current = sim.compute_pressure(virial);
        self.advance(sim, p_current, dt);
    }
}

// ---------------------------------------------------------------------------
// MsdCalculator — mean square displacement
// ---------------------------------------------------------------------------

/// Computes mean square displacement (MSD) relative to reference positions.
///
/// MSD(t) = <|r(t) - r(0)|^2> averaged over all atoms.
#[derive(Debug, Clone)]
pub struct MsdCalculator {
    /// Reference positions at t = 0 (Å).
    pub reference: Vec<[f64; 3]>,
    /// Accumulated MSD values (Å²) indexed by frame.
    pub msd_values: Vec<f64>,
}

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
        let sum: f64 = positions
            .iter()
            .zip(self.reference.iter())
            .take(n)
            .map(|(pos, refp)| {
                pos.iter()
                    .zip(refp.iter())
                    .map(|(&p, &r)| {
                        let mut d = p - r;
                        d -= box_len * (d / box_len).round();
                        d * d
                    })
                    .sum::<f64>()
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
                let r2: f64 = positions[i]
                    .iter()
                    .zip(positions[j].iter())
                    .map(|(&pi, &pj)| {
                        let mut d = pj - pi;
                        d -= box_len * (d / box_len).round();
                        d * d
                    })
                    .sum();
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
    /// The histogram counts unique pairs `i < j`, so each bin is divided by the
    /// number of unique pairs an ideal gas would place in the shell,
    /// `N (N − 1) / 2 · v_shell / V` per frame; an uncorrelated gas gives
    /// g(r) = 1.  Meaningful for `r_max ≤ L/2` (minimum-image convention).
    pub fn gofr(&self) -> Vec<f64> {
        normalise_rdf(
            &self.histogram,
            self.n_frames,
            self.n_atoms,
            self.box_volume,
            self.dr,
        )
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
    /// `forces_fn(positions, box_lengths) -> (forces, virial)` returns the forces
    /// and the scalar pair virial `W = Σ_{i<j} r_ij · F_ij` (kJ mol⁻¹,
    /// minimum-image separations) of the same configuration; a non-interacting gas
    /// returns `0.0`.  The box equation of motion `W L̈ = (V / L)(P − P₀)` is driven
    /// by the pressure `P = (2·KE + W) / (3V)` of the configuration at the start of
    /// the step; the persistent box velocity is stored in `pr`.  Within a step the
    /// particles advance with velocity Verlet, the box and coordinates are then
    /// dilated by the barostat before the new forces are evaluated, and the
    /// Berendsen thermostat (coupling time `tau_t`) acts last.
    pub fn run_npt_parrinello_rahman(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> (Vec<[f64; 3]>, f64),
        pr: &mut ParrinelloRahmanBarostat,
        tau_t: f64,
        record_every: u64,
    ) -> Vec<(f64, f64, f64, f64)> {
        let (initial_forces, mut virial) =
            forces_fn(&self.state.positions, &self.config.box_lengths);
        self.state.forces = initial_forces;

        let mut records = Vec::new();
        let n_steps = self.config.n_steps;

        for s in 0..n_steps {
            let dt = self.config.dt;
            let pressure = self.compute_pressure(virial);
            virial = self.velocity_verlet_step_with(&forces_fn, |sim| {
                pr.advance(sim, pressure, dt);
            });
            self.apply_berendsen_thermostat(self.config.temperature, tau_t);

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
    pub fn compute_msd(&self, reference_positions: &[[f64; 3]]) -> f64 {
        let calc = MsdCalculator::new(reference_positions.to_vec());
        calc.compute(&self.state.positions)
    }

    /// Accumulate one RDF frame from the current positions into `rdf`.
    ///
    /// Assumes a cubic box with side length equal to `box_lengths[0]`.
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
        self.state.kick_velocities(0.5 * dt_used);
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
        self.state.kick_velocities(0.5 * dt_used);
        self.state.step += 1;
        self.state.time += dt_used;
        dt_used
    }
}
