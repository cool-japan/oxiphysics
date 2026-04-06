//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use oxiphysics_core::math::Vec3;

use super::functions::Thermostat;
#[allow(unused_imports)]
use super::functions::*;

/// Computes kinetic temperature per spatial component for diagnostics.
#[allow(dead_code)]
pub struct KineticDiagnostics;
#[allow(dead_code)]
impl KineticDiagnostics {
    /// Component-wise temperatures: (T_x, T_y, T_z).
    ///
    /// T_alpha = sum_i m_i * v_{i,alpha}^2 / (N * k_B)
    pub fn component_temperatures(atoms: &AtomSet, boltzmann_k: f64) -> (f64, f64, f64) {
        let n = atoms.len();
        if n == 0 {
            return (0.0, 0.0, 0.0);
        }
        let nf = n as f64;
        let mut tx = 0.0;
        let mut ty = 0.0;
        let mut tz = 0.0;
        for i in 0..n {
            let v = atoms.velocities[i];
            let m = atoms.masses[i];
            tx += m * v.x * v.x;
            ty += m * v.y * v.y;
            tz += m * v.z * v.z;
        }
        (
            tx / (nf * boltzmann_k),
            ty / (nf * boltzmann_k),
            tz / (nf * boltzmann_k),
        )
    }
    /// Equipartition ratio: max(T_x, T_y, T_z) / min(T_x, T_y, T_z).
    /// A well-thermalised system should have a ratio close to 1.0.
    pub fn equipartition_ratio(atoms: &AtomSet, boltzmann_k: f64) -> f64 {
        let (tx, ty, tz) = Self::component_temperatures(atoms, boltzmann_k);
        let max_t = tx.max(ty).max(tz);
        let min_t = tx.min(ty).min(tz).max(1e-30);
        max_t / min_t
    }
    /// Total kinetic energy.
    pub fn total_kinetic_energy(atoms: &AtomSet) -> f64 {
        atoms.kinetic_energy()
    }
    /// Centre-of-mass velocity.
    pub fn com_velocity(atoms: &AtomSet) -> Vec3 {
        let mut total_mv = Vec3::zeros();
        let mut total_m = 0.0;
        for i in 0..atoms.len() {
            total_mv += atoms.velocities[i] * atoms.masses[i];
            total_m += atoms.masses[i];
        }
        if total_m < 1e-30 {
            return Vec3::zeros();
        }
        total_mv * (1.0 / total_m)
    }
    /// Remove centre-of-mass drift from all velocities.
    pub fn remove_com_drift(atoms: &mut AtomSet) {
        let v_com = Self::com_velocity(atoms);
        for v in &mut atoms.velocities {
            *v -= v_com;
        }
    }
}
/// Bussi-Donadio-Parrinello CSVR thermostat.
///
/// Rescales velocities stochastically so that the system correctly samples
/// the canonical (NVT) ensemble.  Unlike the simple velocity-rescaling
/// thermostat, CSVR adds a stochastic term derived from the Wiener process
/// to the kinetic energy, ensuring that kinetic energy fluctuations match
/// those of the canonical ensemble.
///
/// # Reference
/// Bussi, Donadio, Parrinello. J. Chem. Phys. 126, 014101 (2007).
#[derive(Debug, Clone)]
pub struct CsvrThermostat {
    /// Coupling time constant tau.
    pub tau: f64,
    /// Internal pseudo-random state (LCG).
    pub(super) rng_state: u64,
}
#[allow(dead_code)]
impl CsvrThermostat {
    /// Create a new CSVR thermostat.
    pub fn new(tau: f64) -> Self {
        Self {
            tau,
            rng_state: 98765,
        }
    }
    /// Create with a specific seed.
    pub fn with_seed(tau: f64, seed: u64) -> Self {
        Self {
            tau,
            rng_state: seed,
        }
    }
    /// Compute the target kinetic energy for the canonical (NVT) ensemble.
    ///
    /// In the canonical ensemble the mean kinetic energy is:
    ///
    /// ```text
    /// `KE` = (N_dof / 2) · k_B · T
    /// ```
    ///
    /// where `N_dof = 3 · N_atoms` is the number of translational degrees of
    /// freedom (assuming fixed centre-of-mass removes 3 DOFs in practice; here
    /// we use the raw `3·N` value as in most implementations).
    ///
    /// CSVR drives the kinetic energy towards this target via stochastic
    /// rescaling; this helper returns the reference value for diagnostics and
    /// as input to the rescaling factor calculation.
    ///
    /// # Arguments
    /// * `n_atoms`     – number of atoms.
    /// * `target_temp` – target temperature T.
    /// * `boltzmann_k` – Boltzmann constant k_B.
    ///
    /// # Returns
    /// Target kinetic energy `(3·N/2) · k_B · T`.
    pub fn compute_kinetic_energy_target(
        &self,
        n_atoms: usize,
        target_temp: f64,
        boltzmann_k: f64,
    ) -> f64 {
        let n_dof = (3 * n_atoms) as f64;
        0.5 * n_dof * boltzmann_k * target_temp
    }
    /// LCG uniform in \[0, 1).
    pub(super) fn next_uniform(&mut self) -> f64 {
        self.rng_state = self
            .rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.rng_state >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Box-Muller normal sample.
    pub(super) fn next_normal(&mut self) -> f64 {
        let u1 = self.next_uniform().max(1e-300);
        let u2 = self.next_uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    /// Sum of n_dof independent chi-squared(1) random variables
    /// approximated via the Wilson-Hilferty transformation for n_dof >= 2.
    fn chi_squared(&mut self, n_dof: f64) -> f64 {
        if n_dof < 2.0 {
            let z = self.next_normal();
            return (z * z).max(0.0);
        }
        let z = self.next_normal();
        let a = 2.0 / (9.0 * n_dof);
        let val = n_dof * (1.0 - a + z * a.sqrt()).powi(3);
        val.max(0.0)
    }
}
/// Stochastic velocity rescaling thermostat (full Bussi et al. implementation).
///
/// Unlike `CsvrThermostat` (which uses a simplified approximation), this
/// implementation uses the exact chi-squared sampling via a sum of M
/// standard normals, which is tractable for small systems.
///
/// # Reference
/// Bussi, Donadio, Parrinello, J. Chem. Phys. 126, 014101 (2007).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SvrThermostat {
    /// Coupling time constant τ (ps).
    pub tau: f64,
    /// Number of degrees of freedom (set automatically from atom count).
    pub n_dof: usize,
    /// LCG RNG state.
    pub(super) rng_state: u64,
}
#[allow(dead_code)]
impl SvrThermostat {
    /// Create a new SVR thermostat.
    ///
    /// # Arguments
    /// * `tau`   – coupling time constant.
    /// * `n_dof` – number of degrees of freedom (usually 3·N).
    pub fn new(tau: f64, n_dof: usize) -> Self {
        Self {
            tau,
            n_dof,
            rng_state: 54321,
        }
    }
    /// Create with explicit seed.
    pub fn with_seed(tau: f64, n_dof: usize, seed: u64) -> Self {
        Self {
            tau,
            n_dof,
            rng_state: seed,
        }
    }
    fn next_uniform(&mut self) -> f64 {
        self.rng_state = self
            .rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.rng_state >> 11) as f64 / (1u64 << 53) as f64
    }
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_uniform().max(1e-300);
        let u2 = self.next_uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    /// Sample a chi-squared(n_dof - 1) random variable by summing squares of
    /// `n_dof - 1` standard normals.  For large n_dof (> 50) we fall back to
    /// the Wilson-Hilferty approximation.
    fn sample_chi2_nm1(&mut self) -> f64 {
        let nf = (self.n_dof as f64) - 1.0;
        if nf <= 0.0 {
            return 0.0;
        }
        if nf > 50.0 {
            let z = self.next_normal();
            let a = 2.0 / (9.0 * nf);
            return (nf * (1.0 - a + z * a.sqrt()).powi(3)).max(0.0);
        }
        let count = nf as usize;
        let mut sum = 0.0;
        for _ in 0..count {
            let z = self.next_normal();
            sum += z * z;
        }
        sum
    }
    /// Compute the rescaling factor for one SVR step.
    ///
    /// The new kinetic energy K_new is sampled from:
    ///   K_new = K_target / n_dof * \[chi^2(n_dof - 1) + (n_1 + sqrt(c*(1-c)) * R)^2\]
    ///
    /// where c = exp(-dt/tau), R ~ N(0,1), and n_1 = sqrt(K_cur / K_target).
    pub fn rescale_factor(&mut self, current_ke: f64, target_ke: f64, dt: f64) -> f64 {
        if current_ke < 1e-30 || target_ke < 1e-30 {
            return 1.0;
        }
        let c = (-dt / self.tau).exp();
        let n_f = self.n_dof as f64;
        let r = self.next_normal();
        let chi2 = self.sample_chi2_nm1();
        let n1 = (current_ke / target_ke).sqrt();
        let term1 = c * n1;
        let term2 = (target_ke / n_f) * ((1.0 - c) * chi2);
        let term3 = 2.0 * (c * (1.0 - c) * current_ke / n_f).max(0.0).sqrt() * r;
        let ke_new = (term1 * term1 * target_ke + term2 + term3).max(0.0);
        if ke_new < 1e-30 {
            return 1.0;
        }
        (ke_new / current_ke).sqrt()
    }
}
/// Berendsen (velocity-rescaling) thermostat.
///
/// Rescales velocities by `sqrt(1 + dt/tau * (T_target/T_current - 1))`.
///
/// This is a weak-coupling thermostat that drives the system toward
/// the target temperature on a time scale `tau`.
#[derive(Debug, Clone)]
pub struct BerendsenThermostat {
    /// Coupling time constant.
    pub tau: f64,
}
impl BerendsenThermostat {
    /// Create a new Berendsen thermostat with coupling constant `tau`.
    pub fn new(tau: f64) -> Self {
        Self { tau }
    }
}
/// Switches between thermostats at specified step counts.
///
/// Useful for equilibration protocols: e.g., use Berendsen for the first
/// 10,000 steps, then switch to Nose-Hoover for production.
#[allow(dead_code)]
pub struct ThermostatSwitcher {
    /// List of (step_threshold, thermostat) pairs.
    /// The active thermostat is the one whose threshold is <= current step.
    pub(super) thermostats: Vec<(usize, Box<dyn Thermostat>)>,
    /// Current step counter.
    pub(super) step: usize,
}
#[allow(dead_code)]
impl ThermostatSwitcher {
    /// Create a new switcher starting with a default thermostat.
    pub fn new(initial: Box<dyn Thermostat>) -> Self {
        Self {
            thermostats: vec![(0, initial)],
            step: 0,
        }
    }
    /// Add a thermostat that activates at the given step.
    pub fn add_phase(&mut self, at_step: usize, thermostat: Box<dyn Thermostat>) {
        self.thermostats.push((at_step, thermostat));
        self.thermostats.sort_by_key(|(s, _)| *s);
    }
    /// Apply the currently active thermostat and advance the step counter.
    pub fn apply_step(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let idx = self
            .thermostats
            .iter()
            .enumerate()
            .rev()
            .find(|(_, (threshold, _))| *threshold <= self.step)
            .map(|(i, _)| i)
            .unwrap_or(0);
        self.thermostats[idx]
            .1
            .apply(atoms, target_temp, dt, boltzmann_k);
        self.step += 1;
    }
    /// Current step count.
    pub fn current_step(&self) -> usize {
        self.step
    }
}
/// Linearly interpolates a target temperature from `T_start` to `T_end`
/// over `ramp_steps` steps, then holds at `T_end`.
///
/// Useful for heating/cooling protocols and simulated annealing.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TemperatureRamp {
    /// Starting temperature.
    pub t_start: f64,
    /// Ending temperature.
    pub t_end: f64,
    /// Number of steps over which to ramp.
    pub ramp_steps: usize,
    /// Current step counter.
    pub(super) step: usize,
}
#[allow(dead_code)]
impl TemperatureRamp {
    /// Create a new temperature ramp.
    pub fn new(t_start: f64, t_end: f64, ramp_steps: usize) -> Self {
        Self {
            t_start,
            t_end,
            ramp_steps,
            step: 0,
        }
    }
    /// Get the target temperature at the current step and advance the counter.
    pub fn next_target_temp(&mut self) -> f64 {
        let t = self.current_target_temp();
        self.step += 1;
        t
    }
    /// Get the target temperature at the current step (without advancing).
    pub fn current_target_temp(&self) -> f64 {
        if self.ramp_steps == 0 || self.step >= self.ramp_steps {
            return self.t_end;
        }
        let frac = self.step as f64 / self.ramp_steps as f64;
        self.t_start + frac * (self.t_end - self.t_start)
    }
    /// Whether the ramp has completed.
    pub fn is_complete(&self) -> bool {
        self.step >= self.ramp_steps
    }
    /// Reset the ramp to step 0.
    pub fn reset(&mut self) {
        self.step = 0;
    }
    /// Current step.
    pub fn current_step(&self) -> usize {
        self.step
    }
    /// Target temperature at a specific step (without changing state).
    pub fn temp_at_step(&self, step: usize) -> f64 {
        if self.ramp_steps == 0 || step >= self.ramp_steps {
            return self.t_end;
        }
        let frac = step as f64 / self.ramp_steps as f64;
        self.t_start + frac * (self.t_end - self.t_start)
    }
}
/// Andersen stochastic collision thermostat.
///
/// At each time step, each atom independently undergoes a stochastic
/// collision with the heat bath with probability `nu * dt`, where `nu` is
/// the collision frequency.  Atoms selected for collision have their
/// velocity drawn from a Maxwell-Boltzmann distribution at the target
/// temperature.
///
/// This thermostat correctly samples the canonical (NVT) ensemble for
/// the configurational degrees of freedom, but it interrupts the true
/// dynamics (velocity autocorrelations are disrupted).
///
/// # Reference
/// Andersen, H. C. J. Chem. Phys. 72, 2384 (1980).
#[derive(Debug, Clone)]
pub struct AndersenThermostat {
    /// Collision frequency ν (per unit time).
    pub nu: f64,
    /// Internal pseudo-random state (simple LCG for determinism).
    pub(super) rng_state: u64,
}
impl AndersenThermostat {
    /// Create a new Andersen thermostat with collision frequency `nu`.
    ///
    /// # Arguments
    /// * `nu` – Collision frequency (probability = nu * dt per atom per step).
    pub fn new(nu: f64) -> Self {
        Self {
            nu,
            rng_state: 12345,
        }
    }
    /// Create a new Andersen thermostat with a specific random seed.
    pub fn with_seed(nu: f64, seed: u64) -> Self {
        Self {
            nu,
            rng_state: seed,
        }
    }
    /// Compute the Poisson collision probability for a single atom over time step `dt`.
    ///
    /// In the Andersen thermostat, atoms undergo stochastic collisions with
    /// the heat bath modelled as a Poisson process with rate `ν`.  The
    /// probability that an atom collides at least once in an interval `dt` is:
    ///
    /// ```text
    /// P(collision) = 1 − exp(−ν · dt)
    /// ```
    ///
    /// For small `ν·dt` this approaches `ν·dt` linearly; for large `ν·dt` it
    /// saturates at 1.0.
    ///
    /// # Arguments
    /// * `dt` – integration time step.
    ///
    /// # Returns
    /// Collision probability in `[0, 1]`.
    pub fn compute_collision_probability(&self, dt: f64) -> f64 {
        1.0 - (-self.nu * dt).exp()
    }
    /// Draw the next pseudo-random `f64` in `[0, 1)` using an LCG.
    pub(super) fn next_uniform(&mut self) -> f64 {
        self.rng_state = self
            .rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.rng_state >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Draw a standard normal sample using the Box-Muller transform.
    pub(super) fn next_normal(&mut self) -> f64 {
        let u1 = self.next_uniform().max(1e-300);
        let u2 = self.next_uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}
/// A no-op thermostat that does nothing (NVE ensemble).
#[derive(Debug, Clone, Default)]
pub struct NoThermostat;
impl NoThermostat {
    /// Create a no-op thermostat.
    pub fn new() -> Self {
        Self
    }
}
/// Nose-Hoover thermostat (extended Lagrangian).
///
/// Introduces an additional degree of freedom `xi` that acts as a
/// friction coefficient on velocities. The coupling mass `q` controls
/// the strength of coupling to the heat bath.
///
/// Equations of motion:
/// - dv/dt = F/m - xi*v
/// - dxi/dt = (2*KE - N_dof*k_B*T_target) / Q
#[derive(Debug, Clone)]
pub struct NoseHooverThermostat {
    /// Coupling mass parameter Q.
    pub q: f64,
    /// Friction variable xi.
    pub xi: f64,
}
impl NoseHooverThermostat {
    /// Create a new Nose-Hoover thermostat with coupling mass `q`.
    pub fn new(q: f64) -> Self {
        Self { q, xi: 0.0 }
    }
}
/// Simple velocity-rescaling thermostat.
///
/// Instantaneously rescales all velocities at every step so that the
/// kinetic energy exactly corresponds to the target temperature.
///
/// Unlike the Berendsen thermostat, this applies an immediate (rather than
/// gradual) correction.  It is simple but does not sample the canonical
/// ensemble correctly; use it only when exact temperature maintenance is
/// required (e.g., equilibration).
#[derive(Debug, Clone, Default)]
pub struct VelocityRescalingThermostat;
impl VelocityRescalingThermostat {
    /// Create a new velocity rescaling thermostat.
    pub fn new() -> Self {
        Self
    }
}
/// Tracks temperature over time for monitoring equilibration.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct TemperatureProfile {
    /// Recorded temperature samples.
    pub(super) samples: Vec<f64>,
    /// Corresponding step numbers.
    pub(super) steps: Vec<usize>,
}
#[allow(dead_code)]
impl TemperatureProfile {
    /// Create a new empty profile.
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
            steps: Vec::new(),
        }
    }
    /// Record a temperature sample.
    pub fn record(&mut self, step: usize, temperature: f64) {
        self.steps.push(step);
        self.samples.push(temperature);
    }
    /// Number of recorded samples.
    pub fn len(&self) -> usize {
        self.samples.len()
    }
    /// Whether the profile is empty.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
    /// Mean temperature over all recorded samples.
    pub fn mean(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        self.samples.iter().sum::<f64>() / self.samples.len() as f64
    }
    /// Standard deviation of temperature.
    pub fn std_dev(&self) -> f64 {
        if self.samples.len() < 2 {
            return 0.0;
        }
        let mean = self.mean();
        let var: f64 = self
            .samples
            .iter()
            .map(|&t| (t - mean).powi(2))
            .sum::<f64>()
            / (self.samples.len() - 1) as f64;
        var.sqrt()
    }
    /// Mean temperature over the last `n` samples.
    pub fn mean_last_n(&self, n: usize) -> f64 {
        let start = if self.samples.len() > n {
            self.samples.len() - n
        } else {
            0
        };
        let slice = &self.samples[start..];
        if slice.is_empty() {
            return 0.0;
        }
        slice.iter().sum::<f64>() / slice.len() as f64
    }
    /// Check if temperature has equilibrated (std dev of last n samples < tolerance).
    pub fn is_equilibrated(&self, n: usize, tolerance: f64) -> bool {
        if self.samples.len() < n {
            return false;
        }
        let start = self.samples.len() - n;
        let slice = &self.samples[start..];
        let mean: f64 = slice.iter().sum::<f64>() / slice.len() as f64;
        let var: f64 = slice.iter().map(|&t| (t - mean).powi(2)).sum::<f64>()
            / (slice.len() - 1).max(1) as f64;
        var.sqrt() < tolerance
    }
    /// Maximum temperature recorded.
    pub fn max_temp(&self) -> f64 {
        self.samples
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum temperature recorded.
    pub fn min_temp(&self) -> f64 {
        self.samples.iter().copied().fold(f64::INFINITY, f64::min)
    }
    /// Latest recorded temperature.
    pub fn latest(&self) -> f64 {
        self.samples.last().copied().unwrap_or(0.0)
    }
}
/// Nosé-Hoover Chain (NHC) thermostat.
///
/// Extends the single Nosé-Hoover oscillator to a chain of `M` heat baths,
/// where each bath thermalises the one before it.  This eliminates the
/// quasi-periodic (non-ergodic) behaviour of the single NH oscillator for
/// small systems.
///
/// # Reference
/// Martyna, Klein, Tuckerman. J. Chem. Phys. 97, 2635 (1992).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NoseHooverChain {
    /// Chain length M.
    pub chain_length: usize,
    /// Coupling masses Q\[k\] (kJ mol⁻¹ ps²).
    pub q: Vec<f64>,
    /// Thermostat positions eta\[k\] (dimensionless).
    pub eta: Vec<f64>,
    /// Thermostat velocities xi\[k\] (ps⁻¹).
    pub xi: Vec<f64>,
    /// Number of degrees of freedom.
    pub n_dof: usize,
}
#[allow(dead_code)]
impl NoseHooverChain {
    /// Create a new NHC thermostat.
    ///
    /// # Arguments
    /// * `chain_length` – number of chain members M (recommended: 3–5).
    /// * `q`            – coupling masses for each chain member.
    /// * `n_dof`        – number of physical degrees of freedom (3·N typically).
    pub fn new(chain_length: usize, q: Vec<f64>, n_dof: usize) -> Self {
        assert_eq!(q.len(), chain_length, "q must have chain_length elements");
        Self {
            chain_length,
            q,
            eta: vec![0.0; chain_length],
            xi: vec![0.0; chain_length],
            n_dof,
        }
    }
    /// Create an NHC thermostat with uniform coupling mass.
    ///
    /// A common choice is Q = n_dof * k_B * T_target * tau^2, where tau
    /// is a characteristic thermostat time.  Here we accept the mass
    /// directly so the caller can pre-compute it.
    pub fn with_uniform_mass(chain_length: usize, mass: f64, n_dof: usize) -> Self {
        let q = vec![mass; chain_length];
        Self::new(chain_length, q, n_dof)
    }
    /// Propagate the NHC by one time step `dt`.
    ///
    /// This uses the Yoshida-Suzuki integrator of order 4 with `n_ys = 3`
    /// Yoshida–Suzuki sub-steps.  The velocities are rescaled and the eta
    /// and xi variables are updated.
    ///
    /// Returns the velocity-scaling factor `s` such that `v_new = s * v_old`.
    pub fn propagate(
        &mut self,
        kinetic_energy: f64,
        target_temp: f64,
        boltzmann_k: f64,
        dt: f64,
    ) -> f64 {
        let w1 = 1.0 / (2.0 - 2.0_f64.powf(1.0 / 3.0));
        let w0 = 1.0 - 2.0 * w1;
        let ws = [w1, w0, w1];
        let mut scale = 1.0_f64;
        let kbt = boltzmann_k * target_temp;
        let g0 = 2.0 * kinetic_energy - (self.n_dof as f64) * kbt;
        for &w in &ws {
            let dt_ys = w * dt;
            let dt2 = dt_ys * 0.5;
            if self.chain_length > 1 {
                let m = self.chain_length - 1;
                let g_m = self.q[m - 1] * self.xi[m - 1] * self.xi[m - 1] - kbt;
                self.xi[m] += g_m / self.q[m] * dt2;
            }
            if self.chain_length > 2 {
                for k in (1..self.chain_length - 1).rev() {
                    let g_k = self.q[k - 1] * self.xi[k - 1] * self.xi[k - 1] - kbt;
                    let xi_kp1 = self.xi[k + 1];
                    let aa = (-xi_kp1 * dt2 * 0.5).exp();
                    self.xi[k] = aa * (aa * self.xi[k] + g_k / self.q[k] * dt2);
                }
            }
            if self.chain_length > 1 {
                let aa = (-self.xi[1] * dt2 * 0.5).exp();
                self.xi[0] = aa * (aa * self.xi[0] + g0 / self.q[0] * dt2);
            } else {
                self.xi[0] += g0 / self.q[0] * dt2;
            }
            let vs = (-self.xi[0] * dt_ys).exp();
            scale *= vs;
            for k in 0..self.chain_length {
                self.eta[k] += self.xi[k] * dt_ys;
            }
            let g0_new = 2.0 * kinetic_energy * scale * scale - (self.n_dof as f64) * kbt;
            if self.chain_length > 1 {
                let aa = (-self.xi[1] * dt2 * 0.5).exp();
                self.xi[0] = aa * (aa * self.xi[0] + g0_new / self.q[0] * dt2);
            } else {
                self.xi[0] += g0_new / self.q[0] * dt2;
            }
            if self.chain_length > 2 {
                for k in 1..self.chain_length - 1 {
                    let g_k = self.q[k - 1] * self.xi[k - 1] * self.xi[k - 1] - kbt;
                    let xi_kp1 = self.xi[k + 1];
                    let aa = (-xi_kp1 * dt2 * 0.5).exp();
                    self.xi[k] = aa * (aa * self.xi[k] + g_k / self.q[k] * dt2);
                }
            }
            if self.chain_length > 1 {
                let m = self.chain_length - 1;
                let g_m = self.q[m - 1] * self.xi[m - 1] * self.xi[m - 1] - kbt;
                self.xi[m] += g_m / self.q[m] * dt2;
            }
        }
        scale
    }
    /// Apply the NHC thermostat to an atom set.
    pub fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        if atoms.is_empty() {
            return;
        }
        let ke = atoms.kinetic_energy();
        let scale = self.propagate(ke, target_temp, boltzmann_k, dt);
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
    /// Conserved energy contribution from the NHC variables.
    ///
    /// H_nhc = Σ_k \[Q_k * xi_k^2 / 2 + k_B * T * eta_k\]
    pub fn conserved_energy(&self, target_temp: f64, boltzmann_k: f64) -> f64 {
        let kbt = boltzmann_k * target_temp;
        self.q
            .iter()
            .zip(self.xi.iter())
            .zip(self.eta.iter())
            .map(|((q, xi), eta)| 0.5 * q * xi * xi + kbt * eta)
            .sum()
    }
    /// Chain length accessor.
    pub fn chain_length(&self) -> usize {
        self.chain_length
    }
}
/// Langevin dynamics thermostat using the BAOAB splitting scheme.
///
/// The Langevin equation of motion is:
///   m·dv/dt = F − γ·m·v + σ·ξ(t)
///
/// where γ is the friction coefficient (ps⁻¹), σ = sqrt(2·γ·k_B·T/m) is
/// the amplitude of the random force, and ξ(t) is Gaussian white noise.
///
/// The BAOAB splitting integrates one step as:
///   B: half-step velocity from force
///   A: half-step position from velocity
///   O: full Ornstein-Uhlenbeck step (friction + noise)
///   A: second half-step position
///   B: second half-step velocity from force
///
/// # Reference
/// Leimkuhler & Matthews, Appl. Math. Res. Express. 2013(1), 34 (2013).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LangevinThermostat {
    /// Friction coefficient γ (ps⁻¹).
    pub gamma: f64,
    /// Internal LCG RNG state.
    pub(super) rng_state: u64,
}
#[allow(dead_code)]
impl LangevinThermostat {
    /// Create a new Langevin thermostat with friction coefficient `gamma`.
    pub fn new(gamma: f64) -> Self {
        Self {
            gamma,
            rng_state: 11111,
        }
    }
    /// Create with a specific random seed for reproducibility.
    pub fn with_seed(gamma: f64, seed: u64) -> Self {
        Self {
            gamma,
            rng_state: seed,
        }
    }
    /// LCG uniform in \[0, 1).
    fn next_uniform(&mut self) -> f64 {
        self.rng_state = self
            .rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.rng_state >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Box-Muller standard normal.
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_uniform().max(1e-300);
        let u2 = self.next_uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    /// Apply the stochastic O-step of BAOAB to atom velocities.
    ///
    /// v_new = c1 * v + c2 * sqrt(k_B * T / m) * xi
    ///
    /// where c1 = exp(-gamma*dt), c2 = sqrt(1 - c1^2).
    pub fn apply_ou_step(
        &mut self,
        atoms: &mut AtomSet,
        target_temp: f64,
        dt: f64,
        boltzmann_k: f64,
    ) {
        let c1 = (-self.gamma * dt).exp();
        let c2 = (1.0 - c1 * c1).max(0.0).sqrt();
        let n = atoms.len();
        for i in 0..n {
            let sigma = (boltzmann_k * target_temp / atoms.masses[i]).sqrt();
            let xi_x = self.next_normal();
            let xi_y = self.next_normal();
            let xi_z = self.next_normal();
            let v = atoms.velocities[i];
            atoms.velocities[i] = Vec3::new(
                c1 * v.x + c2 * sigma * xi_x,
                c1 * v.y + c2 * sigma * xi_y,
                c1 * v.z + c2 * sigma * xi_z,
            );
        }
    }
    /// Random force amplitude σ = sqrt(2 · γ · k_B · T / m).
    pub fn random_force_amplitude(&self, mass: f64, target_temp: f64, boltzmann_k: f64) -> f64 {
        (2.0 * self.gamma * boltzmann_k * target_temp / mass).sqrt()
    }
    /// Compute the friction force and stochastic force acting on a single atom.
    ///
    /// In the Langevin equation of motion:
    /// ```text
    /// F_total = F_conservative + F_friction + F_stochastic
    /// F_friction   = -γ · m · v
    /// F_stochastic = σ · ξ,   σ = sqrt(2 · γ · m · k_B · T / dt)
    /// ```
    ///
    /// # Arguments
    /// * `velocity`    – atom velocity \[vx, vy, vz\].
    /// * `mass`        – atom mass.
    /// * `target_temp` – target temperature T.
    /// * `dt`          – integration time step.
    /// * `boltzmann_k` – Boltzmann constant k_B.
    ///
    /// # Returns
    /// `(friction_force, stochastic_force)` both as `[f64; 3]` arrays.
    pub fn compute_friction_force(
        &mut self,
        velocity: [f64; 3],
        mass: f64,
        target_temp: f64,
        dt: f64,
        boltzmann_k: f64,
    ) -> ([f64; 3], [f64; 3]) {
        let f_fric = [
            -self.gamma * mass * velocity[0],
            -self.gamma * mass * velocity[1],
            -self.gamma * mass * velocity[2],
        ];
        let sigma = if dt > 1e-300 {
            (2.0 * self.gamma * mass * boltzmann_k * target_temp / dt).sqrt()
        } else {
            0.0
        };
        let f_stoch = [
            sigma * self.next_normal(),
            sigma * self.next_normal(),
            sigma * self.next_normal(),
        ];
        (f_fric, f_stoch)
    }
}
