//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::simulation::types_sim::*;
use oxiphysics_core::math::Vec3;

/// Phase descriptor for multi-phase SPH simulations.
#[derive(Debug, Clone)]
pub struct SphPhase {
    /// Rest density of the phase \[kg/m³\].
    pub rest_density: f64,
    /// Dynamic viscosity \[Pa·s\].
    pub viscosity: f64,
    /// Speed of sound \[m/s\].
    pub sound_speed: f64,
    /// Human-readable label.
    pub name: String,
}
impl SphPhase {
    /// Preset: liquid water at 20 °C.
    pub fn water() -> Self {
        Self {
            rest_density: 1000.0,
            viscosity: 1.002e-3,
            sound_speed: 1480.0,
            name: "water".into(),
        }
    }
    /// Preset: air at 20 °C, 1 atm.
    pub fn air() -> Self {
        Self {
            rest_density: 1.204,
            viscosity: 1.81e-5,
            sound_speed: 343.0,
            name: "air".into(),
        }
    }
    /// Preset: light mineral oil.
    pub fn oil() -> Self {
        Self {
            rest_density: 870.0,
            viscosity: 3.0e-2,
            sound_speed: 1200.0,
            name: "oil".into(),
        }
    }
    /// Tait equation of state pressure: P = B*((rho/rho0)^gamma - 1).
    pub fn tait_pressure(&self, rho: f64) -> f64 {
        let gamma = 7.0_f64;
        let b = self.rest_density * self.sound_speed * self.sound_speed / gamma;
        b * ((rho / self.rest_density).powf(gamma) - 1.0)
    }
}
/// A lightweight, self-contained SPH simulator using plain `[f64;3]` arrays.
///
/// Maintains its own particle data (positions, velocities, densities, pressures,
/// forces, masses) and runs a simple density → pressure → forces → integrate
/// pipeline each step.
#[derive(Debug, Clone)]
pub struct SphSim {
    /// Simulation configuration.
    pub config: SphSimConfig,
    /// Particle positions \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities \[m/s\].
    pub velocities: Vec<[f64; 3]>,
    /// Current densities \[kg/m³\].
    pub densities: Vec<f64>,
    /// Current pressures \[Pa\].
    pub pressures: Vec<f64>,
    /// Accumulated force per particle \[N\].
    pub forces: Vec<[f64; 3]>,
    /// Particle masses \[kg\].
    pub masses: Vec<f64>,
    /// Simulation state (time, step counter).
    pub state: SphSimState,
}
impl SphSim {
    /// Create a new simulator with `n_particles` particles uniformly distributed
    /// in a unit cube `[0,1)^3` with the given configuration.
    pub fn new(config: SphSimConfig, n_particles: usize) -> Self {
        let n_side = (n_particles as f64).cbrt().ceil() as usize;
        let spacing = 1.0 / n_side.max(1) as f64;
        let mass = config.rest_density * spacing.powi(3);
        let mut positions = Vec::with_capacity(n_particles);
        let mut velocities = Vec::with_capacity(n_particles);
        let mut masses_vec = Vec::with_capacity(n_particles);
        let mut count = 0;
        'outer: for i in 0..n_side {
            for j in 0..n_side {
                for k in 0..n_side {
                    if count >= n_particles {
                        break 'outer;
                    }
                    positions.push([
                        (i as f64 + 0.5) * spacing,
                        (j as f64 + 0.5) * spacing,
                        (k as f64 + 0.5) * spacing,
                    ]);
                    velocities.push([0.0, 0.0, 0.0]);
                    masses_vec.push(mass);
                    count += 1;
                }
            }
        }
        let n = positions.len();
        Self {
            config,
            positions,
            velocities,
            densities: vec![0.0; n],
            pressures: vec![0.0; n],
            forces: vec![[0.0; 3]; n],
            masses: masses_vec,
            state: SphSimState::new(),
        }
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Whether there are no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Run one full step: density → pressure → forces → integrate.
    pub fn step(&mut self) {
        self.compute_density_sum();
        self.compute_pressures_tait();
        self.clear_forces();
        self.apply_gravity();
        self.compute_pressure_forces();
        self.compute_viscosity_forces();
        self.integrate_euler();
        self.state.advance(self.config.dt);
    }
    /// Compute SPH density via kernel summation (cubic spline approximation).
    pub fn compute_density_sum(&mut self) {
        let n = self.len();
        let h = self.config.kernel_radius;
        let h2 = h * h;
        for i in 0..n {
            let mut rho = self.masses[i] * Self::cubic_w(0.0, h);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let d2 = Self::dist2(&self.positions[i], &self.positions[j]);
                if d2 > 4.0 * h2 {
                    continue;
                }
                let r = d2.sqrt();
                rho += self.masses[j] * Self::cubic_w(r, h);
            }
            self.densities[i] = rho;
        }
    }
    /// Compute pressures via Tait EOS: `P = B * ((rho/rho0)^7 - 1)`.
    fn compute_pressures_tait(&mut self) {
        let rho0 = self.config.rest_density;
        let c0 = self.config.estimated_sound_speed().max(1.0);
        let gamma = 7.0_f64;
        let b = rho0 * c0 * c0 / gamma;
        for i in 0..self.len() {
            self.pressures[i] = b * ((self.densities[i] / rho0).powf(gamma) - 1.0);
        }
    }
    /// Compute symmetric SPH pressure forces.
    pub fn compute_pressure_forces(&mut self) {
        let n = self.len();
        let h = self.config.kernel_radius;
        let h2 = h * h;
        for i in 0..n {
            let rhoi = self.densities[i].max(1e-14);
            let pi = self.pressures[i];
            for j in (i + 1)..n {
                let d2 = Self::dist2(&self.positions[i], &self.positions[j]);
                if d2 > 4.0 * h2 {
                    continue;
                }
                let r = d2.sqrt();
                if r < 1e-14 {
                    continue;
                }
                let rhoj = self.densities[j].max(1e-14);
                let pj = self.pressures[j];
                let grad = Self::cubic_grad_w(r, h);
                let factor = -self.masses[j] * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * grad;
                for k in 0..3 {
                    let rhat_k = (self.positions[i][k] - self.positions[j][k]) / r;
                    let fk = factor * rhat_k;
                    self.forces[i][k] += fk * self.masses[i];
                    self.forces[j][k] -= fk * self.masses[j];
                }
            }
        }
    }
    /// Compute artificial viscosity forces (Monaghan style).
    pub fn compute_viscosity_forces(&mut self) {
        let n = self.len();
        let h = self.config.kernel_radius;
        let h2 = h * h;
        let mu = self.config.viscosity;
        for i in 0..n {
            let rhoi = self.densities[i].max(1e-14);
            for j in (i + 1)..n {
                let d2 = Self::dist2(&self.positions[i], &self.positions[j]);
                if d2 > 4.0 * h2 || d2 < 1e-28 {
                    continue;
                }
                let r = d2.sqrt();
                let rhoj = self.densities[j].max(1e-14);
                let rho_avg = 0.5 * (rhoi + rhoj);
                let mut vr_dot = 0.0_f64;
                for k in 0..3 {
                    vr_dot += (self.velocities[i][k] - self.velocities[j][k])
                        * (self.positions[i][k] - self.positions[j][k]);
                }
                if vr_dot >= 0.0 {
                    continue;
                }
                let nu = h * vr_dot / (d2 + 0.01 * h2);
                let pi_visc = -mu * nu / rho_avg;
                let grad = Self::cubic_grad_w(r, h);
                for k in 0..3 {
                    let rhat_k = (self.positions[i][k] - self.positions[j][k]) / r;
                    let fk = -self.masses[j] * pi_visc * grad * rhat_k;
                    self.forces[i][k] += fk * self.masses[i];
                    self.forces[j][k] -= fk * self.masses[j];
                }
            }
        }
    }
    /// Clear all forces to zero.
    fn clear_forces(&mut self) {
        for f in &mut self.forces {
            *f = [0.0; 3];
        }
    }
    /// Add gravity to forces.
    fn apply_gravity(&mut self) {
        let g = self.config.gravity;
        for i in 0..self.len() {
            let m = self.masses[i];
            self.forces[i][0] += m * g[0];
            self.forces[i][1] += m * g[1];
            self.forces[i][2] += m * g[2];
        }
    }
    /// Euler integration: `v += (F/m)*dt`, `r += v*dt`.
    fn integrate_euler(&mut self) {
        let dt = self.config.dt;
        for i in 0..self.len() {
            let m = self.masses[i].max(1e-30);
            for k in 0..3 {
                self.velocities[i][k] += (self.forces[i][k] / m) * dt;
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
        }
    }
    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }
    /// Simplified cubic spline kernel W(r, h).
    pub(crate) fn cubic_w(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q > 2.0 {
            return 0.0;
        }
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
        if q <= 1.0 {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else {
            sigma * 0.25 * (2.0 - q).powi(3)
        }
    }
    /// Gradient magnitude of the cubic spline kernel dW/dr.
    pub(crate) fn cubic_grad_w(r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..=2.0).contains(&q) {
            return 0.0;
        }
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h * h);
        if q <= 1.0 {
            sigma * (-3.0 * q + 2.25 * q * q)
        } else {
            sigma * (-0.75 * (2.0 - q).powi(2))
        }
    }
    /// Squared distance between two points.
    #[inline]
    pub(crate) fn dist2(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
    }
}
impl SphSim {
    /// Compute Continuum Surface Force (CSF) surface tension forces and
    /// accumulate them into `self.forces`.
    ///
    /// Uses the Morris (2000) colour-function formulation:
    ///   F_st_i = sigma * kappa_i * n_i
    /// where `kappa` is the interface curvature and `n_i` the unit outward
    /// normal, both estimated from the colour-function gradient (here we use
    /// density as a proxy).
    pub fn compute_surface_tension_forces(&mut self) {
        let n = self.len();
        let h = self.config.kernel_radius;
        let h2 = h * h;
        let sigma = self.config.surface_tension_coeff;
        if sigma.abs() < 1e-30 {
            return;
        }
        let mut grad_c = vec![[0.0_f64; 3]; n];
        for (i, gc_i) in grad_c.iter_mut().enumerate().take(n) {
            let rhoi = self.densities[i].max(1e-14);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let d2 = Self::dist2(&self.positions[i], &self.positions[j]);
                if d2 > 4.0 * h2 || d2 < 1e-28 {
                    continue;
                }
                let r = d2.sqrt();
                let rhoj = self.densities[j].max(1e-14);
                let ci = rhoi / self.config.rest_density;
                let cj = rhoj / self.config.rest_density;
                let grad_w = Self::cubic_grad_w(r, h);
                let factor = self.masses[j] / rhoj * (cj - ci) * grad_w / r;
                for (k, gck) in gc_i.iter_mut().enumerate().take(3) {
                    *gck += factor * (self.positions[j][k] - self.positions[i][k]);
                }
            }
        }
        for i in 0..n {
            let gc_mag2 = grad_c[i][0] * grad_c[i][0]
                + grad_c[i][1] * grad_c[i][1]
                + grad_c[i][2] * grad_c[i][2];
            if gc_mag2 < 1e-20 {
                continue;
            }
            let gc_mag = gc_mag2.sqrt();
            let n_hat = [
                grad_c[i][0] / gc_mag,
                grad_c[i][1] / gc_mag,
                grad_c[i][2] / gc_mag,
            ];
            let rhoi = self.densities[i].max(1e-14);
            let mut kappa = 0.0_f64;
            for (j, gc_j) in grad_c.iter().enumerate().take(n) {
                if j == i {
                    continue;
                }
                let d2 = Self::dist2(&self.positions[i], &self.positions[j]);
                if d2 > 4.0 * h2 || d2 < 1e-28 {
                    continue;
                }
                let r = d2.sqrt();
                let rhoj = self.densities[j].max(1e-14);
                let gc_j_mag2 = gc_j[0] * gc_j[0] + gc_j[1] * gc_j[1] + gc_j[2] * gc_j[2];
                if gc_j_mag2 < 1e-20 {
                    continue;
                }
                let gc_j_mag = gc_j_mag2.sqrt();
                let n_hat_j = [gc_j[0] / gc_j_mag, gc_j[1] / gc_j_mag, gc_j[2] / gc_j_mag];
                let grad_w = Self::cubic_grad_w(r, h);
                let mut dn_dot_r = 0.0_f64;
                for k in 0..3 {
                    dn_dot_r +=
                        (n_hat_j[k] - n_hat[k]) * (self.positions[j][k] - self.positions[i][k]);
                }
                kappa -= self.masses[j] / rhoj * dn_dot_r * grad_w / r;
            }
            let m_i = self.masses[i];
            let force_scale = sigma * kappa * gc_mag * m_i / rhoi;
            for (k, nk) in n_hat.iter().enumerate() {
                self.forces[i][k] += force_scale * nk;
            }
        }
    }
    /// Run one full step including surface tension:
    /// density → pressure → forces → surface tension → integrate.
    pub fn step_with_surface_tension(&mut self) {
        self.compute_density_sum();
        self.compute_pressures_tait();
        self.clear_forces();
        self.apply_gravity();
        self.compute_pressure_forces();
        self.compute_viscosity_forces();
        self.compute_surface_tension_forces();
        self.integrate_euler();
        self.state.advance(self.config.dt);
    }
}
impl SphSim {
    /// Compute the gravitational potential energy Σ m_i * g_y * y_i.
    pub fn potential_energy(&self) -> f64 {
        let g_mag = (self.config.gravity[0].powi(2)
            + self.config.gravity[1].powi(2)
            + self.config.gravity[2].powi(2))
        .sqrt();
        self.positions
            .iter()
            .zip(self.masses.iter())
            .map(|(p, &m)| m * g_mag * p[1])
            .sum()
    }
    /// Total mechanical energy: KE + PE.
    pub fn total_mechanical_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }
    /// Total linear momentum.
    pub fn linear_momentum(&self) -> [f64; 3] {
        let mut p = [0.0_f64; 3];
        for (v, &m) in self.velocities.iter().zip(self.masses.iter()) {
            p[0] += m * v[0];
            p[1] += m * v[1];
            p[2] += m * v[2];
        }
        p
    }
    /// Mean particle spacing: `V^(1/3) / n^(1/3)` where V = box volume.
    ///
    /// Assumes a unit box `[0,1)^3`.
    pub fn mean_particle_spacing(&self) -> f64 {
        let n = self.len();
        if n == 0 {
            return 0.0;
        }
        (1.0 / n as f64).cbrt()
    }
    /// Return the positions as a flat `Vec`f64` in `\[x0, y0, z0, x1, y1, z1, ...\]` order.
    pub fn flat_positions(&self) -> Vec<f64> {
        self.positions
            .iter()
            .flat_map(|p| [p[0], p[1], p[2]])
            .collect()
    }
    /// Compute the maximum inter-particle distance (brute-force O(n²)).
    pub fn max_inter_particle_distance(&self) -> f64 {
        let n = self.len();
        let mut max_d = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d > max_d {
                    max_d = d;
                }
            }
        }
        max_d
    }
    /// Apply a reflective floor boundary at `y_floor` with coefficient of
    /// restitution `e`: particles below the floor are reflected upward.
    pub fn apply_floor_reflection(&mut self, y_floor: f64, restitution: f64) {
        for i in 0..self.len() {
            if self.positions[i][1] < y_floor {
                self.positions[i][1] = 2.0 * y_floor - self.positions[i][1];
                self.velocities[i][1] *= -restitution;
            }
        }
    }
    /// Scale all particle masses by `factor`.
    pub fn scale_masses(&mut self, factor: f64) {
        for m in &mut self.masses {
            *m *= factor;
        }
    }
    /// Mean pressure across all particles.
    pub fn mean_pressure(&self) -> f64 {
        let n = self.len();
        if n == 0 {
            return 0.0;
        }
        self.pressures.iter().sum::<f64>() / n as f64
    }
    /// Number of particles with negative pressure.
    pub fn count_negative_pressure(&self) -> usize {
        self.pressures.iter().filter(|&&p| p < 0.0).count()
    }
}
/// Tracks kinetic and potential energy samples over a simulation run.
#[derive(Debug, Clone, Default)]
pub struct EnergyTracker {
    /// Sampled kinetic energies.
    pub kinetic: Vec<f64>,
    /// Sampled potential energies.
    pub potential: Vec<f64>,
    /// Sampled simulation times.
    pub times: Vec<f64>,
}
impl EnergyTracker {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record one energy sample at time `t`.
    pub fn record(&mut self, t: f64, ke: f64, pe: f64) {
        self.times.push(t);
        self.kinetic.push(ke);
        self.potential.push(pe);
    }
    /// Total mechanical energies as `KE + PE`.
    pub fn total_energies(&self) -> Vec<f64> {
        self.kinetic
            .iter()
            .zip(self.potential.iter())
            .map(|(k, p)| k + p)
            .collect()
    }
    /// Number of recorded samples.
    pub fn len(&self) -> usize {
        self.times.len()
    }
    /// True if no samples have been recorded.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
    /// Mean kinetic energy over all samples.
    pub fn mean_kinetic(&self) -> f64 {
        if self.kinetic.is_empty() {
            return 0.0;
        }
        self.kinetic.iter().sum::<f64>() / self.kinetic.len() as f64
    }
    /// Maximum kinetic energy over all samples.
    pub fn max_kinetic(&self) -> f64 {
        self.kinetic.iter().copied().fold(0.0_f64, f64::max)
    }
}
/// Fluent builder for [`WcSphSim`].
#[derive(Debug, Clone)]
pub struct WcSphSimBuilder {
    pub(super) h: f64,
    pub(super) rho0: f64,
    pub(super) c0: f64,
    pub(super) gamma: f64,
    pub(super) mu: f64,
    pub(super) sigma: f64,
    pub(super) gravity: [f64; 3],
}
impl WcSphSimBuilder {
    /// Start with water-like defaults (ρ₀ = 1000, c₀ = 100, γ = 7).
    pub fn water() -> Self {
        Self {
            h: 0.1,
            rho0: 1000.0,
            c0: 100.0,
            gamma: 7.0,
            mu: 1e-3,
            sigma: 0.0,
            gravity: [0.0, -9.81, 0.0],
        }
    }
    /// Set the smoothing length.
    pub fn with_h(mut self, h: f64) -> Self {
        self.h = h;
        self
    }
    /// Set the rest density.
    pub fn with_rho0(mut self, rho0: f64) -> Self {
        self.rho0 = rho0;
        self
    }
    /// Set the speed of sound.
    pub fn with_c0(mut self, c0: f64) -> Self {
        self.c0 = c0;
        self
    }
    /// Set the dynamic viscosity.
    pub fn with_mu(mut self, mu: f64) -> Self {
        self.mu = mu;
        self
    }
    /// Set surface tension coefficient.
    pub fn with_sigma(mut self, sigma: f64) -> Self {
        self.sigma = sigma;
        self
    }
    /// Set the gravity vector.
    pub fn with_gravity(mut self, g: [f64; 3]) -> Self {
        self.gravity = g;
        self
    }
    /// Build the [`WcSphSim`].
    pub fn build(self) -> WcSphSim {
        WcSphSim::new(
            self.h,
            self.rho0,
            self.c0,
            self.gamma,
            self.mu,
            self.sigma,
            self.gravity,
        )
    }
}
/// A lightweight harness for benchmarking [`SphSim`] step performance.
#[derive(Debug, Clone)]
pub struct BenchmarkHarness {
    /// Number of warm-up steps (not counted in timing).
    pub warmup_steps: usize,
    /// Number of timed steps.
    pub timed_steps: usize,
    /// Total step-seconds accumulated during timed phase.
    pub step_ns_total: u64,
    /// Number of timing runs performed.
    pub runs: usize,
}
impl BenchmarkHarness {
    /// Create a new harness.
    pub fn new(warmup_steps: usize, timed_steps: usize) -> Self {
        Self {
            warmup_steps,
            timed_steps,
            step_ns_total: 0,
            runs: 0,
        }
    }
    /// Run warm-up + timed steps on `sim`; accumulates timing counts.
    ///
    /// Since we avoid `std::time` in library code, this version just counts
    /// steps and does not measure wall-clock time (timing must be done by the
    /// caller if needed).
    pub fn run(&mut self, sim: &mut SphSim) {
        for _ in 0..self.warmup_steps {
            sim.step();
        }
        for _ in 0..self.timed_steps {
            sim.step();
            self.step_ns_total += 1;
        }
        self.runs += 1;
    }
    /// Total steps executed across all runs.
    pub fn total_steps(&self) -> usize {
        self.runs * (self.warmup_steps + self.timed_steps)
    }
}
/// A lightweight snapshot of simulation state for output callbacks.
#[derive(Debug, Clone)]
pub struct SimSnapshot {
    /// Current simulation time.
    pub time: f64,
    /// Step number.
    pub step: u64,
    /// Particle positions at the time of snapshot.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities at the time of snapshot.
    pub velocities: Vec<[f64; 3]>,
    /// Particle densities at the time of snapshot.
    pub densities: Vec<f64>,
}
impl SimSnapshot {
    /// Create a snapshot from a [`SphSim`].
    pub fn from_sim(sim: &SphSim) -> Self {
        Self {
            time: sim.state.time,
            step: sim.state.step,
            positions: sim.positions.clone(),
            velocities: sim.velocities.clone(),
            densities: sim.densities.clone(),
        }
    }
    /// Number of particles in the snapshot.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// True if no particles are recorded.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Compute the centre of mass position from snapshot data and equal masses.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let n = self.positions.len();
        if n == 0 {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for p in &self.positions {
            com[0] += p[0];
            com[1] += p[1];
            com[2] += p[2];
        }
        let inv_n = 1.0 / n as f64;
        [com[0] * inv_n, com[1] * inv_n, com[2] * inv_n]
    }
}
/// Runs a [`SphSim`] for a prescribed number of steps, collecting statistics
/// and snapshots at user-defined intervals.
#[derive(Debug)]
pub struct SphRunner {
    /// The embedded simulation.
    pub sim: SphSim,
    /// Collected statistics.
    pub stats: SimulationStats,
    /// Snapshots taken at each output interval.
    pub snapshots: Vec<SimSnapshot>,
    /// Output every this many steps (0 = disabled).
    pub snapshot_every: u64,
}
impl SphRunner {
    /// Create a runner wrapping an existing [`SphSim`].
    pub fn new(sim: SphSim, snapshot_every: u64) -> Self {
        Self {
            sim,
            stats: SimulationStats::new(),
            snapshots: Vec::new(),
            snapshot_every,
        }
    }
    /// Compute maximum particle speed.
    fn max_speed(sim: &SphSim) -> f64 {
        sim.velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }
    /// Run for `n_steps` steps.
    pub fn run_steps(&mut self, n_steps: u64) {
        for _ in 0..n_steps {
            self.sim.step();
            let dt = self.sim.config.dt;
            let ke = self.sim.kinetic_energy();
            let ms = Self::max_speed(&self.sim);
            self.stats.record_step(dt, ke, ms);
            if self.snapshot_every > 0 && self.sim.state.step.is_multiple_of(self.snapshot_every) {
                self.snapshots.push(SimSnapshot::from_sim(&self.sim));
            }
        }
    }
    /// Total number of snapshots collected.
    pub fn num_snapshots(&self) -> usize {
        self.snapshots.len()
    }
}
/// Parameters for the SPH simulation.
#[derive(Debug, Clone)]
pub struct SphSimulationParams {
    /// Gravity vector.
    pub gravity: Vec3,
    /// Smoothing length.
    pub smoothing_length: f64,
    /// Numerical speed of sound (for CFL and Tait EOS).
    pub sound_speed: f64,
    /// Rest density (kg/m³).
    pub rest_density: f64,
    /// Kinematic viscosity.
    pub viscosity: f64,
    /// Solver type.
    pub solver_type: SolverType,
    /// Boundary penalty stiffness.
    pub boundary_stiffness: f64,
    /// Boundary penalty damping.
    pub boundary_damping: f64,
    /// Minimum timestep.
    pub min_dt: f64,
    /// Maximum timestep.
    pub max_dt: f64,
}
/// Run `n_replicas` independent SPH simulations with small random velocity
/// perturbations and accumulate statistical summaries.
///
/// This is useful for studying convergence, statistical reproducibility, and
/// sensitivity to initial conditions.
#[derive(Debug, Clone)]
pub struct EnsembleSimulation {
    /// Number of replicas.
    pub n_replicas: usize,
    /// Number of steps to advance each replica.
    pub n_steps: usize,
    /// Shared base configuration.
    pub config: SphSimConfig,
    /// Initial number of particles per replica.
    pub n_particles: usize,
    /// Mean kinetic energy measured at the end of each replica run.
    pub final_ke: Vec<f64>,
    /// Final total momenta (x-component) per replica.
    pub final_mom_x: Vec<f64>,
}
impl EnsembleSimulation {
    /// Create a new ensemble runner.
    pub fn new(
        config: SphSimConfig,
        n_particles: usize,
        n_replicas: usize,
        n_steps: usize,
    ) -> Self {
        Self {
            n_replicas,
            n_steps,
            config,
            n_particles,
            final_ke: Vec::new(),
            final_mom_x: Vec::new(),
        }
    }
    /// Run all replicas.  The velocity perturbation amplitude `v_perturb`
    /// is added to a deterministic pseudo-random offset per particle index.
    pub fn run(&mut self, v_perturb: f64) {
        self.final_ke.clear();
        self.final_mom_x.clear();
        for replica in 0..self.n_replicas {
            let mut sim = SphSim::new(self.config.clone(), self.n_particles);
            for (i, v) in sim.velocities.iter_mut().enumerate() {
                let t = ((i + 1 + replica * 137) as f64).sin();
                v[0] += v_perturb * t;
            }
            for _ in 0..self.n_steps {
                sim.step();
            }
            let ke: f64 = sim
                .velocities
                .iter()
                .zip(sim.masses.iter())
                .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
                .sum();
            let mom_x: f64 = sim
                .velocities
                .iter()
                .zip(sim.masses.iter())
                .map(|(v, &m)| m * v[0])
                .sum();
            self.final_ke.push(ke);
            self.final_mom_x.push(mom_x);
        }
    }
    /// Mean kinetic energy across replicas.
    pub fn mean_ke(&self) -> f64 {
        if self.final_ke.is_empty() {
            return 0.0;
        }
        self.final_ke.iter().sum::<f64>() / self.final_ke.len() as f64
    }
    /// Standard deviation of final kinetic energy across replicas.
    pub fn std_ke(&self) -> f64 {
        let n = self.final_ke.len();
        if n < 2 {
            return 0.0;
        }
        let mean = self.mean_ke();
        let var = self
            .final_ke
            .iter()
            .map(|&e| (e - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        var.sqrt()
    }
    /// Whether all replicas produced finite kinetic energies.
    pub fn all_finite(&self) -> bool {
        self.final_ke.iter().all(|e| e.is_finite())
    }
}
/// A particle set that carries per-particle phase identifiers, enabling
/// multi-phase SPH simulations.
#[derive(Debug, Clone)]
pub struct MultiPhaseParticleSet {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// Per-particle densities.
    pub densities: Vec<f64>,
    /// Per-particle pressures.
    pub pressures: Vec<f64>,
    /// Phase index for each particle.
    pub phase_ids: Vec<usize>,
    /// Registered phases.
    pub phases: Vec<SphPhase>,
}
impl MultiPhaseParticleSet {
    /// Create an empty multi-phase particle set.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            masses: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            phase_ids: Vec::new(),
            phases: Vec::new(),
        }
    }
    /// Register a phase and return its index.
    pub fn add_phase(&mut self, phase: SphPhase) -> usize {
        let idx = self.phases.len();
        self.phases.push(phase);
        idx
    }
    /// Add a particle belonging to `phase_id`.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64, phase_id: usize) {
        let rho0 = if phase_id < self.phases.len() {
            self.phases[phase_id].rest_density
        } else {
            1000.0
        };
        self.positions.push(pos);
        self.velocities.push(vel);
        self.masses.push(mass);
        self.densities.push(rho0);
        self.pressures.push(0.0);
        self.phase_ids.push(phase_id);
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// True if there are no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Count particles per phase. Returns a vector of length `phases.len()`.
    pub fn count_by_phase(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.phases.len()];
        for &pid in &self.phase_ids {
            if pid < counts.len() {
                counts[pid] += 1;
            }
        }
        counts
    }
    /// Update pressures using the Tait EOS for each particle's phase.
    pub fn update_pressures(&mut self) {
        for i in 0..self.len() {
            let pid = self.phase_ids[i];
            if pid < self.phases.len() {
                self.pressures[i] = self.phases[pid].tait_pressure(self.densities[i]);
            }
        }
    }
    /// Compute total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }
    /// Compute the centre of mass.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for i in 0..self.len() {
            let m = self.masses[i];
            com[0] += m * self.positions[i][0];
            com[1] += m * self.positions[i][1];
            com[2] += m * self.positions[i][2];
        }
        let inv = 1.0 / total_mass;
        [com[0] * inv, com[1] * inv, com[2] * inv]
    }
}
/// Lightweight configuration for array-based SPH mini-simulations.
///
/// This is independent of the higher-level [`SphSimulationParams`] and is
/// intended for quick prototyping and educational use with [`SphSim`].
#[derive(Debug, Clone)]
pub struct SphSimConfig {
    /// Smoothing / kernel radius (m).
    pub kernel_radius: f64,
    /// Rest density (kg/m³).
    pub rest_density: f64,
    /// Dynamic viscosity coefficient (Pa·s).
    pub viscosity: f64,
    /// Surface tension coefficient (N/m).
    pub surface_tension_coeff: f64,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
    /// Time step size (s).
    pub dt: f64,
}
impl SphSimConfig {
    /// Create a default water-like configuration.
    pub fn default_water() -> Self {
        Self {
            kernel_radius: 0.1,
            rest_density: 1000.0,
            viscosity: 0.01,
            surface_tension_coeff: 0.0728,
            gravity: [0.0, -9.81, 0.0],
            dt: 0.001,
        }
    }
    /// Speed of sound estimate: `c0 = 10 * sqrt(2 * g * H)` with a default
    /// domain height H = 1 m.
    pub fn estimated_sound_speed(&self) -> f64 {
        let g_mag =
            (self.gravity[0].powi(2) + self.gravity[1].powi(2) + self.gravity[2].powi(2)).sqrt();
        10.0 * (2.0 * g_mag * 1.0).sqrt()
    }
    /// CFL-limited time step: `dt = 0.4 * h / c0`.
    pub fn cfl_dt(&self) -> f64 {
        let c0 = self.estimated_sound_speed().max(1.0);
        0.4 * self.kernel_radius / c0
    }
}
/// Choice of SPH pressure solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverType {
    /// Weakly compressible SPH with Tait equation of state.
    Wcsph,
    /// Divergence-free SPH with iterative pressure solve.
    Dfsph,
}
/// Tracks the current time and step number for a simulation.
#[derive(Debug, Clone)]
pub struct SphSimState {
    /// Current simulation time (s).
    pub time: f64,
    /// Current step number.
    pub step: u64,
}
impl SphSimState {
    /// Create a fresh state at t = 0, step = 0.
    pub fn new() -> Self {
        Self { time: 0.0, step: 0 }
    }
    /// Advance the state by one step of size `dt`.
    pub fn advance(&mut self, dt: f64) {
        self.time += dt;
        self.step += 1;
    }
}
/// Adaptive timestep controller that tracks step history and adjusts `dt`.
#[derive(Debug, Clone)]
pub struct AdaptiveTimestep {
    /// Current timestep.
    pub dt: f64,
    /// Minimum allowed timestep.
    pub dt_min: f64,
    /// Maximum allowed timestep.
    pub dt_max: f64,
    /// CFL number target (default 0.4).
    pub cfl_target: f64,
    /// Viscous safety factor (default 0.125).
    pub viscous_factor: f64,
    /// History of accepted steps.
    pub history: Vec<f64>,
}
impl AdaptiveTimestep {
    /// Create a new adaptive timestep controller.
    pub fn new(dt_initial: f64, dt_min: f64, dt_max: f64) -> Self {
        Self {
            dt: dt_initial,
            dt_min,
            dt_max,
            cfl_target: 0.4,
            viscous_factor: 0.125,
            history: Vec::new(),
        }
    }
    /// Compute the CFL-limited timestep given the maximum velocity magnitude and kernel radius.
    pub fn cfl_limit(&self, max_velocity: f64, h: f64, sound_speed: f64) -> f64 {
        let u_max = max_velocity + sound_speed;
        if u_max < 1e-14 {
            return self.dt_max;
        }
        self.cfl_target * h / u_max
    }
    /// Compute the viscosity-limited timestep.
    pub fn viscous_limit(&self, h: f64, kinematic_viscosity: f64) -> f64 {
        if kinematic_viscosity < 1e-30 {
            return self.dt_max;
        }
        self.viscous_factor * h * h / kinematic_viscosity
    }
    /// Select the next timestep as the minimum of CFL, viscous, and max limits.
    pub fn select(&mut self, max_velocity: f64, h: f64, sound_speed: f64, nu: f64) -> f64 {
        let dt_cfl = self.cfl_limit(max_velocity, h, sound_speed);
        let dt_visc = self.viscous_limit(h, nu);
        let dt_new = dt_cfl.min(dt_visc).clamp(self.dt_min, self.dt_max);
        self.dt = dt_new;
        self.history.push(dt_new);
        dt_new
    }
    /// Number of steps taken so far.
    pub fn steps_taken(&self) -> usize {
        self.history.len()
    }
    /// Mean timestep over all recorded steps.
    pub fn mean_dt(&self) -> f64 {
        if self.history.is_empty() {
            return self.dt;
        }
        self.history.iter().sum::<f64>() / self.history.len() as f64
    }
}
/// Statistics collected during an SPH simulation run.
#[derive(Debug, Clone, Default)]
pub struct SimulationStats {
    /// Number of timesteps executed.
    pub num_steps: u64,
    /// Total simulated time elapsed.
    pub total_time: f64,
    /// Minimum timestep ever used.
    pub min_dt: f64,
    /// Maximum timestep ever used.
    pub max_dt: f64,
    /// Last recorded kinetic energy.
    pub last_kinetic_energy: f64,
    /// Last recorded maximum particle speed.
    pub last_max_speed: f64,
}
impl SimulationStats {
    /// Create a fresh statistics record.
    pub fn new() -> Self {
        Self {
            min_dt: f64::INFINITY,
            ..Default::default()
        }
    }
    /// Record one accepted timestep and update bookkeeping.
    pub fn record_step(&mut self, dt: f64, ke: f64, max_speed: f64) {
        self.num_steps += 1;
        self.total_time += dt;
        if dt < self.min_dt {
            self.min_dt = dt;
        }
        if dt > self.max_dt {
            self.max_dt = dt;
        }
        self.last_kinetic_energy = ke;
        self.last_max_speed = max_speed;
    }
    /// Mean simulated time per step.
    pub fn mean_dt(&self) -> f64 {
        if self.num_steps == 0 {
            return 0.0;
        }
        self.total_time / self.num_steps as f64
    }
}
