//! Simulation types: WcSphSim, SphSimulation, TwoPhaseSimState, PeriodicSphSim, CflController.

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use crate::boundary_sph::BoundarySet;
use crate::dfsph::{DfsphParams, DfsphSolver};
use crate::kernel::SphKernel;
use crate::neighbor::SpatialHash;
use crate::particle::ParticleSet;
use crate::timestep;
use crate::wcsph::WcsphParams;

#[allow(unused_imports)]
use super::types::{SolverType, SphSim, SphSimulationParams};

/// Weakly-compressible SPH simulation using `SphParticleSet` (SoA layout).
///
/// Implements the full pipeline:
///   density_step → pressure_step → force_step → integrate
///
/// Physics:
/// - WCSPH Tait EOS: `P = B * ((ρ/ρ₀)^γ - 1)`
/// - Morris et al. viscosity
/// - CSF surface tension (optional)
/// - CFL + viscous adaptive time stepping
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WcSphSim {
    /// Particle data (SoA).
    pub particles: crate::particle::SphParticleSet,
    /// Smoothing length \[m\].
    pub h: f64,
    /// Rest density \[kg/m³\].
    pub rho0: f64,
    /// Numerical speed of sound \[m/s\].
    pub c0: f64,
    /// Tait EOS exponent (typically 7 for water).
    pub gamma: f64,
    /// Dynamic viscosity \[Pa·s\].
    pub mu: f64,
    /// Surface tension coefficient \[N/m\] (0 = disabled).
    pub sigma: f64,
    /// Gravity \[m/s²\].
    pub gravity: [f64; 3],
    /// Current simulation time \[s\].
    pub time: f64,
    /// Total kinetic energy tracker.
    pub kinetic_energy_history: Vec<f64>,
    /// Total potential energy tracker.
    pub potential_energy_history: Vec<f64>,
}
#[allow(dead_code)]
impl WcSphSim {
    /// Create a new WCSPH simulation.
    pub fn new(
        h: f64,
        rho0: f64,
        c0: f64,
        gamma: f64,
        mu: f64,
        sigma: f64,
        gravity: [f64; 3],
    ) -> Self {
        Self {
            particles: crate::particle::SphParticleSet::new(),
            h,
            rho0,
            c0,
            gamma,
            mu,
            sigma,
            gravity,
            time: 0.0,
            kinetic_energy_history: Vec::new(),
            potential_energy_history: Vec::new(),
        }
    }
    /// Tait EOS stiffness constant B = ρ₀ c₀² / γ.
    pub fn tait_b(&self) -> f64 {
        self.rho0 * self.c0 * self.c0 / self.gamma
    }
    /// **Density step**: recompute ρ via SPH kernel summation.
    pub fn density_step(&mut self) {
        self.particles.compute_sph_densities(self.h);
    }
    /// **Pressure step**: apply WCSPH Tait EOS.
    pub fn pressure_step(&mut self) {
        self.particles
            .update_wcsph_pressure(self.rho0, self.c0, self.gamma);
    }
    /// **Force step**: clear, apply gravity, pressure forces, viscosity,
    /// and optionally surface tension.
    pub fn force_step(&mut self) {
        self.particles.clear_accelerations();
        self.particles.apply_gravity(self.gravity);
        self.particles.accumulate_pressure_acceleration(self.h);
        self.particles.accumulate_morris_viscosity(self.h, self.mu);
        if self.sigma.abs() > 1e-30 {
            self.accumulate_csf_surface_tension();
        }
    }
    /// **Integrate** with Euler (v += a*dt, r += v*dt).
    pub fn integrate(&mut self, dt: f64) {
        self.particles.integrate(dt);
        self.time += dt;
    }
    /// **Full step**: density → pressure → forces → integrate.
    pub fn step(&mut self, dt: f64) {
        self.density_step();
        self.pressure_step();
        self.force_step();
        self.integrate(dt);
        let ke = self.kinetic_energy();
        let pe = self.potential_energy();
        self.kinetic_energy_history.push(ke);
        self.potential_energy_history.push(pe);
    }
    /// **Full step** using leapfrog integrator for better energy conservation.
    pub fn step_leapfrog(&mut self, dt: f64) {
        self.particles.leapfrog_kick_drift(dt);
        self.density_step();
        self.pressure_step();
        self.force_step();
        self.particles.finish_verlet(dt);
        self.time += dt;
    }
    /// **Full step** using Velocity Verlet integrator.
    pub fn step_verlet(&mut self, dt: f64) {
        self.particles.integrate_verlet_half(dt);
        self.density_step();
        self.pressure_step();
        self.force_step();
        self.particles.finish_verlet(dt);
        self.time += dt;
    }
    /// Compute the CFL-limited adaptive time step.
    pub fn cfl_dt(&self, cfl: f64, dt_min: f64, dt_max: f64) -> f64 {
        let nu = if self.rho0 > 1e-20 {
            self.mu / self.rho0
        } else {
            0.0
        };
        self.particles
            .adaptive_timestep(self.h, self.c0, nu, cfl, dt_min, dt_max)
    }
    /// Run the simulation adaptively for `total_time` seconds.
    pub fn run_adaptive(&mut self, total_time: f64, cfl: f64, dt_min: f64, dt_max: f64) {
        let t_end = self.time + total_time;
        while self.time < t_end {
            let mut dt = self.cfl_dt(cfl, dt_min, dt_max);
            dt = dt.min(t_end - self.time);
            if dt <= 0.0 {
                break;
            }
            self.step(dt);
        }
    }
    /// Total kinetic energy: Σ ½ mᵢ |vᵢ|².
    pub fn kinetic_energy(&self) -> f64 {
        self.particles.kinetic_energy()
    }
    /// Total gravitational potential energy: Σ −mᵢ (g · rᵢ).
    pub fn potential_energy(&self) -> f64 {
        self.particles.potential_energy(self.gravity)
    }
    /// Total mechanical energy.
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }
    /// Total linear momentum.
    pub fn total_momentum(&self) -> [f64; 3] {
        self.particles.total_momentum()
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.particles.len()
    }
    /// True if no particles.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }
    /// CSF surface tension force accumulation (colour-function proxy = density).
    ///
    /// Adds `F_st = σ * κ * n̂ * mᵢ / ρᵢ` to accelerations.
    fn accumulate_csf_surface_tension(&mut self) {
        let n = self.particles.len();
        let h = self.h;
        let h2 = h * h;
        let sigma = self.sigma;
        let rho0 = self.rho0;
        let mut grad_c = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            let rho_i = self.particles.densities[i].max(1e-14);
            let ci = rho_i / rho0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.particles.positions[i][0] - self.particles.positions[j][0];
                let dy = self.particles.positions[i][1] - self.particles.positions[j][1];
                let dz = self.particles.positions[i][2] - self.particles.positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 > 4.0 * h2 || r2 < 1e-28 {
                    continue;
                }
                let r = r2.sqrt();
                let rho_j = self.particles.densities[j].max(1e-14);
                let cj = rho_j / rho0;
                let dw = crate::particle::cubic_spline_gradient_pub(r, h);
                let factor = self.particles.masses[j] / rho_j * (cj - ci) * dw / r;
                grad_c[i][0] += factor * (-dx);
                grad_c[i][1] += factor * (-dy);
                grad_c[i][2] += factor * (-dz);
            }
        }
        for i in 0..n {
            if self.particles.is_boundary[i] {
                continue;
            }
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
            let rho_i = self.particles.densities[i].max(1e-14);
            let mut kappa = 0.0_f64;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.particles.positions[i][0] - self.particles.positions[j][0];
                let dy = self.particles.positions[i][1] - self.particles.positions[j][1];
                let dz = self.particles.positions[i][2] - self.particles.positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 > 4.0 * h2 || r2 < 1e-28 {
                    continue;
                }
                let r = r2.sqrt();
                let rho_j = self.particles.densities[j].max(1e-14);
                let gc_j_mag2 = grad_c[j][0] * grad_c[j][0]
                    + grad_c[j][1] * grad_c[j][1]
                    + grad_c[j][2] * grad_c[j][2];
                if gc_j_mag2 < 1e-20 {
                    continue;
                }
                let gc_j_mag = gc_j_mag2.sqrt();
                let n_hat_j = [
                    grad_c[j][0] / gc_j_mag,
                    grad_c[j][1] / gc_j_mag,
                    grad_c[j][2] / gc_j_mag,
                ];
                let dw = crate::particle::cubic_spline_gradient_pub(r, h);
                let dn_dot_r = (n_hat_j[0] - n_hat[0]) * (-dx)
                    + (n_hat_j[1] - n_hat[1]) * (-dy)
                    + (n_hat_j[2] - n_hat[2]) * (-dz);
                kappa -= self.particles.masses[j] / rho_j * dn_dot_r * dw / r;
            }
            let force_scale = sigma * kappa * gc_mag / rho_i;
            for k in 0..3 {
                self.particles.accelerations[i][k] += force_scale * n_hat[k];
            }
        }
    }
}
/// High-level SPH simulation.
pub struct SphSimulation {
    /// Fluid particles.
    pub particles: ParticleSet,
    /// Boundary elements.
    pub boundaries: BoundarySet,
    /// Simulation parameters.
    pub params: SphSimulationParams,
    /// DFSPH solver state (only used when solver_type == Dfsph).
    dfsph_solver: DfsphSolver,
    /// Current simulation time.
    pub time: f64,
}
impl SphSimulation {
    /// Create a new SPH simulation.
    pub fn new(params: SphSimulationParams) -> Self {
        Self {
            particles: ParticleSet::new(),
            boundaries: BoundarySet::new(),
            params,
            dfsph_solver: DfsphSolver::new(0),
            time: 0.0,
        }
    }
    /// Advance the simulation by one timestep `dt` using the specified kernel.
    pub fn step(&mut self, dt: f64, kernel: &dyn SphKernel) {
        let h = self.params.smoothing_length;
        let support_radius = 2.0 * h;
        let neighbors = SpatialHash::find_all_neighbors(&self.particles.positions, support_radius);
        match self.params.solver_type {
            SolverType::Wcsph => {
                let wcsph_params = WcsphParams {
                    rest_density: self.params.rest_density,
                    stiffness: self.params.sound_speed
                        * self.params.sound_speed
                        * self.params.rest_density
                        / 7.0,
                    gamma: 7.0,
                    viscosity: self.params.viscosity,
                    smoothing_length: h,
                };
                crate::wcsph::step(
                    &mut self.particles,
                    &neighbors,
                    kernel,
                    &wcsph_params,
                    dt,
                    self.params.gravity,
                );
            }
            SolverType::Dfsph => {
                let dfsph_params = DfsphParams {
                    rest_density: self.params.rest_density,
                    smoothing_length: h,
                    viscosity: self.params.viscosity,
                    ..Default::default()
                };
                self.dfsph_solver.resize(self.particles.len());
                crate::dfsph::step(
                    &mut self.particles,
                    &neighbors,
                    kernel,
                    &mut self.dfsph_solver,
                    &dfsph_params,
                    dt,
                    self.params.gravity,
                );
            }
        }
        crate::boundary_sph::apply_boundary_forces(
            &mut self.particles,
            &self.boundaries,
            kernel,
            h,
            self.params.rest_density,
            self.params.boundary_stiffness,
            self.params.boundary_damping,
        );
        self.time += dt;
    }
    /// Run the simulation for `total_time` with adaptive time stepping.
    pub fn run(&mut self, total_time: f64, kernel: &dyn SphKernel) {
        let end_time = self.time + total_time;
        while self.time < end_time {
            let dt = timestep::adaptive_timestep(
                &self.particles,
                self.params.smoothing_length,
                self.params.sound_speed,
                self.params.viscosity,
                self.params.min_dt,
                self.params.max_dt,
            );
            let dt = dt.min(end_time - self.time);
            if dt <= 0.0 {
                break;
            }
            self.step(dt, kernel);
        }
    }
}
impl SphSimulation {
    /// Wrap all particle positions into the periodic box `[0, box_size)^3`.
    ///
    /// This implements periodic boundary conditions (PBC): any coordinate
    /// that drifts outside the box is folded back using modular arithmetic.
    #[allow(dead_code)]
    pub fn apply_pbc(&mut self, box_size: [f64; 3]) {
        for p in &mut self.particles.positions {
            for k in 0..3 {
                let l = box_size[k];
                if l <= 0.0 {
                    continue;
                }
                p[k] = p[k].rem_euclid(l);
            }
        }
    }
    /// Total mechanical energy = kinetic energy + gravitational potential energy.
    ///
    /// Potential energy is computed relative to the y = 0 plane using
    /// `E_pot = sum_i m_i * |g| * y_i` where |g| is the magnitude of gravity.
    #[allow(dead_code)]
    pub fn total_energy(&self) -> f64 {
        let ke = self
            .particles
            .velocities
            .iter()
            .zip(self.particles.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v.x * v.x + v.y * v.y + v.z * v.z))
            .sum::<f64>();
        let g_mag = (self.params.gravity.x * self.params.gravity.x
            + self.params.gravity.y * self.params.gravity.y
            + self.params.gravity.z * self.params.gravity.z)
            .sqrt();
        let pe = self
            .particles
            .positions
            .iter()
            .zip(self.particles.masses.iter())
            .map(|(p, &m)| m * g_mag * p.y)
            .sum::<f64>();
        ke + pe
    }
    /// Total linear momentum vector `p = sum_i m_i * v_i`.
    ///
    /// Returns `[px, py, pz]`.
    #[allow(dead_code)]
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut mom = [0.0_f64; 3];
        for (v, &m) in self
            .particles
            .velocities
            .iter()
            .zip(self.particles.masses.iter())
        {
            mom[0] += m * v.x;
            mom[1] += m * v.y;
            mom[2] += m * v.z;
        }
        mom
    }
    /// Compute the Courant number for the current state.
    ///
    /// The Courant–Friedrichs–Lewy (CFL) condition for SPH is:
    ///   `C = (v_max + c0) * dt / h`
    ///
    /// where `v_max` is the maximum particle speed, `c0` is the speed of
    /// sound, `dt` is the current time step and `h` is the smoothing length.
    ///
    /// A value C ≤ 0.4 is typically required for stability.
    ///
    /// # Arguments
    /// * `dt` – Current time step \[s\].
    #[allow(dead_code)]
    pub fn compute_courant_number(&self, dt: f64) -> f64 {
        let v_max = self
            .particles
            .velocities
            .iter()
            .map(|v| (v.x * v.x + v.y * v.y + v.z * v.z).sqrt())
            .fold(0.0_f64, f64::max);
        let h = self.params.smoothing_length.max(1e-30);
        (v_max + self.params.sound_speed) * dt / h
    }
    /// Compute an adaptive time step satisfying the CFL condition.
    ///
    /// The CFL-stable time step is:
    ///   `dt = cfl_factor * h / (v_max + c0)`
    ///
    /// clamped to `[min_dt, max_dt]` from the simulation parameters.
    ///
    /// # Arguments
    /// * `cfl_factor` – Safety factor (commonly 0.4).
    #[allow(dead_code)]
    pub fn adaptive_timestep(&self, cfl_factor: f64) -> f64 {
        let v_max = self
            .particles
            .velocities
            .iter()
            .map(|v| (v.x * v.x + v.y * v.y + v.z * v.z).sqrt())
            .fold(0.0_f64, f64::max);
        let h = self.params.smoothing_length.max(1e-30);
        let denom = (v_max + self.params.sound_speed).max(1e-30);
        let dt_cfl = cfl_factor * h / denom;
        dt_cfl.clamp(self.params.min_dt, self.params.max_dt)
    }
    /// Compute total mechanical energy: kinetic + gravitational potential.
    ///
    /// Kinetic energy: `KE = Σ 0.5 * m_i * |v_i|²`
    /// Potential energy: `PE = Σ m_i * |g| * y_i`
    ///
    /// This is an alias exposed on `SphSimulation` for the educational API;
    /// see also [`total_energy`](Self::total_energy).
    #[allow(dead_code)]
    pub fn compute_mechanical_energy(&self) -> f64 {
        self.total_energy()
    }
}
impl SphSimulation {
    /// Apply periodic boundary conditions to all fluid particles.
    ///
    /// Wraps each coordinate into `[0, L)` for the given box sizes `box_size`.
    #[allow(dead_code)]
    pub fn apply_periodic_bc(&mut self, box_size: [f64; 3]) {
        for pos in &mut self.particles.positions {
            for k in 0..3 {
                let l = box_size[k];
                if l <= 0.0 {
                    continue;
                }
                while pos[k] < 0.0 {
                    pos[k] += l;
                }
                while pos[k] >= l {
                    pos[k] -= l;
                }
            }
        }
    }
    /// Compute mean density across all fluid particles.
    #[allow(dead_code)]
    pub fn mean_density(&self) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        self.particles.densities.iter().sum::<f64>() / n as f64
    }
    /// Compute density standard deviation across all particles.
    #[allow(dead_code)]
    pub fn density_std(&self) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        let mean = self.mean_density();
        let var = self
            .particles
            .densities
            .iter()
            .map(|&d| (d - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        var.sqrt()
    }
    /// Maximum speed among all fluid particles.
    #[allow(dead_code)]
    pub fn max_speed_simulation(&self) -> f64 {
        self.particles
            .positions
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let v = &self.particles.velocities[i];
                (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
            })
            .fold(0.0_f64, f64::max)
    }
    /// Number of fluid particles that have exceeded a speed threshold `v_max`.
    #[allow(dead_code)]
    pub fn count_fast_particles(&self, v_max: f64) -> usize {
        let v2_max = v_max * v_max;
        self.particles
            .velocities
            .iter()
            .filter(|v| v.x * v.x + v.y * v.y + v.z * v.z > v2_max)
            .count()
    }
    /// Rescale all particle velocities by `factor` (e.g. for temperature
    /// rescaling or damping).
    #[allow(dead_code)]
    pub fn rescale_velocities(&mut self, factor: f64) {
        for v in &mut self.particles.velocities {
            v.x *= factor;
            v.y *= factor;
            v.z *= factor;
        }
    }
    /// Sum of all particle masses.
    #[allow(dead_code)]
    pub fn total_mass(&self) -> f64 {
        self.particles.masses.iter().sum()
    }
    /// Centre of mass of all particles.
    #[allow(dead_code)]
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_m = self.total_mass();
        if total_m == 0.0 {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for (pos, &m) in self
            .particles
            .positions
            .iter()
            .zip(self.particles.masses.iter())
        {
            com[0] += m * pos.x;
            com[1] += m * pos.y;
            com[2] += m * pos.z;
        }
        com[0] /= total_m;
        com[1] /= total_m;
        com[2] /= total_m;
        com
    }
    /// Total kinetic energy Σ ½ mᵢ |vᵢ|² of the simulation.
    #[allow(dead_code)]
    pub fn kinetic_energy_sim(&self) -> f64 {
        self.particles
            .velocities
            .iter()
            .zip(self.particles.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v.x * v.x + v.y * v.y + v.z * v.z))
            .sum()
    }
    /// Add a single fluid particle at `pos` with velocity `vel` and `mass`.
    #[allow(dead_code)]
    pub fn add_particle_raw(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64) {
        use crate::particle::SphParticle;
        use oxiphysics_core::math::Vec3;
        let p = SphParticle::new(
            Vec3::new(pos[0], pos[1], pos[2]),
            Vec3::new(vel[0], vel[1], vel[2]),
            mass,
        );
        self.particles.add_particle(&p);
    }
    /// Advance the simulation by exactly `n` steps with a fixed time step `dt`
    /// and the given `kernel`.
    #[allow(dead_code)]
    pub fn step_n(&mut self, n: usize, dt: f64, kernel: &dyn crate::kernel::SphKernel) {
        for _ in 0..n {
            self.step(dt, kernel);
        }
    }
    /// Bounding box of all particle positions: `(lo, hi)`.
    #[allow(dead_code)]
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.particles.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for pos in &self.particles.positions {
            if pos.x < lo[0] {
                lo[0] = pos.x;
            }
            if pos.y < lo[1] {
                lo[1] = pos.y;
            }
            if pos.z < lo[2] {
                lo[2] = pos.z;
            }
            if pos.x > hi[0] {
                hi[0] = pos.x;
            }
            if pos.y > hi[1] {
                hi[1] = pos.y;
            }
            if pos.z > hi[2] {
                hi[2] = pos.z;
            }
        }
        (lo, hi)
    }
}
/// Simulation state tracking two immiscible SPH fluid phases.
///
/// Phase A is stored in `particles_a` and phase B in `particles_b`; a
/// combined step advances both according to their respective configurations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TwoPhaseSimState {
    /// Phase-A fluid particles.
    pub particles_a: SphSim,
    /// Phase-B fluid particles.
    pub particles_b: SphSim,
    /// Interfacial surface-tension coefficient (N/m).
    pub sigma_ab: f64,
    /// Current time.
    pub time: f64,
}
#[allow(dead_code)]
impl TwoPhaseSimState {
    /// Construct from two existing [`SphSim`] instances.
    pub fn new(a: SphSim, b: SphSim, sigma_ab: f64) -> Self {
        Self {
            particles_a: a,
            particles_b: b,
            sigma_ab,
            time: 0.0,
        }
    }
    /// Advance both phases by one step.
    pub fn step(&mut self) {
        let dt = self.particles_a.config.dt;
        self.particles_a.step();
        self.particles_b.step();
        self.time += dt;
    }
    /// Total number of particles across both phases.
    pub fn total_particles(&self) -> usize {
        self.particles_a.len() + self.particles_b.len()
    }
    /// Total kinetic energy of both phases combined.
    pub fn total_ke(&self) -> f64 {
        let ke_a: f64 = self
            .particles_a
            .velocities
            .iter()
            .zip(self.particles_a.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum();
        let ke_b: f64 = self
            .particles_b
            .velocities
            .iter()
            .zip(self.particles_b.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum();
        ke_a + ke_b
    }
    /// Total momentum of both phases combined.
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut mom = [0.0_f64; 3];
        for (v, &m) in self
            .particles_a
            .velocities
            .iter()
            .zip(self.particles_a.masses.iter())
            .chain(
                self.particles_b
                    .velocities
                    .iter()
                    .zip(self.particles_b.masses.iter()),
            )
        {
            mom[0] += m * v[0];
            mom[1] += m * v[1];
            mom[2] += m * v[2];
        }
        mom
    }
}
/// A minimal periodic-box SPH simulation using plain arrays.
///
/// Positions wrap with PBC; provides `total_energy()` and `total_momentum()`
/// for conservation diagnostics. Inherits the SPH pipeline from [`SphSim`].
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PeriodicSphSim {
    /// Inner simulation state.
    pub inner: SphSim,
    /// Simulation box lengths `[Lx, Ly, Lz]` \[m\].
    pub box_size: [f64; 3],
}
#[allow(dead_code)]
impl PeriodicSphSim {
    /// Create a periodic simulation wrapping an existing [`SphSim`].
    pub fn new(inner: SphSim, box_size: [f64; 3]) -> Self {
        Self { inner, box_size }
    }
    /// Advance one step with periodic boundary conditions.
    ///
    /// Runs `SphSim::step()` then wraps positions into the periodic box.
    pub fn step(&mut self) {
        self.inner.step();
        self.apply_pbc();
    }
    /// Advance one step including surface tension and PBC.
    pub fn step_with_surface_tension(&mut self) {
        self.inner.step_with_surface_tension();
        self.apply_pbc();
    }
    /// Wrap all positions into `[0, box_size)^3`.
    pub fn apply_pbc(&mut self) {
        for pos in &mut self.inner.positions {
            for k in 0..3 {
                let l = self.box_size[k];
                if l > 0.0 {
                    pos[k] = pos[k].rem_euclid(l);
                }
            }
        }
    }
    /// Total mechanical energy = KE + gravitational PE.
    ///
    /// PE is computed as `m * |g| * y` relative to the bottom of the box.
    pub fn total_energy(&self) -> f64 {
        let sim = &self.inner;
        let ke: f64 = sim
            .velocities
            .iter()
            .zip(sim.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum();
        let g = sim.config.gravity;
        let g_mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        let pe: f64 = sim
            .positions
            .iter()
            .zip(sim.masses.iter())
            .map(|(p, &m)| m * g_mag * p[1])
            .sum();
        ke + pe
    }
    /// Total linear momentum `[px, py, pz]`.
    pub fn total_momentum(&self) -> [f64; 3] {
        let sim = &self.inner;
        let mut mom = [0.0_f64; 3];
        for (v, &m) in sim.velocities.iter().zip(sim.masses.iter()) {
            mom[0] += m * v[0];
            mom[1] += m * v[1];
            mom[2] += m * v[2];
        }
        mom
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    /// Whether empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}
/// Standalone CFL time step controller for SPH simulations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CflController {
    /// CFL number (Courant number).  Typical value: 0.4.
    pub cfl: f64,
    /// Viscous safety factor.  Typical value: 0.125.
    pub viscous_factor: f64,
    /// Minimum allowed timestep.
    pub dt_min: f64,
    /// Maximum allowed timestep.
    pub dt_max: f64,
}
#[allow(dead_code)]
impl CflController {
    /// Create a controller with standard SPH defaults.
    pub fn standard(dt_min: f64, dt_max: f64) -> Self {
        Self {
            cfl: 0.4,
            viscous_factor: 0.125,
            dt_min,
            dt_max,
        }
    }
    /// CFL time step: `dt = cfl * h / (c0 + v_max)`.
    pub fn cfl_dt(&self, h: f64, c0: f64, v_max: f64) -> f64 {
        let denom = c0 + v_max;
        if denom < 1e-30 {
            return self.dt_max;
        }
        self.cfl * h / denom
    }
    /// Viscous time step: `dt = viscous_factor * h² / ν`.
    pub fn viscous_dt(&self, h: f64, nu: f64) -> f64 {
        if nu < 1e-30 {
            return self.dt_max;
        }
        self.viscous_factor * h * h / nu
    }
    /// Combined time step: `min(cfl_dt, viscous_dt)` clamped to `[dt_min, dt_max]`.
    pub fn combined_dt(&self, h: f64, c0: f64, v_max: f64, nu: f64) -> f64 {
        let dt = self.cfl_dt(h, c0, v_max).min(self.viscous_dt(h, nu));
        dt.clamp(self.dt_min, self.dt_max)
    }
}
