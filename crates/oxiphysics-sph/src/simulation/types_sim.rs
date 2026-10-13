//! Simulation types: WcSphSim, SphSimulation, TwoPhaseSimState, PeriodicSphSim, CflController.

use crate::boundary_sph::BoundarySet;
use crate::dfsph::{DfsphParams, DfsphSolver};
use crate::kernel::SphKernel;
use crate::neighbor::SpatialHash;
use crate::particle::ParticleSet;
use crate::timestep;
use crate::wcsph::WcsphParams;

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
    /// Neighbour lists (indices `j ≠ i` with `|xᵢ − xⱼ| ≤ 2h`, boundary
    /// particles included) for the cubic-spline support of this simulation.
    ///
    /// One list per step serves the density, pressure, viscosity and surface
    /// tension sums; see [`density_step_with_neighbors`](Self::density_step_with_neighbors)
    /// and [`force_step_with_neighbors`](Self::force_step_with_neighbors).
    pub fn neighbor_lists(&self) -> Vec<Vec<usize>> {
        self.particles.support_neighbor_lists(self.h)
    }
    /// **Density step**: recompute ρ via SPH kernel summation.
    ///
    /// Builds its own neighbour lists; [`step`](Self::step) shares one list
    /// between the density and force sums instead.
    pub fn density_step(&mut self) {
        self.particles.compute_sph_densities(self.h);
    }
    /// **Density step** from neighbour lists built by
    /// [`neighbor_lists`](Self::neighbor_lists) for the current positions.
    pub fn density_step_with_neighbors(&mut self, neighbors: &[Vec<usize>]) {
        self.particles
            .compute_sph_densities_with_neighbors(self.h, neighbors);
    }
    /// **Pressure step**: apply WCSPH Tait EOS.
    pub fn pressure_step(&mut self) {
        self.particles
            .update_wcsph_pressure(self.rho0, self.c0, self.gamma);
    }
    /// **Force step**: clear, apply gravity, pressure forces, viscosity,
    /// and optionally surface tension.
    ///
    /// Builds its own neighbour lists; see
    /// [`force_step_with_neighbors`](Self::force_step_with_neighbors).
    pub fn force_step(&mut self) {
        let neighbors = self.neighbor_lists();
        self.force_step_with_neighbors(&neighbors);
    }
    /// **Force step** from neighbour lists built by
    /// [`neighbor_lists`](Self::neighbor_lists) for the current positions:
    /// clear, apply gravity, pressure forces, viscosity and (when `sigma ≠ 0`)
    /// the continuum surface force.  Every sum costs `O(N · k)` for `k`
    /// neighbours per particle.
    pub fn force_step_with_neighbors(&mut self, neighbors: &[Vec<usize>]) {
        self.particles.clear_accelerations();
        self.particles.apply_gravity(self.gravity);
        self.particles
            .accumulate_pressure_acceleration_with_neighbors(self.h, neighbors);
        self.particles
            .accumulate_morris_viscosity_with_neighbors(self.h, self.mu, neighbors);
        if self.sigma.abs() > 1e-30 {
            self.accumulate_csf_surface_tension(neighbors);
        }
    }
    /// **Integrate** with Euler (v += a*dt, r += v*dt).
    pub fn integrate(&mut self, dt: f64) {
        self.particles.integrate(dt);
        self.time += dt;
    }
    /// One density → pressure → force evaluation at the current positions with
    /// a single set of neighbour lists.
    fn evaluate_forces(&mut self) {
        let neighbors = self.neighbor_lists();
        self.density_step_with_neighbors(&neighbors);
        self.pressure_step();
        self.force_step_with_neighbors(&neighbors);
    }
    /// **Full step**: density → pressure → forces → integrate.
    ///
    /// The neighbour lists are built once and shared by all pair sums.
    pub fn step(&mut self, dt: f64) {
        self.evaluate_forces();
        self.integrate(dt);
        let ke = self.kinetic_energy();
        let pe = self.potential_energy();
        self.kinetic_energy_history.push(ke);
        self.potential_energy_history.push(pe);
    }
    /// **Full step** using leapfrog integrator for better energy conservation.
    pub fn step_leapfrog(&mut self, dt: f64) {
        self.particles.leapfrog_kick_drift(dt);
        self.evaluate_forces();
        self.particles.finish_verlet(dt);
        self.time += dt;
    }
    /// **Full step** using Velocity Verlet integrator.
    pub fn step_verlet(&mut self, dt: f64) {
        self.particles.integrate_verlet_half(dt);
        self.evaluate_forces();
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
    /// Brackbill-Kothe-Zemach continuum surface force with SPH operators over
    /// the neighbour lists (cubic-spline support `2h`):
    ///
    /// 1. colour `cᵢ = ρᵢ/ρ₀` and its gradient
    ///    `∇cᵢ = Σⱼ (mⱼ/ρⱼ)(cⱼ − cᵢ) ∇ᵢWᵢⱼ` (points towards larger `c`, i.e.
    ///    into the fluid);
    /// 2. unit normal `n̂ᵢ = ∇cᵢ/|∇cᵢ|` wherever `|∇cᵢ|` is non-negligible;
    /// 3. curvature `κᵢ = −∇·n̂ᵢ = −Σⱼ (mⱼ/ρⱼ)(n̂ⱼ − n̂ᵢ)·∇ᵢWᵢⱼ`, positive for a
    ///    convex fluid region (`2/R` on a sphere of radius `R`);
    /// 4. acceleration `aᵢ += σ κᵢ ∇cᵢ / ρᵢ` for every fluid particle, which
    ///    pulls convex surfaces inwards.
    ///
    /// Boundary particles take part in the sums but receive no acceleration.
    fn accumulate_csf_surface_tension(&mut self, neighbors: &[Vec<usize>]) {
        const NORMAL_EPS2: f64 = 1e-20;
        let n = self.particles.len();
        let h = self.h;
        let h2 = h * h;
        let sigma = self.sigma;
        let rho0 = self.rho0;
        let positions = &self.particles.positions;
        let densities = &self.particles.densities;
        let masses = &self.particles.masses;
        let colour: Vec<f64> = densities.iter().map(|&rho| rho.max(1e-14) / rho0).collect();
        // Geometry of the pair (i, j) inside the support: the separation
        // d = x_i − x_j and the weight V_j (dW/dr) / r, so that
        // ∇_i W_ij = (dW/dr) d / r and V_j ∇_i W_ij = weight · d.
        let pair_weight = |i: usize, j: usize| -> Option<([f64; 3], f64)> {
            if j == i || j >= n {
                return None;
            }
            let d = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if r2 > 4.0 * h2 || r2 < 1e-28 {
                return None;
            }
            let r = r2.sqrt();
            let rho_j = densities[j].max(1e-14);
            let dw = crate::particle::cubic_spline_gradient_pub(r, h);
            Some((d, masses[j] / rho_j * dw / r))
        };
        // ∇c_i = Σ_j V_j (c_j − c_i) ∇_i W_ij
        let mut grad_c = vec![[0.0_f64; 3]; n];
        for (i, g) in grad_c.iter_mut().enumerate() {
            let Some(list) = neighbors.get(i) else {
                continue;
            };
            for &j in list {
                if let Some((d, weight)) = pair_weight(i, j) {
                    let factor = weight * (colour[j] - colour[i]);
                    g[0] += factor * d[0];
                    g[1] += factor * d[1];
                    g[2] += factor * d[2];
                }
            }
        }
        let normals: Vec<Option<[f64; 3]>> = grad_c
            .iter()
            .map(|g| {
                let mag2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
                if mag2 < NORMAL_EPS2 {
                    None
                } else {
                    let mag = mag2.sqrt();
                    Some([g[0] / mag, g[1] / mag, g[2] / mag])
                }
            })
            .collect();
        for (i, list) in neighbors.iter().enumerate().take(n) {
            if self.particles.is_boundary[i] {
                continue;
            }
            let Some(n_i) = normals[i] else {
                continue;
            };
            // κ_i = −∇·n̂_i = −Σ_j V_j (n̂_j − n̂_i) · ∇_i W_ij
            let mut kappa = 0.0_f64;
            for &j in list {
                let (Some((d, weight)), Some(n_j)) =
                    (pair_weight(i, j), normals.get(j).copied().flatten())
                else {
                    continue;
                };
                let dn_dot_d =
                    (n_j[0] - n_i[0]) * d[0] + (n_j[1] - n_i[1]) * d[1] + (n_j[2] - n_i[2]) * d[2];
                kappa -= weight * dn_dot_d;
            }
            let rho_i = densities[i].max(1e-14);
            let scale = sigma * kappa / rho_i;
            for (k, &gk) in grad_c[i].iter().enumerate() {
                self.particles.accelerations[i][k] += scale * gk;
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
    pub fn compute_mechanical_energy(&self) -> f64 {
        self.total_energy()
    }
}
impl SphSimulation {
    /// Apply periodic boundary conditions to all fluid particles.
    ///
    /// Wraps each coordinate into `[0, L)` for the given box sizes `box_size`.
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
    pub fn mean_density(&self) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        self.particles.densities.iter().sum::<f64>() / n as f64
    }
    /// Compute density standard deviation across all particles.
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
    pub fn rescale_velocities(&mut self, factor: f64) {
        for v in &mut self.particles.velocities {
            v.x *= factor;
            v.y *= factor;
            v.z *= factor;
        }
    }
    /// Sum of all particle masses.
    pub fn total_mass(&self) -> f64 {
        self.particles.masses.iter().sum()
    }
    /// Centre of mass of all particles.
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
    pub fn kinetic_energy_sim(&self) -> f64 {
        self.particles
            .velocities
            .iter()
            .zip(self.particles.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v.x * v.x + v.y * v.y + v.z * v.z))
            .sum()
    }
    /// Add a single fluid particle at `pos` with velocity `vel` and `mass`.
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
    pub fn step_n(&mut self, n: usize, dt: f64, kernel: &dyn crate::kernel::SphKernel) {
        for _ in 0..n {
            self.step(dt, kernel);
        }
    }
    /// Bounding box of all particle positions: `(lo, hi)`.
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
pub struct PeriodicSphSim {
    /// Inner simulation state.
    pub inner: SphSim,
    /// Simulation box lengths `[Lx, Ly, Lz]` \[m\].
    pub box_size: [f64; 3],
}
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
            for (k, p) in pos.iter_mut().enumerate().take(3) {
                let l = self.box_size[k];
                if l > 0.0 {
                    *p = p.rem_euclid(l);
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
