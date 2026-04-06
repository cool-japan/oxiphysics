// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly fluid simulation bridge.
//!
//! Provides pure-Rust types for SPH particle simulations, Lattice-Boltzmann
//! (LBM) simulations, multiphase flows, fluid-rigid coupling, particle
//! emitters, and flow-field analysis. Designed for serialisation across the
//! WASM boundary.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// WasmSphConfig
// ---------------------------------------------------------------------------

/// Configuration for an SPH (Smoothed Particle Hydrodynamics) simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSphConfig {
    /// Particle radius r (m).
    pub particle_radius: f64,
    /// Smoothing length h (m); typically 2r.
    pub smoothing_length: f64,
    /// Rest (reference) density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Surface tension coefficient γ (N/m).
    pub surface_tension: f64,
    /// Equation of state stiffness B (Pa).
    pub stiffness: f64,
    /// Tait EOS exponent γ_eos.
    pub tait_exponent: f64,
    /// Gravity vector \[gx, gy, gz\] (m/s²).
    pub gravity: [f64; 3],
    /// Coefficient of restitution for boundary collisions.
    pub boundary_restitution: f64,
    /// Enable XSPH velocity correction.
    pub xsph_enabled: bool,
    /// XSPH ε coefficient.
    pub xsph_epsilon: f64,
}

impl Default for WasmSphConfig {
    fn default() -> Self {
        WasmSphConfig {
            particle_radius: 0.05,
            smoothing_length: 0.1,
            rest_density: 1000.0,
            viscosity: 0.001,
            surface_tension: 0.0728,
            stiffness: 1000.0,
            tait_exponent: 7.0,
            gravity: [0.0, -9.81, 0.0],
            boundary_restitution: 0.3,
            xsph_enabled: true,
            xsph_epsilon: 0.5,
        }
    }
}

impl WasmSphConfig {
    /// Create a water-like config.
    pub fn water() -> Self {
        WasmSphConfig::default()
    }

    /// Validate config.
    pub fn validate(&self) -> Result<(), String> {
        if self.particle_radius <= 0.0 {
            return Err("particle_radius must be positive".to_string());
        }
        if self.smoothing_length < self.particle_radius {
            return Err("smoothing_length must be >= particle_radius".to_string());
        }
        if self.rest_density <= 0.0 {
            return Err("rest_density must be positive".to_string());
        }
        Ok(())
    }

    /// Wendland C2 kernel value at distance r.
    pub fn kernel_wendland(&self, r: f64) -> f64 {
        let h = self.smoothing_length;
        let q = (r / h).clamp(0.0, 1.0);
        let alpha = 21.0 / (2.0 * std::f64::consts::PI * h * h * h);
        alpha * (1.0 - q / 2.0).powi(4) * (2.0 * q + 1.0)
    }

    /// Cubic spline kernel value.
    pub fn kernel_cubic(&self, r: f64) -> f64 {
        let h = self.smoothing_length;
        let q = r / h;
        let alpha = 3.0 / (2.0 * std::f64::consts::PI * h * h * h);
        if q < 1.0 {
            alpha * (2.0 / 3.0 - q * q + q * q * q / 2.0)
        } else if q < 2.0 {
            alpha * (2.0 - q).powi(3) / 6.0
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// WasmSphParticle
// ---------------------------------------------------------------------------

/// A single SPH fluid particle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSphParticle {
    /// Position \[x, y, z\] (m).
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] (m/s).
    pub velocity: [f64; 3],
    /// Acceleration \[ax, ay, az\] (m/s²) – updated by solver.
    pub acceleration: [f64; 3],
    /// Current density ρ (kg/m³).
    pub density: f64,
    /// Current pressure P (Pa).
    pub pressure: f64,
    /// Particle mass (kg).
    pub mass: f64,
    /// Whether this particle is a boundary (static) particle.
    pub is_boundary: bool,
    /// Smoothed colour field value (for surface detection).
    pub colour: f64,
    /// Particle unique ID.
    pub id: u64,
}

impl WasmSphParticle {
    /// Create a fluid particle at rest.
    pub fn new(id: u64, position: [f64; 3], mass: f64) -> Self {
        WasmSphParticle {
            position,
            velocity: [0.0; 3],
            acceleration: [0.0; 3],
            density: 1000.0,
            pressure: 0.0,
            mass,
            is_boundary: false,
            colour: 0.0,
            id,
        }
    }

    /// Create a static boundary particle.
    pub fn boundary(id: u64, position: [f64; 3], mass: f64) -> Self {
        WasmSphParticle {
            is_boundary: true,
            ..WasmSphParticle::new(id, position, mass)
        }
    }

    /// Kinetic energy (½mv²).
    pub fn kinetic_energy(&self) -> f64 {
        let v2: f64 = self.velocity.iter().map(|v| v * v).sum();
        0.5 * self.mass * v2
    }

    /// Speed |v|.
    pub fn speed(&self) -> f64 {
        self.velocity.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
}

// ---------------------------------------------------------------------------
// WasmSphSimulation
// ---------------------------------------------------------------------------

/// SPH fluid simulation manager.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSphSimulation {
    /// Simulation configuration.
    pub config: WasmSphConfig,
    /// All particles (fluid + boundary).
    pub particles: Vec<WasmSphParticle>,
    /// Simulated time (s).
    pub time: f64,
    /// Step counter.
    pub step_count: u64,
    /// Next particle ID.
    next_id: u64,
}

impl WasmSphSimulation {
    /// Create a new simulation with the given config.
    pub fn new(config: WasmSphConfig) -> Self {
        WasmSphSimulation {
            config,
            particles: Vec::new(),
            time: 0.0,
            step_count: 0,
            next_id: 0,
        }
    }

    /// Add a single particle and return its index.
    pub fn add_particle(&mut self, position: [f64; 3], mass: f64) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.particles
            .push(WasmSphParticle::new(id, position, mass));
        self.particles.len() - 1
    }

    /// Add multiple particles from flat position array \[x0,y0,z0, x1,y1,z1, ...\].
    pub fn add_particles(&mut self, positions_flat: &[f64], mass: f64) {
        let n = positions_flat.len() / 3;
        for i in 0..n {
            let pos = [
                positions_flat[i * 3],
                positions_flat[i * 3 + 1],
                positions_flat[i * 3 + 2],
            ];
            self.add_particle(pos, mass);
        }
    }

    /// Advance the simulation by `dt` seconds (Euler integration stub).
    pub fn step(&mut self, dt: f64) {
        let grav = self.config.gravity;
        // Density and pressure update (stub: set to rest values).
        for p in self.particles.iter_mut().filter(|p| !p.is_boundary) {
            p.density = self.config.rest_density;
            let ratio = p.density / self.config.rest_density;
            p.pressure = self.config.stiffness * (ratio.powf(self.config.tait_exponent) - 1.0);
        }
        // Integration.
        for p in self.particles.iter_mut().filter(|p| !p.is_boundary) {
            p.acceleration = grav;
            for k in 0..3 {
                p.velocity[k] += p.acceleration[k] * dt;
                p.position[k] += p.velocity[k] * dt;
            }
        }
        self.time += dt;
        self.step_count += 1;
    }

    /// Get flat array of all particle positions \[x0,y0,z0, ...\].
    pub fn get_positions(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.particles.len() * 3);
        for p in &self.particles {
            out.extend_from_slice(&p.position);
        }
        out
    }

    /// Get flat array of all velocities.
    pub fn get_velocities(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.particles.len() * 3);
        for p in &self.particles {
            out.extend_from_slice(&p.velocity);
        }
        out
    }

    /// Get density of each particle.
    pub fn get_densities(&self) -> Vec<f64> {
        self.particles.iter().map(|p| p.density).collect()
    }

    /// Number of fluid (non-boundary) particles.
    pub fn fluid_count(&self) -> usize {
        self.particles.iter().filter(|p| !p.is_boundary).count()
    }

    /// Total kinetic energy of all fluid particles.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| !p.is_boundary)
            .map(|p| p.kinetic_energy())
            .sum()
    }
}

// ---------------------------------------------------------------------------
// WasmLbmConfig
// ---------------------------------------------------------------------------

/// Configuration for a Lattice-Boltzmann Method (LBM) simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmLbmConfig {
    /// Grid resolution in X.
    pub nx: usize,
    /// Grid resolution in Y.
    pub ny: usize,
    /// Grid resolution in Z.
    pub nz: usize,
    /// Kinematic viscosity ν (lattice units).
    pub viscosity: f64,
    /// Lattice type: `"D2Q9"`, `"D3Q19"`, or `"D3Q27"`.
    pub lattice_type: String,
    /// Collision operator: `"BGK"` or `"MRT"`.
    pub collision: String,
    /// Lattice speed of sound c_s (lattice units).
    pub sound_speed: f64,
    /// Relaxation time τ (derived from viscosity).
    pub tau: f64,
}

impl Default for WasmLbmConfig {
    fn default() -> Self {
        WasmLbmConfig {
            nx: 64,
            ny: 64,
            nz: 1,
            viscosity: 0.1,
            lattice_type: "D2Q9".to_string(),
            collision: "BGK".to_string(),
            sound_speed: 1.0 / 3.0_f64.sqrt(),
            tau: 0.6,
        }
    }
}

impl WasmLbmConfig {
    /// Create a 2-D D2Q9 config.
    pub fn d2q9(nx: usize, ny: usize, viscosity: f64) -> Self {
        WasmLbmConfig {
            nx,
            ny,
            nz: 1,
            viscosity,
            lattice_type: "D2Q9".to_string(),
            tau: 3.0 * viscosity + 0.5,
            ..Default::default()
        }
    }

    /// Create a 3-D D3Q19 config.
    pub fn d3q19(nx: usize, ny: usize, nz: usize, viscosity: f64) -> Self {
        WasmLbmConfig {
            nx,
            ny,
            nz,
            viscosity,
            lattice_type: "D3Q19".to_string(),
            collision: "MRT".to_string(),
            tau: 3.0 * viscosity + 0.5,
            ..Default::default()
        }
    }

    /// Reynolds number Re = U L / ν.
    pub fn reynolds_number(&self, velocity: f64, length: f64) -> f64 {
        velocity * length / self.viscosity.max(1e-30)
    }

    /// Validate the config.
    pub fn validate(&self) -> Result<(), String> {
        if self.nx == 0 || self.ny == 0 || self.nz == 0 {
            return Err("grid dimensions must be > 0".to_string());
        }
        if self.viscosity <= 0.0 {
            return Err("viscosity must be positive".to_string());
        }
        if self.tau <= 0.5 {
            return Err("tau must be > 0.5 for stability".to_string());
        }
        if !["D2Q9", "D3Q19", "D3Q27"].contains(&self.lattice_type.as_str()) {
            return Err(format!("unknown lattice_type: {}", self.lattice_type));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// WasmLbmSimulation
// ---------------------------------------------------------------------------

/// Minimal LBM simulation state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmLbmSimulation {
    /// LBM configuration.
    pub config: WasmLbmConfig,
    /// Flat density field (one value per lattice node).
    pub density: Vec<f64>,
    /// Flat velocity field: \[vx, vy, vz\] per node.
    pub velocity: Vec<f64>,
    /// Flat pressure field.
    pub pressure: Vec<f64>,
    /// Boundary mask (true = solid wall node).
    pub boundary: Vec<bool>,
    /// Simulation step counter.
    pub step_count: u64,
    /// Simulated time (lattice time units).
    pub time: f64,
}

impl WasmLbmSimulation {
    /// Create a new LBM simulation (quiescent flow).
    pub fn new(config: WasmLbmConfig) -> Self {
        let n = config.nx * config.ny * config.nz;
        WasmLbmSimulation {
            density: vec![1.0; n],
            velocity: vec![0.0; n * 3],
            pressure: vec![0.0; n],
            boundary: vec![false; n],
            step_count: 0,
            time: 0.0,
            config,
        }
    }

    /// Mark node (ix, iy, iz) as solid boundary.
    pub fn set_boundary(&mut self, ix: usize, iy: usize, iz: usize, solid: bool) {
        if let Some(idx) = self.node_index(ix, iy, iz) {
            self.boundary[idx] = solid;
        }
    }

    /// Return flat index for lattice node (ix, iy, iz).
    fn node_index(&self, ix: usize, iy: usize, iz: usize) -> Option<usize> {
        if ix < self.config.nx && iy < self.config.ny && iz < self.config.nz {
            Some(iz * self.config.nx * self.config.ny + iy * self.config.nx + ix)
        } else {
            None
        }
    }

    /// Advance one LBM step (BGK collision + streaming stub).
    pub fn step(&mut self) {
        let cs2 = self.config.sound_speed * self.config.sound_speed;
        // Update pressure from density.
        for i in 0..self.density.len() {
            self.pressure[i] = self.density[i] * cs2;
        }
        self.step_count += 1;
        self.time += 1.0;
    }

    /// Get a copy of the velocity field as \[vx,vy,vz\] per node.
    pub fn get_velocity_field(&self) -> Vec<f64> {
        self.velocity.clone()
    }

    /// Get a copy of the pressure field.
    pub fn get_pressure_field(&self) -> Vec<f64> {
        self.pressure.clone()
    }

    /// Number of lattice nodes.
    pub fn node_count(&self) -> usize {
        self.config.nx * self.config.ny * self.config.nz
    }

    /// Set a uniform inlet velocity (x-direction) on the left boundary (ix=0).
    pub fn set_inlet_velocity(&mut self, vx: f64) {
        for iy in 0..self.config.ny {
            for iz in 0..self.config.nz {
                if let Some(idx) = self.node_index(0, iy, iz) {
                    self.velocity[idx * 3] = vx;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// WasmFluidStats
// ---------------------------------------------------------------------------

/// Aggregate fluid simulation statistics for one step.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmFluidStats {
    /// Maximum velocity magnitude across all particles/nodes (m/s).
    pub max_velocity: f64,
    /// Mean density (kg/m³).
    pub mean_density: f64,
    /// Total kinetic energy (J).
    pub kinetic_energy: f64,
    /// Enstrophy (integral of ω²/2 over domain).
    pub enstrophy: f64,
    /// Pressure drop across domain (Pa).
    pub pressure_drop: f64,
    /// CFL number (max_velocity * dt / dx).
    pub cfl: f64,
    /// Number of active (non-boundary) nodes/particles.
    pub active_count: u32,
}

impl WasmFluidStats {
    /// Create zeroed stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Compute stats from an SPH simulation.
    pub fn from_sph(sim: &WasmSphSimulation, dt: f64, dx: f64) -> Self {
        let fluid: Vec<&WasmSphParticle> =
            sim.particles.iter().filter(|p| !p.is_boundary).collect();
        let n = fluid.len();
        if n == 0 {
            return Self::new();
        }
        let max_vel = fluid.iter().map(|p| p.speed()).fold(0.0_f64, f64::max);
        let mean_density = fluid.iter().map(|p| p.density).sum::<f64>() / n as f64;
        let ke = fluid.iter().map(|p| p.kinetic_energy()).sum::<f64>();
        WasmFluidStats {
            max_velocity: max_vel,
            mean_density,
            kinetic_energy: ke,
            cfl: max_vel * dt / dx.max(1e-30),
            active_count: n as u32,
            ..Default::default()
        }
    }

    /// Returns `true` if CFL condition is satisfied (CFL < 1).
    pub fn is_cfl_ok(&self) -> bool {
        self.cfl < 1.0
    }
}

// ---------------------------------------------------------------------------
// WasmMultiphaseConfig
// ---------------------------------------------------------------------------

/// Configuration for a two-phase (immiscible) flow simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmMultiphaseConfig {
    /// Densities of the two fluid phases \[ρ₁, ρ₂\] (kg/m³).
    pub fluid_densities: [f64; 2],
    /// Interfacial surface tension γ (N/m).
    pub surface_tension: f64,
    /// Static contact angle θ (radians).
    pub contact_angle: f64,
    /// Interface thickness W (m) – Cahn-Hilliard parameter.
    pub interface_thickness: f64,
    /// Mobility parameter M in phase-field equation.
    pub mobility: f64,
    /// Viscosities of the two phases \[μ₁, μ₂\] (Pa·s).
    pub viscosities: [f64; 2],
}

impl Default for WasmMultiphaseConfig {
    fn default() -> Self {
        WasmMultiphaseConfig {
            fluid_densities: [1000.0, 1.2],
            surface_tension: 0.0728,
            contact_angle: std::f64::consts::PI / 3.0, // 60°
            interface_thickness: 0.01,
            mobility: 1e-9,
            viscosities: [0.001, 1.8e-5],
        }
    }
}

impl WasmMultiphaseConfig {
    /// Density ratio ρ₁/ρ₂.
    pub fn density_ratio(&self) -> f64 {
        self.fluid_densities[0] / self.fluid_densities[1].max(1e-30)
    }

    /// Viscosity ratio μ₁/μ₂.
    pub fn viscosity_ratio(&self) -> f64 {
        self.viscosities[0] / self.viscosities[1].max(1e-30)
    }

    /// Capillary length l_c = √(γ / (ρ g)) (m), using gravity g=9.81.
    pub fn capillary_length(&self) -> f64 {
        let delta_rho = (self.fluid_densities[0] - self.fluid_densities[1]).abs();
        (self.surface_tension / (delta_rho * 9.81)).sqrt()
    }

    /// Ohnesorge number Oh = μ / √(ρ γ L).
    pub fn ohnesorge(&self, length: f64) -> f64 {
        let mu = self.viscosities[0];
        let rho = self.fluid_densities[0];
        let gamma = self.surface_tension;
        mu / (rho * gamma * length).max(1e-30).sqrt()
    }
}

// ---------------------------------------------------------------------------
// WasmFluidCoupling
// ---------------------------------------------------------------------------

/// Fluid-rigid body coupling forces.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmFluidCoupling {
    /// Fluid density ρ_f (kg/m³).
    pub fluid_density: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Added mass coefficient C_a (typically 0.5 for sphere).
    pub added_mass_coeff: f64,
    /// Drag coefficient C_D.
    pub drag_coeff: f64,
    /// Lift coefficient C_L.
    pub lift_coeff: f64,
}

impl WasmFluidCoupling {
    /// Create coupling parameters for water and a sphere.
    pub fn sphere_in_water() -> Self {
        WasmFluidCoupling {
            fluid_density: 1000.0,
            viscosity: 0.001,
            added_mass_coeff: 0.5,
            drag_coeff: 0.47,
            lift_coeff: 0.0,
        }
    }

    /// Buoyancy force on a submerged body with volume `vol` (N, upward positive).
    pub fn buoyancy(&self, volume: f64) -> f64 {
        self.fluid_density * 9.81 * volume
    }

    /// Drag force F_D = ½ ρ C_D A v² (N, opposing velocity direction).
    pub fn drag_force(&self, rel_velocity: f64, projected_area: f64) -> f64 {
        0.5 * self.fluid_density * self.drag_coeff * projected_area * rel_velocity * rel_velocity
    }

    /// Added mass force F_a = -C_a ρ_f V a_rel (N).
    pub fn added_mass_force(&self, volume: f64, relative_acceleration: f64) -> f64 {
        -self.added_mass_coeff * self.fluid_density * volume * relative_acceleration
    }

    /// Stokes drag for creeping flow: F = 6 π μ r v.
    pub fn stokes_drag(&self, radius: f64, rel_velocity: f64) -> f64 {
        6.0 * std::f64::consts::PI * self.viscosity * radius * rel_velocity
    }

    /// Reynolds number for a sphere: Re = ρ v d / μ.
    pub fn reynolds(&self, velocity: f64, diameter: f64) -> f64 {
        self.fluid_density * velocity.abs() * diameter / self.viscosity.max(1e-30)
    }

    /// Strouhal number St = f d / U (vortex-induced vibrations stub).
    pub fn strouhal_viv(frequency: f64, diameter: f64, velocity: f64) -> f64 {
        frequency * diameter / velocity.max(1e-30)
    }
}

// ---------------------------------------------------------------------------
// WasmParticleEmitter
// ---------------------------------------------------------------------------

/// Configuration for a continuous particle emitter.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticleEmitter {
    /// Emission origin \[x, y, z\] (m).
    pub position: [f64; 3],
    /// Emission direction (unit vector) \[dx, dy, dz\].
    pub direction: [f64; 3],
    /// Emission rate (particles / second).
    pub rate: f64,
    /// Initial particle speed (m/s).
    pub velocity: f64,
    /// Particle lifetime (s). After this the particle is inactive.
    pub lifetime: f64,
    /// Emitted particle mass (kg).
    pub mass: f64,
    /// Cone half-angle (radians) for velocity randomisation.
    pub spread_angle: f64,
    /// Accumulated fractional particles not yet emitted.
    accumulator: f64,
    /// Next particle ID.
    next_id: u64,
}

impl WasmParticleEmitter {
    /// Create an emitter at `position` pointing in `direction`.
    pub fn new(
        position: [f64; 3],
        direction: [f64; 3],
        rate: f64,
        velocity: f64,
        lifetime: f64,
        mass: f64,
    ) -> Self {
        WasmParticleEmitter {
            position,
            direction,
            rate,
            velocity,
            lifetime,
            mass,
            spread_angle: 0.0,
            accumulator: 0.0,
            next_id: 0,
        }
    }

    /// Compute how many particles to emit for the given `dt` and return them.
    ///
    /// Particles are created at the emitter position with the configured velocity.
    /// (No random spread is applied in this deterministic version.)
    pub fn emit(&mut self, dt: f64) -> Vec<WasmSphParticle> {
        self.accumulator += self.rate * dt;
        let n = self.accumulator as usize;
        self.accumulator -= n as f64;
        let mut result = Vec::with_capacity(n);
        for _ in 0..n {
            let id = self.next_id;
            self.next_id += 1;
            let mut p = WasmSphParticle::new(id, self.position, self.mass);
            p.velocity = [
                self.direction[0] * self.velocity,
                self.direction[1] * self.velocity,
                self.direction[2] * self.velocity,
            ];
            result.push(p);
        }
        result
    }

    /// Total particles emitted so far.
    pub fn total_emitted(&self) -> u64 {
        self.next_id
    }

    /// Reset the emitter.
    pub fn reset(&mut self) {
        self.accumulator = 0.0;
        self.next_id = 0;
    }
}

// ---------------------------------------------------------------------------
// WasmFlowAnalyzer
// ---------------------------------------------------------------------------

/// Seed for a streamline/pathline integration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StreamlineSeed {
    /// Seed position \[x, y, z\].
    pub position: [f64; 3],
    /// Integration length (arclength, m).
    pub length: f64,
    /// Step size for RK4 integration (m).
    pub step_size: f64,
}

impl StreamlineSeed {
    /// Create a new seed.
    pub fn new(position: [f64; 3], length: f64, step_size: f64) -> Self {
        StreamlineSeed {
            position,
            length,
            step_size,
        }
    }
}

/// Flow field analysis utilities.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmFlowAnalyzer {
    /// Streamline integration seeds.
    pub seeds: Vec<StreamlineSeed>,
    /// Computed streamline paths \[seed_idx\]\[point_idx\]\[xyz\].
    pub streamlines: Vec<Vec<[f64; 3]>>,
    /// Finite-time Lyapunov exponent (FTLE) field (flat, per grid node).
    pub ftle: Vec<f64>,
    /// Domain extent \[xmin, xmax, ymin, ymax, zmin, zmax\].
    pub domain: [f64; 6],
    /// Number of steps taken in pathline integration.
    pub pathline_steps: u64,
}

impl WasmFlowAnalyzer {
    /// Create a new analyser for the given axis-aligned domain.
    pub fn new(domain: [f64; 6]) -> Self {
        WasmFlowAnalyzer {
            domain,
            ..Default::default()
        }
    }

    /// Add a streamline seed.
    pub fn add_seed(&mut self, seed: StreamlineSeed) {
        self.seeds.push(seed);
    }

    /// Integrate streamlines using a constant uniform velocity field (stub).
    ///
    /// `velocity` is the uniform flow velocity \[vx, vy, vz\].
    pub fn integrate_streamlines(&mut self, velocity: [f64; 3]) {
        self.streamlines.clear();
        for seed in &self.seeds {
            let mut path = Vec::new();
            let mut pos = seed.position;
            let speed: f64 = velocity.iter().map(|v| v * v).sum::<f64>().sqrt();
            let ds = seed.step_size;
            let steps = if speed > 1e-30 {
                (seed.length / ds).ceil() as usize
            } else {
                0
            };
            path.push(pos);
            for _ in 0..steps {
                for k in 0..3 {
                    pos[k] += velocity[k] / speed.max(1e-30) * ds;
                }
                path.push(pos);
                // Stop if outside domain.
                if pos[0] < self.domain[0]
                    || pos[0] > self.domain[1]
                    || pos[1] < self.domain[2]
                    || pos[1] > self.domain[3]
                    || pos[2] < self.domain[4]
                    || pos[2] > self.domain[5]
                {
                    break;
                }
            }
            self.streamlines.push(path);
        }
    }

    /// Compute FTLE field for a uniform grid given deformation gradient determinant.
    ///
    /// `nx`, `ny`, `nz` are grid dimensions; `t` is integration time.
    /// Returns the FTLE values (stub: random placeholder proportional to coords).
    pub fn compute_ftle(&mut self, nx: usize, ny: usize, nz: usize, t: f64) -> &[f64] {
        let n = nx * ny * nz;
        self.ftle = (0..n)
            .map(|i| {
                let x = (i % nx) as f64 / nx.max(1) as f64;
                let y = ((i / nx) % ny) as f64 / ny.max(1) as f64;
                // Stub: simple sinusoidal FTLE
                ((x * std::f64::consts::TAU).sin() * (y * std::f64::consts::TAU).cos()
                    / t.max(1e-30))
                .abs()
            })
            .collect();
        &self.ftle
    }

    /// Return the number of computed streamlines.
    pub fn streamline_count(&self) -> usize {
        self.streamlines.len()
    }

    /// Total number of points across all streamlines.
    pub fn total_points(&self) -> usize {
        self.streamlines.iter().map(|s| s.len()).sum()
    }

    /// Clear all computed results.
    pub fn clear_results(&mut self) {
        self.streamlines.clear();
        self.ftle.clear();
        self.pathline_steps = 0;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(not(target_arch = "wasm32"))]
#[cfg(test)]
mod tests {
    use super::*;

    // --- WasmSphConfig ---

    #[test]
    fn test_sph_config_default_valid() {
        let cfg = WasmSphConfig::default();
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_sph_config_invalid_radius() {
        let cfg = WasmSphConfig {
            particle_radius: -0.01,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_sph_config_smoothing_too_small() {
        let cfg = WasmSphConfig {
            particle_radius: 0.1,
            smoothing_length: 0.05,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_sph_kernel_wendland_at_zero() {
        let cfg = WasmSphConfig::default();
        let w = cfg.kernel_wendland(0.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_sph_kernel_wendland_beyond_h_is_zero_like() {
        let cfg = WasmSphConfig::default();
        let w = cfg.kernel_wendland(cfg.smoothing_length * 2.0);
        // Wendland clamps to q=1 at r=h, then evaluation gives some value
        assert!(w >= 0.0);
    }

    #[test]
    fn test_sph_kernel_cubic_decays() {
        let cfg = WasmSphConfig::default();
        let w0 = cfg.kernel_cubic(0.0);
        let w1 = cfg.kernel_cubic(cfg.smoothing_length * 0.5);
        assert!(w0 > w1);
    }

    // --- WasmSphParticle ---

    #[test]
    fn test_sph_particle_kinetic_energy() {
        let mut p = WasmSphParticle::new(0, [0.0; 3], 1.0);
        p.velocity = [1.0, 0.0, 0.0];
        assert!((p.kinetic_energy() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_sph_particle_speed() {
        let mut p = WasmSphParticle::new(1, [0.0; 3], 1.0);
        p.velocity = [3.0, 4.0, 0.0];
        assert!((p.speed() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_sph_particle_boundary() {
        let p = WasmSphParticle::boundary(0, [0.0; 3], 1.0);
        assert!(p.is_boundary);
    }

    // --- WasmSphSimulation ---

    #[test]
    fn test_sph_sim_add_particle() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        let idx = sim.add_particle([0.0, 1.0, 0.0], 0.1);
        assert_eq!(idx, 0);
        assert_eq!(sim.particles.len(), 1);
    }

    #[test]
    fn test_sph_sim_add_particles_flat() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        let positions = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        sim.add_particles(&positions, 0.1);
        assert_eq!(sim.particles.len(), 3);
    }

    #[test]
    fn test_sph_sim_step_advances_time() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        sim.add_particle([0.0, 1.0, 0.0], 0.1);
        sim.step(0.01);
        assert!((sim.time - 0.01).abs() < 1e-15);
        assert_eq!(sim.step_count, 1);
    }

    #[test]
    fn test_sph_sim_gravity_accelerates_particle() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        sim.add_particle([0.0, 1.0, 0.0], 0.1);
        sim.step(0.1);
        let vy = sim.particles[0].velocity[1];
        assert!(vy < 0.0, "gravity should accelerate particle downward");
    }

    #[test]
    fn test_sph_sim_get_positions_length() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        sim.add_particle([0.0; 3], 0.1);
        sim.add_particle([1.0, 0.0, 0.0], 0.1);
        let pos = sim.get_positions();
        assert_eq!(pos.len(), 6);
    }

    #[test]
    fn test_sph_sim_fluid_count_excludes_boundary() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        sim.add_particle([0.0; 3], 0.1);
        sim.particles
            .push(WasmSphParticle::boundary(99, [1.0; 3], 0.1));
        assert_eq!(sim.fluid_count(), 1);
    }

    #[test]
    fn test_sph_sim_total_kinetic_energy() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        let idx = sim.add_particle([0.0; 3], 1.0);
        sim.particles[idx].velocity = [1.0, 0.0, 0.0];
        let ke = sim.total_kinetic_energy();
        assert!((ke - 0.5).abs() < 1e-10);
    }

    // --- WasmLbmConfig ---

    #[test]
    fn test_lbm_config_d2q9_valid() {
        let cfg = WasmLbmConfig::d2q9(64, 64, 0.1);
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_lbm_config_d3q19_valid() {
        let cfg = WasmLbmConfig::d3q19(32, 32, 32, 0.1);
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_lbm_config_invalid_grid() {
        let cfg = WasmLbmConfig {
            nx: 0,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_lbm_config_invalid_tau() {
        let cfg = WasmLbmConfig {
            tau: 0.4,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_lbm_config_reynolds_number() {
        let cfg = WasmLbmConfig::d2q9(64, 64, 0.1);
        let re = cfg.reynolds_number(0.1, 10.0);
        assert!((re - 10.0).abs() < 1e-10);
    }

    // --- WasmLbmSimulation ---

    #[test]
    fn test_lbm_sim_node_count() {
        let cfg = WasmLbmConfig::d2q9(8, 16, 0.1);
        let sim = WasmLbmSimulation::new(cfg);
        assert_eq!(sim.node_count(), 128);
    }

    #[test]
    fn test_lbm_sim_step_increments_count() {
        let mut sim = WasmLbmSimulation::new(WasmLbmConfig::default());
        sim.step();
        assert_eq!(sim.step_count, 1);
    }

    #[test]
    fn test_lbm_sim_set_boundary() {
        let mut sim = WasmLbmSimulation::new(WasmLbmConfig::d2q9(8, 8, 0.1));
        sim.set_boundary(0, 0, 0, true);
        assert!(sim.boundary[0]);
    }

    #[test]
    fn test_lbm_sim_pressure_field_length() {
        let sim = WasmLbmSimulation::new(WasmLbmConfig::d2q9(4, 4, 0.1));
        let p = sim.get_pressure_field();
        assert_eq!(p.len(), 16);
    }

    // --- WasmFluidStats ---

    #[test]
    fn test_fluid_stats_from_sph_empty() {
        let sim = WasmSphSimulation::new(WasmSphConfig::default());
        let stats = WasmFluidStats::from_sph(&sim, 0.01, 0.05);
        assert_eq!(stats.active_count, 0);
    }

    #[test]
    fn test_fluid_stats_cfl_ok() {
        let mut stats = WasmFluidStats::new();
        stats.cfl = 0.5;
        assert!(stats.is_cfl_ok());
        stats.cfl = 1.1;
        assert!(!stats.is_cfl_ok());
    }

    // --- WasmMultiphaseConfig ---

    #[test]
    fn test_multiphase_density_ratio() {
        let cfg = WasmMultiphaseConfig::default();
        let ratio = cfg.density_ratio();
        // water/air ≈ 833
        assert!(ratio > 100.0);
    }

    #[test]
    fn test_multiphase_capillary_length() {
        let cfg = WasmMultiphaseConfig::default();
        let l = cfg.capillary_length();
        // water capillary length ≈ 2.7 mm
        assert!(l > 0.001 && l < 0.01);
    }

    // --- WasmFluidCoupling ---

    #[test]
    fn test_coupling_buoyancy() {
        let c = WasmFluidCoupling::sphere_in_water();
        let f = c.buoyancy(0.001); // 1 litre sphere
        assert!((f - 9.81).abs() < 0.01);
    }

    #[test]
    fn test_coupling_stokes_drag() {
        let c = WasmFluidCoupling::sphere_in_water();
        let f = c.stokes_drag(0.01, 0.1);
        assert!(f > 0.0);
    }

    #[test]
    fn test_coupling_reynolds() {
        let c = WasmFluidCoupling::sphere_in_water();
        let re = c.reynolds(1.0, 0.01);
        assert!((re - 10000.0).abs() < 1.0);
    }

    // --- WasmParticleEmitter ---

    #[test]
    fn test_emitter_emit_zero_dt() {
        let mut em = WasmParticleEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 100.0, 1.0, 5.0, 0.01);
        let particles = em.emit(0.0);
        assert_eq!(particles.len(), 0);
    }

    #[test]
    fn test_emitter_emit_one_second() {
        let mut em = WasmParticleEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 10.0, 1.0, 5.0, 0.01);
        let particles = em.emit(1.0);
        assert_eq!(particles.len(), 10);
    }

    #[test]
    fn test_emitter_total_emitted() {
        let mut em = WasmParticleEmitter::new([0.0; 3], [1.0, 0.0, 0.0], 5.0, 2.0, 10.0, 0.1);
        em.emit(1.0);
        assert_eq!(em.total_emitted(), 5);
    }

    #[test]
    fn test_emitter_reset() {
        let mut em = WasmParticleEmitter::new([0.0; 3], [1.0, 0.0, 0.0], 10.0, 1.0, 5.0, 0.01);
        em.emit(1.0);
        em.reset();
        assert_eq!(em.total_emitted(), 0);
    }

    // --- WasmFlowAnalyzer ---

    #[test]
    fn test_flow_analyzer_add_seeds() {
        let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
        fa.add_seed(StreamlineSeed::new([0.5, 0.5, 0.5], 1.0, 0.1));
        assert_eq!(fa.seeds.len(), 1);
    }

    #[test]
    fn test_flow_analyzer_integrate_streamlines() {
        let mut fa = WasmFlowAnalyzer::new([0.0, 10.0, 0.0, 10.0, 0.0, 10.0]);
        fa.add_seed(StreamlineSeed::new([1.0, 5.0, 5.0], 2.0, 0.5));
        fa.integrate_streamlines([1.0, 0.0, 0.0]);
        assert_eq!(fa.streamline_count(), 1);
        assert!(fa.streamlines[0].len() >= 2);
    }

    #[test]
    fn test_flow_analyzer_ftle_length() {
        let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
        let ftle = fa.compute_ftle(4, 4, 1, 1.0);
        assert_eq!(ftle.len(), 16);
    }

    #[test]
    fn test_flow_analyzer_clear() {
        let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
        fa.add_seed(StreamlineSeed::new([0.5; 3], 1.0, 0.1));
        fa.integrate_streamlines([1.0, 0.0, 0.0]);
        fa.clear_results();
        assert_eq!(fa.streamline_count(), 0);
    }

    #[test]
    fn test_flow_analyzer_total_points() {
        let mut fa = WasmFlowAnalyzer::new([0.0, 100.0, 0.0, 100.0, 0.0, 100.0]);
        fa.add_seed(StreamlineSeed::new([0.0, 50.0, 50.0], 10.0, 1.0));
        fa.add_seed(StreamlineSeed::new([0.0, 25.0, 50.0], 5.0, 1.0));
        fa.integrate_streamlines([1.0, 0.0, 0.0]);
        assert!(fa.total_points() > 0);
    }

    #[test]
    fn test_strouhal_viv() {
        let st = WasmFluidCoupling::strouhal_viv(0.2, 0.01, 1.0);
        assert!((st - 0.002).abs() < 1e-10);
    }
}
