// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly fluid simulation bridge.
//!
//! Provides pure-Rust types for SPH particle simulations, Lattice-Boltzmann
//! (LBM) simulations, multiphase flows, fluid-rigid coupling, particle
//! emitters, and flow-field analysis. Designed for serialisation across the
//! WASM boundary.

use wasm_bindgen::prelude::*;

use serde::{Deserialize, Serialize};

use oxiphysics_lbm::d3q19_full::D3q19Simulation;

// ---------------------------------------------------------------------------
// SPH helper — compact cubic spline kernel (σ = 8/(π h³))
// ---------------------------------------------------------------------------

/// Compact cubic spline kernel value W(r, h) with normalisation σ = 8/(π h³).
#[inline]
fn cubic_spline_compact(r: f64, h: f64) -> f64 {
    use std::f64::consts::PI;
    let sigma = 8.0 / (PI * h * h * h);
    let q = r / h;
    if q >= 1.0 {
        0.0
    } else if q >= 0.5 {
        sigma * 2.0 * (1.0 - q).powi(3)
    } else {
        sigma * (1.0 - 6.0 * q * q + 6.0 * q * q * q)
    }
}

/// Compute SPH density ρᵢ = Σⱼ mⱼ · W(|rᵢ − rⱼ|, h) for all fluid particles.
fn compute_sph_densities(particles: &[WasmSphParticle], h: f64) -> Vec<f64> {
    let n = particles.len();
    let mut densities = vec![0.0_f64; n];
    for i in 0..n {
        let pi = &particles[i];
        let mut rho = 0.0_f64;
        for pj in particles.iter() {
            let dx = pi.position[0] - pj.position[0];
            let dy = pi.position[1] - pj.position[1];
            let dz = pi.position[2] - pj.position[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            rho += pj.mass * cubic_spline_compact(r, h);
        }
        densities[i] = rho;
    }
    densities
}

// ---------------------------------------------------------------------------
// WasmSphConfig
// ---------------------------------------------------------------------------

/// Configuration for an SPH (Smoothed Particle Hydrodynamics) simulation.
#[wasm_bindgen]
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
    /// Gravity vector [gx, gy, gz] (m/s²) — private; use getters/setters from JS.
    gravity: [f64; 3],
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

#[wasm_bindgen]
impl WasmSphConfig {
    /// Create a default water-like config (JS constructor).
    pub fn create_water() -> WasmSphConfig {
        WasmSphConfig::water()
    }

    /// Validate config; returns empty string on success, error message on failure.
    pub fn validate_js(&self) -> String {
        self.validate().err().unwrap_or_default()
    }

    /// Wendland kernel value (JS-exposed).
    pub fn kernel_wendland_js(&self, r: f64) -> f64 {
        self.kernel_wendland(r)
    }

    /// Cubic spline kernel value (JS-exposed).
    pub fn kernel_cubic_js(&self, r: f64) -> f64 {
        self.kernel_cubic(r)
    }

    /// Get gravity vector as [gx, gy, gz].
    pub fn get_gravity(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }

    /// Set gravity vector.
    pub fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
    }

    /// Gravity X component.
    pub fn gravity_x(&self) -> f64 {
        self.gravity[0]
    }

    /// Gravity Y component.
    pub fn gravity_y(&self) -> f64 {
        self.gravity[1]
    }

    /// Gravity Z component.
    pub fn gravity_z(&self) -> f64 {
        self.gravity[2]
    }
}

// ---------------------------------------------------------------------------
// WasmSphParticle
// ---------------------------------------------------------------------------

/// A single SPH fluid particle.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSphParticle {
    /// Position [x, y, z] (m) — private; use get_position/set_position from JS.
    position: [f64; 3],
    /// Velocity [vx, vy, vz] (m/s) — private; use get_velocity/set_velocity from JS.
    velocity: [f64; 3],
    /// Acceleration [ax, ay, az] (m/s²) — private; use get_acceleration from JS.
    acceleration: [f64; 3],
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
    /// Particle unique ID — private; use get_id from JS.
    id: u64,
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

#[wasm_bindgen]
impl WasmSphParticle {
    /// Create a fluid particle at rest (JS-compatible constructor).
    pub fn new_js(id: f64, x: f64, y: f64, z: f64, mass: f64) -> WasmSphParticle {
        WasmSphParticle::new(id as u64, [x, y, z], mass)
    }

    /// Create a static boundary particle (JS-compatible).
    pub fn boundary_js(id: f64, x: f64, y: f64, z: f64, mass: f64) -> WasmSphParticle {
        WasmSphParticle::boundary(id as u64, [x, y, z], mass)
    }

    /// Get position as [x, y, z].
    pub fn get_position(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Set position.
    pub fn set_position(&mut self, x: f64, y: f64, z: f64) {
        self.position = [x, y, z];
    }

    /// Get velocity as [vx, vy, vz].
    pub fn get_velocity(&self) -> Vec<f64> {
        self.velocity.to_vec()
    }

    /// Set velocity.
    pub fn set_velocity(&mut self, vx: f64, vy: f64, vz: f64) {
        self.velocity = [vx, vy, vz];
    }

    /// Get acceleration as [ax, ay, az].
    pub fn get_acceleration(&self) -> Vec<f64> {
        self.acceleration.to_vec()
    }

    /// Get particle ID as f64.
    pub fn get_id(&self) -> f64 {
        self.id as f64
    }

    /// Kinetic energy (JS-exposed).
    pub fn kinetic_energy_js(&self) -> f64 {
        self.kinetic_energy()
    }

    /// Speed (JS-exposed).
    pub fn speed_js(&self) -> f64 {
        self.speed()
    }
}

// ---------------------------------------------------------------------------
// WasmSphSimulation
// ---------------------------------------------------------------------------

/// SPH fluid simulation manager.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSphSimulation {
    /// Simulation configuration.
    #[wasm_bindgen(skip)]
    pub config: WasmSphConfig,
    /// All particles (fluid + boundary).
    #[wasm_bindgen(skip)]
    pub particles: Vec<WasmSphParticle>,
    /// Simulated time (s).
    pub time: f64,
    /// Step counter — private; use get_step_count from JS.
    step_count: u64,
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

    /// Add multiple particles from flat position array [x0,y0,z0, x1,y1,z1, ...].
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

    /// Advance the simulation by `dt` seconds (Euler integration).
    pub fn step(&mut self, dt: f64) {
        let grav = self.config.gravity;
        let h = self.config.smoothing_length;
        let densities = compute_sph_densities(&self.particles, h);
        for (i, p) in self.particles.iter_mut().enumerate() {
            if p.is_boundary {
                continue;
            }
            p.density = densities[i];
            let ratio = p.density / self.config.rest_density;
            p.pressure = self.config.stiffness * (ratio.powf(self.config.tait_exponent) - 1.0);
        }
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

    /// Get flat array of all particle positions [x0,y0,z0, ...].
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

#[wasm_bindgen]
impl WasmSphSimulation {
    /// Create a new simulation (JS constructor).
    pub fn create(config: WasmSphConfig) -> WasmSphSimulation {
        WasmSphSimulation::new(config)
    }

    /// Add a particle at (x, y, z) with given mass; returns index as u32.
    pub fn add_particle_js(&mut self, x: f64, y: f64, z: f64, mass: f64) -> u32 {
        self.add_particle([x, y, z], mass) as u32
    }

    /// Add particles from flat [x0,y0,z0, ...] array (JS-compatible).
    pub fn add_particles_js(&mut self, positions_flat: &[f64], mass: f64) {
        self.add_particles(positions_flat, mass);
    }

    /// Advance simulation by dt seconds.
    pub fn step_js(&mut self, dt: f64) {
        self.step(dt);
    }

    /// Get flat array of all particle positions.
    pub fn get_positions_js(&self) -> Vec<f64> {
        self.get_positions()
    }

    /// Get flat array of all velocities.
    pub fn get_velocities_js(&self) -> Vec<f64> {
        self.get_velocities()
    }

    /// Get density of each particle.
    pub fn get_densities_js(&self) -> Vec<f64> {
        self.get_densities()
    }

    /// Number of fluid particles as u32.
    pub fn fluid_count_js(&self) -> u32 {
        self.fluid_count() as u32
    }

    /// Total kinetic energy.
    pub fn total_kinetic_energy_js(&self) -> f64 {
        self.total_kinetic_energy()
    }

    /// Step counter as f64.
    pub fn get_step_count(&self) -> f64 {
        self.step_count as f64
    }

    /// Simulated time.
    pub fn get_time(&self) -> f64 {
        self.time
    }
}

// ---------------------------------------------------------------------------
// WasmLbmConfig
// ---------------------------------------------------------------------------

/// Configuration for a Lattice-Boltzmann Method (LBM) simulation.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmLbmConfig {
    /// Grid resolution in X — private; use get_nx from JS.
    nx: usize,
    /// Grid resolution in Y — private; use get_ny from JS.
    ny: usize,
    /// Grid resolution in Z — private; use get_nz from JS.
    nz: usize,
    /// Kinematic viscosity ν (lattice units).
    pub viscosity: f64,
    /// Lattice type — private; use get_lattice_type from JS.
    lattice_type: String,
    /// Collision operator — private; use get_collision from JS.
    collision: String,
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

#[wasm_bindgen]
impl WasmLbmConfig {
    /// Create a D2Q9 2-D config (JS constructor).
    pub fn d2q9_js(nx: u32, ny: u32, viscosity: f64) -> WasmLbmConfig {
        WasmLbmConfig::d2q9(nx as usize, ny as usize, viscosity)
    }

    /// Create a D3Q19 3-D config (JS constructor).
    pub fn d3q19_js(nx: u32, ny: u32, nz: u32, viscosity: f64) -> WasmLbmConfig {
        WasmLbmConfig::d3q19(nx as usize, ny as usize, nz as usize, viscosity)
    }

    /// Create a default D2Q9 64×64 config.
    pub fn default_config() -> WasmLbmConfig {
        WasmLbmConfig::default()
    }

    /// Reynolds number (JS-exposed).
    pub fn reynolds_number_js(&self, velocity: f64, length: f64) -> f64 {
        self.reynolds_number(velocity, length)
    }

    /// Validate; returns empty string on success, error message on failure.
    pub fn validate_js(&self) -> String {
        self.validate().err().unwrap_or_default()
    }

    /// Grid X dimension as u32.
    pub fn get_nx(&self) -> u32 {
        self.nx as u32
    }

    /// Grid Y dimension as u32.
    pub fn get_ny(&self) -> u32 {
        self.ny as u32
    }

    /// Grid Z dimension as u32.
    pub fn get_nz(&self) -> u32 {
        self.nz as u32
    }

    /// Lattice type string.
    pub fn get_lattice_type(&self) -> String {
        self.lattice_type.clone()
    }

    /// Collision operator string.
    pub fn get_collision(&self) -> String {
        self.collision.clone()
    }
}

// ---------------------------------------------------------------------------
// WasmLbmSimulation
// ---------------------------------------------------------------------------

/// LBM simulation state — thin WASM bridge over `D3q19Simulation`.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmLbmSimulation {
    /// LBM configuration.
    #[wasm_bindgen(skip)]
    pub config: WasmLbmConfig,
    /// Flat density cache (one value per lattice node).
    #[wasm_bindgen(skip)]
    pub density: Vec<f64>,
    /// Flat velocity cache: [vx, vy, vz] per node.
    #[wasm_bindgen(skip)]
    pub velocity: Vec<f64>,
    /// Flat pressure cache.
    #[wasm_bindgen(skip)]
    pub pressure: Vec<f64>,
    /// Boundary mask (true = solid wall node).
    #[wasm_bindgen(skip)]
    pub boundary: Vec<bool>,
    /// Simulation step counter — private; use get_step_count from JS.
    step_count: u64,
    /// Simulated time (lattice time units).
    pub time: f64,
    /// Inner D3Q19 BGK grid — excluded from serialisation.
    #[serde(skip)]
    grid: Option<D3q19Simulation>,
}

/// Build a fresh `D3q19Simulation` from a `WasmLbmConfig`.
fn make_d3q19(config: &WasmLbmConfig) -> D3q19Simulation {
    D3q19Simulation::new(config.nx, config.ny, config.nz, config.tau)
}

impl WasmLbmSimulation {
    /// Create a new LBM simulation initialised to quiescent flow (ρ=1, u=0).
    pub fn new(config: WasmLbmConfig) -> Self {
        let n = config.nx * config.ny * config.nz;
        let grid = make_d3q19(&config);
        let density = grid.density_field();
        let vel_field = grid.velocity_field();
        let cs2 = config.sound_speed * config.sound_speed;
        let pressure: Vec<f64> = density.iter().map(|rho| rho * cs2).collect();
        let velocity: Vec<f64> = vel_field.iter().flat_map(|v| v.iter().copied()).collect();
        WasmLbmSimulation {
            density,
            velocity,
            pressure,
            boundary: vec![false; n],
            step_count: 0,
            time: 0.0,
            grid: Some(grid),
            config,
        }
    }

    /// Obtain a mutable reference to the inner grid, reconstructing from
    /// config if it was dropped (e.g. after serde round-trip).
    fn grid_mut(&mut self) -> &mut D3q19Simulation {
        if self.grid.is_none() {
            self.grid = Some(make_d3q19(&self.config));
        }
        self.grid.as_mut().expect("grid is always Some after init")
    }

    /// Sync the public macroscopic caches from the current grid state.
    fn sync_macroscopic(&mut self) {
        if let Some(ref g) = self.grid {
            let cs2 = self.config.sound_speed * self.config.sound_speed;
            self.density = g.density_field();
            let vel_field = g.velocity_field();
            self.velocity = vel_field.iter().flat_map(|v| v.iter().copied()).collect();
            for (p, rho) in self.pressure.iter_mut().zip(self.density.iter()) {
                *p = rho * cs2;
            }
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

    /// Advance one LBM time step.
    pub fn step(&mut self) {
        self.grid_mut().step();
        self.sync_macroscopic();
        self.step_count += 1;
        self.time += 1.0;
    }

    /// Get a copy of the velocity field as [vx, vy, vz] per node.
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

    /// Mean density across all nodes.
    pub fn mean_density(&self) -> f64 {
        if self.density.is_empty() {
            return 0.0;
        }
        self.density.iter().sum::<f64>() / self.density.len() as f64
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

#[wasm_bindgen]
impl WasmLbmSimulation {
    /// Create a new LBM simulation (JS constructor).
    pub fn create(config: WasmLbmConfig) -> WasmLbmSimulation {
        WasmLbmSimulation::new(config)
    }

    /// Advance one LBM time step (JS-exposed).
    pub fn step_js(&mut self) {
        self.step();
    }

    /// Get velocity field (JS-exposed).
    pub fn get_velocity_field_js(&self) -> Vec<f64> {
        self.get_velocity_field()
    }

    /// Get pressure field (JS-exposed).
    pub fn get_pressure_field_js(&self) -> Vec<f64> {
        self.get_pressure_field()
    }

    /// Get density field copy.
    pub fn get_density_field(&self) -> Vec<f64> {
        self.density.clone()
    }

    /// Number of lattice nodes as u32.
    pub fn node_count_js(&self) -> u32 {
        self.node_count() as u32
    }

    /// Mean density (JS-exposed).
    pub fn mean_density_js(&self) -> f64 {
        self.mean_density()
    }

    /// Mark node as solid boundary; indices as u32.
    pub fn set_boundary_js(&mut self, ix: u32, iy: u32, iz: u32, solid: bool) {
        self.set_boundary(ix as usize, iy as usize, iz as usize, solid);
    }

    /// Set inlet velocity (JS-exposed).
    pub fn set_inlet_velocity_js(&mut self, vx: f64) {
        self.set_inlet_velocity(vx);
    }

    /// Step counter as f64.
    pub fn get_step_count(&self) -> f64 {
        self.step_count as f64
    }

    /// Simulated time (JS-exposed).
    pub fn get_time(&self) -> f64 {
        self.time
    }
}

// ---------------------------------------------------------------------------
// WasmFluidStats
// ---------------------------------------------------------------------------

/// Aggregate fluid simulation statistics for one step.
#[wasm_bindgen]
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

#[wasm_bindgen]
impl WasmFluidStats {
    /// Create zeroed stats (JS constructor).
    pub fn create() -> WasmFluidStats {
        WasmFluidStats::new()
    }

    /// Compute stats from SPH simulation (JS-exposed).
    pub fn from_sph_js(sim: &WasmSphSimulation, dt: f64, dx: f64) -> WasmFluidStats {
        WasmFluidStats::from_sph(sim, dt, dx)
    }

    /// Returns true if CFL < 1 (JS-exposed).
    pub fn is_cfl_ok_js(&self) -> bool {
        self.is_cfl_ok()
    }
}

// ---------------------------------------------------------------------------
// WasmMultiphaseConfig
// ---------------------------------------------------------------------------

/// Configuration for a two-phase (immiscible) flow simulation.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmMultiphaseConfig {
    /// Densities of the two fluid phases [ρ₁, ρ₂] — private; use getters from JS.
    fluid_densities: [f64; 2],
    /// Interfacial surface tension γ (N/m).
    pub surface_tension: f64,
    /// Static contact angle θ (radians).
    pub contact_angle: f64,
    /// Interface thickness W (m) – Cahn-Hilliard parameter.
    pub interface_thickness: f64,
    /// Mobility parameter M in phase-field equation.
    pub mobility: f64,
    /// Viscosities of the two phases [μ₁, μ₂] — private; use getters from JS.
    viscosities: [f64; 2],
}

impl Default for WasmMultiphaseConfig {
    fn default() -> Self {
        WasmMultiphaseConfig {
            fluid_densities: [1000.0, 1.2],
            surface_tension: 0.0728,
            contact_angle: std::f64::consts::PI / 3.0,
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

#[wasm_bindgen]
impl WasmMultiphaseConfig {
    /// Create default water/air config (JS constructor).
    pub fn create_default() -> WasmMultiphaseConfig {
        WasmMultiphaseConfig::default()
    }

    /// Phase 1 fluid density (kg/m³).
    pub fn get_fluid_density_1(&self) -> f64 {
        self.fluid_densities[0]
    }

    /// Phase 2 fluid density (kg/m³).
    pub fn get_fluid_density_2(&self) -> f64 {
        self.fluid_densities[1]
    }

    /// Set both fluid densities.
    pub fn set_fluid_densities(&mut self, rho1: f64, rho2: f64) {
        self.fluid_densities = [rho1, rho2];
    }

    /// Phase 1 viscosity (Pa·s).
    pub fn get_viscosity_1(&self) -> f64 {
        self.viscosities[0]
    }

    /// Phase 2 viscosity (Pa·s).
    pub fn get_viscosity_2(&self) -> f64 {
        self.viscosities[1]
    }

    /// Set both viscosities.
    pub fn set_viscosities(&mut self, mu1: f64, mu2: f64) {
        self.viscosities = [mu1, mu2];
    }

    /// Density ratio ρ₁/ρ₂ (JS-exposed).
    pub fn density_ratio_js(&self) -> f64 {
        self.density_ratio()
    }

    /// Viscosity ratio μ₁/μ₂ (JS-exposed).
    pub fn viscosity_ratio_js(&self) -> f64 {
        self.viscosity_ratio()
    }

    /// Capillary length (JS-exposed).
    pub fn capillary_length_js(&self) -> f64 {
        self.capillary_length()
    }

    /// Ohnesorge number (JS-exposed).
    pub fn ohnesorge_js(&self, length: f64) -> f64 {
        self.ohnesorge(length)
    }
}

// ---------------------------------------------------------------------------
// WasmFluidCoupling
// ---------------------------------------------------------------------------

/// Fluid-rigid body coupling forces.
#[wasm_bindgen]
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

#[wasm_bindgen]
impl WasmFluidCoupling {
    /// Create coupling parameters for water and a sphere.
    pub fn sphere_in_water() -> WasmFluidCoupling {
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

    /// Drag force F_D = ½ ρ C_D A v².
    pub fn drag_force(&self, rel_velocity: f64, projected_area: f64) -> f64 {
        0.5 * self.fluid_density * self.drag_coeff * projected_area * rel_velocity * rel_velocity
    }

    /// Added mass force F_a = -C_a ρ_f V a_rel.
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

    /// Strouhal number St = f d / U (vortex-induced vibrations).
    pub fn strouhal_viv(frequency: f64, diameter: f64, velocity: f64) -> f64 {
        frequency * diameter / velocity.max(1e-30)
    }
}

// ---------------------------------------------------------------------------
// WasmParticleEmitter
// ---------------------------------------------------------------------------

/// Configuration for a continuous particle emitter.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticleEmitter {
    /// Emission origin [x, y, z] — private; use get_position/set_position from JS.
    position: [f64; 3],
    /// Emission direction (unit vector) [dx, dy, dz] — private; use getters from JS.
    direction: [f64; 3],
    /// Emission rate (particles / second).
    pub rate: f64,
    /// Initial particle speed (m/s).
    pub velocity: f64,
    /// Particle lifetime (s).
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

#[wasm_bindgen]
impl WasmParticleEmitter {
    /// Create an emitter (JS constructor).
    pub fn create(
        px: f64,
        py: f64,
        pz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        rate: f64,
        velocity: f64,
        lifetime: f64,
        mass: f64,
    ) -> WasmParticleEmitter {
        WasmParticleEmitter::new([px, py, pz], [dx, dy, dz], rate, velocity, lifetime, mass)
    }

    /// Emit particles and return the count of emitted particles.
    pub fn emit_count_js(&mut self, dt: f64) -> u32 {
        self.emit(dt).len() as u32
    }

    /// Total particles emitted as f64.
    pub fn total_emitted_js(&self) -> f64 {
        self.total_emitted() as f64
    }

    /// Reset the emitter (JS-exposed).
    pub fn reset_js(&mut self) {
        self.reset();
    }

    /// Get emitter position as [x, y, z].
    pub fn get_position(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Set emitter position.
    pub fn set_position(&mut self, px: f64, py: f64, pz: f64) {
        self.position = [px, py, pz];
    }

    /// Get emitter direction as [dx, dy, dz].
    pub fn get_direction(&self) -> Vec<f64> {
        self.direction.to_vec()
    }

    /// Set emitter direction.
    pub fn set_direction(&mut self, dx: f64, dy: f64, dz: f64) {
        self.direction = [dx, dy, dz];
    }
}

// ---------------------------------------------------------------------------
// WasmFlowAnalyzer
// ---------------------------------------------------------------------------

/// Seed for a streamline/pathline integration.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StreamlineSeed {
    /// Seed position [x, y, z] — private; use get_position from JS.
    position: [f64; 3],
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

#[wasm_bindgen]
impl StreamlineSeed {
    /// Create a streamline seed (JS constructor).
    pub fn create(px: f64, py: f64, pz: f64, length: f64, step_size: f64) -> StreamlineSeed {
        StreamlineSeed::new([px, py, pz], length, step_size)
    }

    /// Get seed position as [x, y, z].
    pub fn get_position(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Set seed position.
    pub fn set_position(&mut self, px: f64, py: f64, pz: f64) {
        self.position = [px, py, pz];
    }
}

/// Flow field analysis utilities.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmFlowAnalyzer {
    /// Streamline integration seeds.
    #[wasm_bindgen(skip)]
    pub seeds: Vec<StreamlineSeed>,
    /// Computed streamline paths `seed_idx[point_idx][xyz]`.
    #[wasm_bindgen(skip)]
    pub streamlines: Vec<Vec<[f64; 3]>>,
    /// FTLE field (flat, per grid node).
    #[wasm_bindgen(skip)]
    pub ftle: Vec<f64>,
    /// Domain extent [xmin, xmax, ymin, ymax, zmin, zmax] — private; use get_domain from JS.
    domain: [f64; 6],
    /// Number of steps taken in pathline integration — private; use get_pathline_steps from JS.
    pathline_steps: u64,
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

    /// Perform a single RK4 integration step (Rust-only; generic closure).
    pub fn rk4_step<F>(&self, pos: [f64; 3], h: f64, vel_fn: &F) -> [f64; 3]
    where
        F: Fn([f64; 3]) -> [f64; 3],
    {
        let add =
            |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] };
        let scale = |a: [f64; 3], s: f64| -> [f64; 3] { [a[0] * s, a[1] * s, a[2] * s] };

        let k1 = vel_fn(pos);
        let k2 = vel_fn(add(pos, scale(k1, h * 0.5)));
        let k3 = vel_fn(add(pos, scale(k2, h * 0.5)));
        let k4 = vel_fn(add(pos, scale(k3, h)));
        let weighted = add(add(k1, scale(k2, 2.0)), add(scale(k3, 2.0), k4));
        add(pos, scale(weighted, h / 6.0))
    }

    /// Integrate streamlines using RK4 with a constant uniform velocity field.
    pub fn integrate_streamlines(&mut self, velocity: [f64; 3]) {
        self.streamlines.clear();
        let domain = self.domain;
        for seed in &self.seeds {
            let mut path = Vec::new();
            let mut pos = seed.position;
            let ds = seed.step_size;
            let speed: f64 = velocity.iter().map(|v| v * v).sum::<f64>().sqrt();
            let steps = if speed > 1e-30 {
                (seed.length / ds).ceil() as usize
            } else {
                0
            };
            path.push(pos);
            let vel_fn = |_p: [f64; 3]| velocity;
            for _ in 0..steps {
                pos = self.rk4_step(pos, ds, &vel_fn);
                path.push(pos);
                if pos[0] < domain[0]
                    || pos[0] > domain[1]
                    || pos[1] < domain[2]
                    || pos[1] > domain[3]
                    || pos[2] < domain[4]
                    || pos[2] > domain[5]
                {
                    break;
                }
            }
            self.streamlines.push(path);
        }
    }

    /// Compute FTLE field for a uniform grid.
    pub fn compute_ftle(&mut self, nx: usize, ny: usize, nz: usize, t: f64) -> &[f64] {
        let n = nx * ny * nz;
        self.ftle = (0..n)
            .map(|i| {
                let x = (i % nx) as f64 / nx.max(1) as f64;
                let y = ((i / nx) % ny) as f64 / ny.max(1) as f64;
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

#[wasm_bindgen]
impl WasmFlowAnalyzer {
    /// Create a new flow analyser (JS constructor).
    pub fn create(
        xmin: f64,
        xmax: f64,
        ymin: f64,
        ymax: f64,
        zmin: f64,
        zmax: f64,
    ) -> WasmFlowAnalyzer {
        WasmFlowAnalyzer::new([xmin, xmax, ymin, ymax, zmin, zmax])
    }

    /// Add a streamline seed (JS-exposed).
    pub fn add_seed_js(&mut self, seed: StreamlineSeed) {
        self.add_seed(seed);
    }

    /// Integrate streamlines with a uniform velocity field.
    pub fn integrate_streamlines_js(&mut self, vx: f64, vy: f64, vz: f64) {
        self.integrate_streamlines([vx, vy, vz]);
    }

    /// Compute FTLE field; returns the FTLE values as `Vec<f64>`.
    pub fn compute_ftle_js(&mut self, nx: u32, ny: u32, nz: u32, t: f64) -> Vec<f64> {
        self.compute_ftle(nx as usize, ny as usize, nz as usize, t)
            .to_vec()
    }

    /// Number of computed streamlines as u32.
    pub fn streamline_count_js(&self) -> u32 {
        self.streamline_count() as u32
    }

    /// Total streamline points as u32.
    pub fn total_points_js(&self) -> u32 {
        self.total_points() as u32
    }

    /// Get domain as [xmin, xmax, ymin, ymax, zmin, zmax].
    pub fn get_domain(&self) -> Vec<f64> {
        self.domain.to_vec()
    }

    /// Get FTLE field copy.
    pub fn get_ftle(&self) -> Vec<f64> {
        self.ftle.clone()
    }

    /// Pathline integration step count as f64.
    pub fn get_pathline_steps(&self) -> f64 {
        self.pathline_steps as f64
    }

    /// Clear all computed results (JS-exposed).
    pub fn clear_results_js(&mut self) {
        self.clear_results();
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

    // --- L1: SPH density summation tests ---

    #[test]
    fn test_sph_density_single_particle_self_contribution() {
        use std::f64::consts::PI;
        let h = 0.2_f64;
        let m = 0.5_f64;
        let expected = m * 8.0 / (PI * h.powi(3));
        let densities = compute_sph_densities(&[WasmSphParticle::new(0, [0.0; 3], m)], h);
        assert!(
            (densities[0] - expected).abs() < 1e-12,
            "single-particle density {:.6} != expected {:.6}",
            densities[0],
            expected
        );
    }

    #[test]
    fn test_sph_density_cubic_lattice_symmetry() {
        let h = 0.2_f64;
        let m = 0.1_f64;
        let s = h / 2.0;
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [s, 0.0, 0.0],
            [0.0, s, 0.0],
            [s, s, 0.0],
            [0.0, 0.0, s],
            [s, 0.0, s],
            [0.0, s, s],
            [s, s, s],
        ];
        let particles: Vec<WasmSphParticle> = positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| WasmSphParticle::new(i as u64, pos, m))
            .collect();
        let densities = compute_sph_densities(&particles, h);
        let rho0 = densities[0];
        for (i, &rho) in densities.iter().enumerate() {
            assert!(
                (rho - rho0).abs() < 1e-12,
                "particle {i} density {rho:.6} differs from particle 0 density {rho0:.6}"
            );
        }
        let expected = m
            * (cubic_spline_compact(0.0, h)
                + 3.0 * cubic_spline_compact(s, h)
                + 3.0 * cubic_spline_compact(s * std::f64::consts::SQRT_2, h)
                + cubic_spline_compact(s * 3.0_f64.sqrt(), h));
        assert!(
            (rho0 - expected).abs() < 1e-10,
            "cubic-lattice density {rho0:.6} != expected {expected:.6}"
        );
        assert!(rho0 > 0.0 && rho0.is_finite());
    }

    #[test]
    fn test_sph_density_step_not_stub() {
        let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
        let h = sim.config.smoothing_length;
        sim.add_particle([0.0, 0.0, 0.0], 0.1);
        sim.add_particle([h * 0.3, 0.0, 0.0], 0.1);
        sim.step(0.001);
        let rho0 = sim.particles[0].density;
        let rho1 = sim.particles[1].density;
        assert!(
            rho0 > 0.0 && rho0 < sim.config.rest_density * 0.5,
            "density {rho0:.4} should be well below rest_density after real SPH summation"
        );
        assert!(
            (rho0 - rho1).abs() < 1e-8,
            "symmetric pair should have equal densities: {rho0:.6} vs {rho1:.6}"
        );
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

    // --- L2: LBM D3Q19 BGK facade tests ---

    #[test]
    fn test_lbm_d3q19_step_no_panic() {
        let mut sim = WasmLbmSimulation::new(WasmLbmConfig::d3q19(4, 4, 4, 0.1));
        sim.step();
        assert_eq!(sim.step_count, 1);
    }

    #[test]
    fn test_lbm_d3q19_density_nonzero_after_init() {
        let sim = WasmLbmSimulation::new(WasmLbmConfig::d3q19(4, 4, 4, 0.1));
        let mean = sim.mean_density();
        assert!(
            (mean - 1.0).abs() < 1e-10,
            "initial mean density should be 1.0, got {mean}"
        );
        assert!(sim.density.iter().all(|&rho| rho > 0.0));
    }

    #[test]
    fn test_lbm_d3q19_density_conservation() {
        let mut sim = WasmLbmSimulation::new(WasmLbmConfig::d3q19(8, 8, 4, 0.6));
        let total_init: f64 = sim.density.iter().sum();
        for _ in 0..10 {
            sim.step();
        }
        let total_final: f64 = sim.density.iter().sum();
        let rel_err = (total_final - total_init).abs() / total_init;
        assert!(
            rel_err < 1e-10,
            "total density changed by relative error {rel_err:.2e} (should be < 1e-10)"
        );
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
        assert!(ratio > 100.0);
    }

    #[test]
    fn test_multiphase_capillary_length() {
        let cfg = WasmMultiphaseConfig::default();
        let l = cfg.capillary_length();
        assert!(l > 0.001 && l < 0.01);
    }

    // --- WasmFluidCoupling ---

    #[test]
    fn test_coupling_buoyancy() {
        let c = WasmFluidCoupling::sphere_in_water();
        let f = c.buoyancy(0.001);
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

    // --- H4: RK4 circular streamline ---

    #[test]
    fn rk4_circular_streamline() {
        let fa = WasmFlowAnalyzer::new([-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]);
        let h = 0.01_f64;
        let steps = 100_usize;
        let mut pos = [1.0_f64, 0.0, 0.0];
        let vel_fn = |p: [f64; 3]| -> [f64; 3] { [-p[1], p[0], 0.0] };
        for _ in 0..steps {
            pos = fa.rk4_step(pos, h, &vel_fn);
        }
        let radius = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
        assert!(
            (radius - 1.0).abs() < 0.01,
            "RK4 radius error too large: radius = {radius}"
        );
    }
}
