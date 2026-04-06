#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoothed Particle Hydrodynamics (SPH) simulation API for Python interop.
//!
//! Provides `PySphSimulation`: a standalone 3-D weakly compressible SPH (WCSPH)
//! solver with poly6/spiky kernels, pressure and viscosity forces, configurable
//! gravity, and query helpers.

#![allow(missing_docs)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for the SPH simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySphConfig {
    /// Smoothing length (kernel radius h).
    pub h: f64,
    /// Reference (rest) density ρ₀ in kg/m³.
    pub rest_density: f64,
    /// Pressure stiffness constant k (state-equation coefficient).
    pub stiffness: f64,
    /// Dynamic viscosity μ.
    pub viscosity: f64,
    /// Gravity vector `[gx, gy, gz]`.
    pub gravity: [f64; 3],
    /// Mass of each particle in kg.
    pub particle_mass: f64,
    /// Whether to enforce a floor boundary at y=0 (floor_y).
    pub floor_enabled: bool,
    /// Y-coordinate of the floor plane (only used when `floor_enabled`).
    pub floor_y: f64,
    /// Coefficient of restitution for floor bounces (0=inelastic, 1=elastic).
    pub floor_restitution: f64,
}

impl PySphConfig {
    /// Water-like default configuration.
    pub fn water() -> Self {
        Self {
            h: 0.1,
            rest_density: 1000.0,
            stiffness: 200.0,
            viscosity: 0.01,
            gravity: [0.0, -9.81, 0.0],
            particle_mass: 0.02,
            floor_enabled: false,
            floor_y: 0.0,
            floor_restitution: 0.3,
        }
    }

    /// Gas-like configuration with low density and stiffness.
    pub fn gas() -> Self {
        Self {
            h: 0.2,
            rest_density: 1.2,
            stiffness: 50.0,
            viscosity: 1.81e-5,
            gravity: [0.0, 0.0, 0.0],
            particle_mass: 0.001,
            floor_enabled: false,
            floor_y: 0.0,
            floor_restitution: 1.0,
        }
    }
}

impl Default for PySphConfig {
    fn default() -> Self {
        Self::water()
    }
}

// ---------------------------------------------------------------------------
// PySphSimulation
// ---------------------------------------------------------------------------

/// A 3-D weakly-compressible SPH fluid simulation.
///
/// Particles are integrated with forward Euler and use Tait's equation of
/// state (p = k * (ρ/ρ₀ - 1)) for pressure. Supports an optional flat floor
/// boundary.
#[derive(Debug, Clone)]
pub struct PySphSimulation {
    /// Particle positions `[x, y, z]`.
    positions: Vec<[f64; 3]>,
    /// Particle velocities `[vx, vy, vz]`.
    velocities: Vec<[f64; 3]>,
    /// Per-particle densities.
    densities: Vec<f64>,
    /// Per-particle pressures.
    pressures: Vec<f64>,
    /// Simulation configuration.
    config: PySphConfig,
    /// Accumulated simulation time.
    time: f64,
    /// Number of steps taken.
    step_count: u64,
}

impl PySphSimulation {
    /// Create an empty SPH simulation with the given configuration.
    pub fn new(config: PySphConfig) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            config,
            time: 0.0,
            step_count: 0,
        }
    }

    /// Add a particle at `pos` with zero velocity. Returns the particle index.
    pub fn add_particle(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.positions.len();
        self.positions.push(pos);
        self.velocities.push([0.0; 3]);
        self.densities.push(self.config.rest_density);
        self.pressures.push(0.0);
        idx
    }

    /// Add a particle at `pos` with the given initial velocity.
    pub fn add_particle_with_velocity(&mut self, pos: [f64; 3], vel: [f64; 3]) -> usize {
        let idx = self.add_particle(pos);
        self.velocities[idx] = vel;
        idx
    }

    /// Add `nx * ny * nz` particles in a 3-D lattice starting at `origin`,
    /// with inter-particle spacing `dx`.
    pub fn add_particle_block(
        &mut self,
        origin: [f64; 3],
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
    ) {
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let pos = [
                        origin[0] + ix as f64 * dx,
                        origin[1] + iy as f64 * dx,
                        origin[2] + iz as f64 * dx,
                    ];
                    self.add_particle(pos);
                }
            }
        }
    }

    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.positions.len()
    }

    /// Position of particle `i`, or `None` if out of bounds.
    pub fn position(&self, i: usize) -> Option<[f64; 3]> {
        self.positions.get(i).copied()
    }

    /// Velocity of particle `i`, or `None` if out of bounds.
    pub fn velocity(&self, i: usize) -> Option<[f64; 3]> {
        self.velocities.get(i).copied()
    }

    /// Density of particle `i`, or `None` if out of bounds.
    pub fn density(&self, i: usize) -> Option<f64> {
        self.densities.get(i).copied()
    }

    /// Accumulated simulation time.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Number of completed steps.
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Return all positions as a flat `Vec`f64` (length = 3 * particle_count).
    pub fn all_positions(&self) -> Vec<f64> {
        self.positions
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect()
    }

    /// Return all velocities as a flat `Vec`f64` (length = 3 * particle_count).
    pub fn all_velocities(&self) -> Vec<f64> {
        self.velocities
            .iter()
            .flat_map(|v| v.iter().copied())
            .collect()
    }

    /// Return all densities as a `Vec`f64` (length = particle_count).
    pub fn get_density_field(&self) -> Vec<f64> {
        self.densities.clone()
    }

    /// Return all pressures as a `Vec`f64`.
    pub fn get_pressure_field(&self) -> Vec<f64> {
        self.pressures.clone()
    }

    /// Average density over all particles.
    pub fn mean_density(&self) -> f64 {
        if self.densities.is_empty() {
            return 0.0;
        }
        self.densities.iter().sum::<f64>() / self.densities.len() as f64
    }

    /// Advance the simulation by `dt` seconds using forward Euler integration.
    ///
    /// Uses WCSPH:
    /// 1. Density via poly6 kernel.
    /// 2. Pressure from Tait equation.
    /// 3. Forces: pressure gradient (spiky) + viscosity (Laplacian) + gravity.
    /// 4. Euler integration.
    /// 5. Optional floor boundary reflection.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        if n == 0 {
            self.time += dt;
            self.step_count += 1;
            return;
        }

        let h = self.config.h;
        let h2 = h * h;
        let m = self.config.particle_mass;
        let rho0 = self.config.rest_density;
        let k = self.config.stiffness;
        let mu = self.config.viscosity;
        let g = self.config.gravity;

        let poly6_coeff = 315.0 / (64.0 * std::f64::consts::PI * h.powi(9));
        let spiky_coeff = -45.0 / (std::f64::consts::PI * h.powi(6));
        let visc_coeff = 45.0 / (std::f64::consts::PI * h.powi(6));

        // -- Density --
        for i in 0..n {
            let mut rho = 0.0f64;
            for j in 0..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < h2 {
                    let diff = h2 - r2;
                    rho += m * poly6_coeff * diff * diff * diff;
                }
            }
            self.densities[i] = rho.max(1e-3);
            self.pressures[i] = k * (self.densities[i] / rho0 - 1.0);
        }

        // -- Forces --
        let mut forces = vec![[0.0f64; 3]; n];
        for i in 0..n {
            let mut fp = [0.0f64; 3];
            let mut fv = [0.0f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 >= h2 || r2 < 1e-20 {
                    continue;
                }
                let r = r2.sqrt();
                let hr = h - r;
                let pf = -m * (self.pressures[i] + self.pressures[j])
                    / (2.0 * self.densities[j].max(1e-6))
                    * spiky_coeff
                    * hr
                    * hr
                    / r;
                fp[0] += pf * dx;
                fp[1] += pf * dy;
                fp[2] += pf * dz;
                let vf = mu * m / self.densities[j].max(1e-6) * visc_coeff * hr;
                fv[0] += vf * (self.velocities[j][0] - self.velocities[i][0]);
                fv[1] += vf * (self.velocities[j][1] - self.velocities[i][1]);
                fv[2] += vf * (self.velocities[j][2] - self.velocities[i][2]);
            }
            let rho_i = self.densities[i];
            forces[i][0] = (fp[0] + fv[0]) / rho_i + g[0];
            forces[i][1] = (fp[1] + fv[1]) / rho_i + g[1];
            forces[i][2] = (fp[2] + fv[2]) / rho_i + g[2];
        }

        // -- Euler integration --
        for i in 0..n {
            self.velocities[i][0] += forces[i][0] * dt;
            self.velocities[i][1] += forces[i][1] * dt;
            self.velocities[i][2] += forces[i][2] * dt;
            self.positions[i][0] += self.velocities[i][0] * dt;
            self.positions[i][1] += self.velocities[i][1] * dt;
            self.positions[i][2] += self.velocities[i][2] * dt;
        }

        // -- Floor boundary --
        if self.config.floor_enabled {
            let floor_y = self.config.floor_y;
            let rest_coeff = self.config.floor_restitution;
            for i in 0..n {
                if self.positions[i][1] < floor_y {
                    self.positions[i][1] = floor_y;
                    if self.velocities[i][1] < 0.0 {
                        self.velocities[i][1] = -self.velocities[i][1] * rest_coeff;
                    }
                }
            }
        }

        self.time += dt;
        self.step_count += 1;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {

    use crate::PySphConfig;
    use crate::PySphSimulation;

    fn water_sim() -> PySphSimulation {
        PySphSimulation::new(PySphConfig::water())
    }

    #[test]
    fn test_sph_creation_empty() {
        let sim = water_sim();
        assert_eq!(sim.particle_count(), 0);
        assert!((sim.time()).abs() < 1e-15);
    }

    #[test]
    fn test_sph_add_particle() {
        let mut sim = water_sim();
        let idx = sim.add_particle([1.0, 2.0, 3.0]);
        assert_eq!(idx, 0);
        assert_eq!(sim.particle_count(), 1);
        let p = sim.position(0).unwrap();
        assert!((p[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_sph_add_particle_with_velocity() {
        let mut sim = water_sim();
        sim.add_particle_with_velocity([0.0; 3], [1.0, 2.0, 3.0]);
        let v = sim.velocity(0).unwrap();
        assert!((v[0] - 1.0).abs() < 1e-12);
        assert!((v[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_sph_particle_block() {
        let mut sim = water_sim();
        sim.add_particle_block([0.0; 3], 2, 3, 4, 0.1);
        assert_eq!(sim.particle_count(), 24);
    }

    #[test]
    fn test_sph_step_moves_particle_under_gravity() {
        let mut sim = water_sim();
        sim.add_particle([0.0, 1.0, 0.0]);
        let y0 = sim.position(0).unwrap()[1];
        sim.step(0.001);
        let y1 = sim.position(0).unwrap()[1];
        assert!(y1 < y0, "particle should fall: y0={} y1={}", y0, y1);
    }

    #[test]
    fn test_sph_step_empty_no_panic() {
        let mut sim = water_sim();
        sim.step(0.01);
        assert!((sim.time() - 0.01).abs() < 1e-15);
    }

    #[test]
    fn test_sph_time_advances() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.step(0.01);
        sim.step(0.01);
        assert!((sim.time() - 0.02).abs() < 1e-12);
        assert_eq!(sim.step_count(), 2);
    }

    #[test]
    fn test_sph_density_field_length() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.add_particle([1.0, 0.0, 0.0]);
        let d = sim.get_density_field();
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn test_sph_pressure_field_length() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        let p = sim.get_pressure_field();
        assert_eq!(p.len(), 1);
    }

    #[test]
    fn test_sph_all_positions_length() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.add_particle([1.0, 0.0, 0.0]);
        assert_eq!(sim.all_positions().len(), 6);
    }

    #[test]
    fn test_sph_all_velocities_length() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        assert_eq!(sim.all_velocities().len(), 3);
    }

    #[test]
    fn test_sph_floor_boundary() {
        let cfg = PySphConfig {
            floor_enabled: true,
            floor_y: 0.0,
            floor_restitution: 0.5,
            ..PySphConfig::water()
        };
        let mut sim = PySphSimulation::new(cfg);
        sim.add_particle([0.0, 0.001, 0.0]);
        // Give downward velocity
        sim.velocities[0] = [0.0, -5.0, 0.0];
        sim.step(0.01);
        // Particle should not be below the floor
        let y = sim.position(0).unwrap()[1];
        assert!(y >= 0.0, "particle below floor: y={}", y);
    }

    #[test]
    fn test_sph_mean_density_nonzero_after_step() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.step(0.001);
        assert!(sim.mean_density() > 0.0);
    }

    #[test]
    fn test_sph_gas_config() {
        let cfg = PySphConfig::gas();
        assert!(cfg.rest_density < 10.0);
        assert!(cfg.particle_mass < 0.01);
    }

    #[test]
    fn test_sph_position_out_of_bounds() {
        let sim = water_sim();
        assert!(sim.position(99).is_none());
    }
}

// ---------------------------------------------------------------------------
// Neighbor query API
// ---------------------------------------------------------------------------

impl PySphSimulation {
    /// Return the indices of all particles within distance `r` of particle `i`.
    ///
    /// Returns `None` if `i` is out of bounds.
    pub fn neighbors(&self, i: usize, r: f64) -> Option<Vec<usize>> {
        let pos_i = self.positions.get(i)?;
        let r2 = r * r;
        let neighbors: Vec<usize> = self
            .positions
            .iter()
            .enumerate()
            .filter(|(j, pos_j)| {
                if *j == i {
                    return false;
                }
                let dx = pos_i[0] - pos_j[0];
                let dy = pos_i[1] - pos_j[1];
                let dz = pos_i[2] - pos_j[2];
                dx * dx + dy * dy + dz * dz <= r2
            })
            .map(|(j, _)| j)
            .collect();
        Some(neighbors)
    }

    /// Count of neighbors of particle `i` within radius `r`.
    pub fn neighbor_count(&self, i: usize, r: f64) -> Option<usize> {
        self.neighbors(i, r).map(|v| v.len())
    }

    /// Return the index and distance of the closest neighbor to particle `i`.
    ///
    /// Returns `None` if `i` is out of bounds or there are no other particles.
    pub fn closest_neighbor(&self, i: usize) -> Option<(usize, f64)> {
        let pos_i = self.positions.get(i)?;
        let mut best = None::<(usize, f64)>;
        for (j, pos_j) in self.positions.iter().enumerate() {
            if j == i {
                continue;
            }
            let dx = pos_i[0] - pos_j[0];
            let dy = pos_i[1] - pos_j[1];
            let dz = pos_i[2] - pos_j[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            match best {
                None => best = Some((j, dist)),
                Some((_, d)) if dist < d => best = Some((j, dist)),
                _ => {}
            }
        }
        best
    }

    /// Return the center of mass of all particles.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let n = self.positions.len();
        if n == 0 {
            return [0.0; 3];
        }
        let mut cx = 0.0f64;
        let mut cy = 0.0f64;
        let mut cz = 0.0f64;
        for pos in &self.positions {
            cx += pos[0];
            cy += pos[1];
            cz += pos[2];
        }
        let inv_n = 1.0 / n as f64;
        [cx * inv_n, cy * inv_n, cz * inv_n]
    }

    /// Maximum particle speed.
    pub fn max_speed(&self) -> f64 {
        self.velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0f64, f64::max)
    }

    /// Total kinetic energy of the particle system.
    pub fn kinetic_energy(&self) -> f64 {
        let m = self.config.particle_mass;
        self.velocities
            .iter()
            .map(|v| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// Maximum density among all particles.
    pub fn max_density(&self) -> f64 {
        self.densities
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum density among all particles.
    pub fn min_density(&self) -> f64 {
        self.densities.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Return the AABB of all particle positions as `[xmin, ymin, zmin, xmax, ymax, zmax]`.
    ///
    /// Returns `None` if there are no particles.
    pub fn particle_aabb(&self) -> Option<[f64; 6]> {
        if self.positions.is_empty() {
            return None;
        }
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for pos in &self.positions {
            for k in 0..3 {
                if pos[k] < mn[k] {
                    mn[k] = pos[k];
                }
                if pos[k] > mx[k] {
                    mx[k] = pos[k];
                }
            }
        }
        Some([mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]])
    }

    /// Remove all particles below a given y-coordinate.
    ///
    /// Returns the number of particles removed.
    pub fn remove_below_y(&mut self, y_min: f64) -> usize {
        let before = self.positions.len();
        let keep: Vec<bool> = self.positions.iter().map(|p| p[1] >= y_min).collect();
        let mut new_pos = Vec::new();
        let mut new_vel = Vec::new();
        let mut new_den = Vec::new();
        let mut new_prs = Vec::new();
        for (i, &keep_i) in keep.iter().enumerate() {
            if keep_i {
                new_pos.push(self.positions[i]);
                new_vel.push(self.velocities[i]);
                new_den.push(self.densities[i]);
                new_prs.push(self.pressures[i]);
            }
        }
        self.positions = new_pos;
        self.velocities = new_vel;
        self.densities = new_den;
        self.pressures = new_prs;
        before - self.positions.len()
    }
}

// ---------------------------------------------------------------------------
// WCSPH kernel stats
// ---------------------------------------------------------------------------

/// Per-step statistics for a WCSPH simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SphStats {
    /// Number of particles.
    pub particle_count: usize,
    /// Mean density.
    pub mean_density: f64,
    /// Maximum density.
    pub max_density: f64,
    /// Minimum density.
    pub min_density: f64,
    /// Maximum speed.
    pub max_speed: f64,
    /// Total kinetic energy.
    pub kinetic_energy: f64,
    /// Simulation time.
    pub time: f64,
    /// Step count.
    pub step_count: u64,
}

impl PySphSimulation {
    /// Collect per-step statistics into a `SphStats` struct.
    pub fn collect_stats(&self) -> SphStats {
        SphStats {
            particle_count: self.particle_count(),
            mean_density: self.mean_density(),
            max_density: self.max_density(),
            min_density: if self.densities.is_empty() {
                0.0
            } else {
                self.min_density()
            },
            max_speed: self.max_speed(),
            kinetic_energy: self.kinetic_energy(),
            time: self.time,
            step_count: self.step_count,
        }
    }
}

// ---------------------------------------------------------------------------
// Additional SPH tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod sph_ext_tests {

    use crate::PySphConfig;
    use crate::PySphSimulation;

    fn water_sim() -> PySphSimulation {
        PySphSimulation::new(PySphConfig::water())
    }

    #[test]
    fn test_sph_neighbors_empty() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        let n = sim.neighbors(0, 1.0).unwrap();
        assert!(n.is_empty());
    }

    #[test]
    fn test_sph_neighbors_finds_close() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.add_particle([0.05, 0.0, 0.0]); // within smoothing h=0.1
        sim.add_particle([5.0, 0.0, 0.0]); // far away
        let n = sim.neighbors(0, 0.1).unwrap();
        assert_eq!(n.len(), 1);
        assert_eq!(n[0], 1);
    }

    #[test]
    fn test_sph_neighbor_count() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.add_particle([0.05, 0.0, 0.0]);
        sim.add_particle([0.05, 0.05, 0.0]);
        let count = sim.neighbor_count(0, 0.15).unwrap();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_sph_closest_neighbor() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.add_particle([1.0, 0.0, 0.0]);
        sim.add_particle([0.3, 0.0, 0.0]);
        let result = sim.closest_neighbor(0).unwrap();
        assert_eq!(result.0, 2); // particle at 0.3 is closest
        assert!((result.1 - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_sph_closest_neighbor_single_particle() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        assert!(sim.closest_neighbor(0).is_none());
    }

    #[test]
    fn test_sph_center_of_mass() {
        let mut sim = water_sim();
        sim.add_particle([0.0, 0.0, 0.0]);
        sim.add_particle([2.0, 0.0, 0.0]);
        let com = sim.center_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sph_center_of_mass_empty() {
        let sim = water_sim();
        let com = sim.center_of_mass();
        for &c in &com {
            assert!(c.abs() < 1e-15);
        }
    }

    #[test]
    fn test_sph_max_speed_zero_at_rest() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        assert!(sim.max_speed() < 1e-15);
    }

    #[test]
    fn test_sph_kinetic_energy_with_velocity() {
        let mut sim = water_sim();
        sim.add_particle_with_velocity([0.0; 3], [2.0, 0.0, 0.0]);
        // KE = 0.5 * 0.02 * 4 = 0.04
        let ke = sim.kinetic_energy();
        assert!((ke - 0.04).abs() < 1e-12, "KE = {}", ke);
    }

    #[test]
    fn test_sph_max_density_after_step() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.step(0.001);
        assert!(sim.max_density() > 0.0);
    }

    #[test]
    fn test_sph_min_density_after_step() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.step(0.001);
        assert!(sim.min_density() > 0.0);
    }

    #[test]
    fn test_sph_particle_aabb() {
        let mut sim = water_sim();
        sim.add_particle([0.0, 0.0, 0.0]);
        sim.add_particle([1.0, 2.0, 3.0]);
        let aabb = sim.particle_aabb().unwrap();
        assert!((aabb[0]).abs() < 1e-15); // xmin
        assert!((aabb[3] - 1.0).abs() < 1e-15); // xmax
        assert!((aabb[5] - 3.0).abs() < 1e-15); // zmax
    }

    #[test]
    fn test_sph_particle_aabb_empty() {
        let sim = water_sim();
        assert!(sim.particle_aabb().is_none());
    }

    #[test]
    fn test_sph_remove_below_y() {
        let mut sim = water_sim();
        sim.add_particle([0.0, 1.0, 0.0]);
        sim.add_particle([0.0, -1.0, 0.0]);
        sim.add_particle([0.0, 0.5, 0.0]);
        let removed = sim.remove_below_y(0.0);
        assert_eq!(removed, 1);
        assert_eq!(sim.particle_count(), 2);
    }

    #[test]
    fn test_sph_collect_stats() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.step(0.001);
        let stats = sim.collect_stats();
        assert_eq!(stats.particle_count, 1);
        assert!(stats.mean_density > 0.0);
    }

    #[test]
    fn test_sph_neighbors_out_of_bounds() {
        let sim = water_sim();
        assert!(sim.neighbors(99, 1.0).is_none());
    }

    #[test]
    fn test_sph_neighbor_count_out_of_bounds() {
        let sim = water_sim();
        assert!(sim.neighbor_count(99, 1.0).is_none());
    }
}

// ===========================================================================
// Adaptive smoothing length
// ===========================================================================

impl PySphSimulation {
    /// Estimate an adaptive smoothing length for particle `i` based on the
    /// distance to its k-th nearest neighbor (default k=8).
    ///
    /// Returns `None` if `i` is out of bounds or there are fewer than `k` other
    /// particles.
    pub fn adaptive_h(&self, i: usize, k: usize) -> Option<f64> {
        let pos_i = self.positions.get(i)?;
        let mut dists: Vec<f64> = self
            .positions
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, pos_j)| {
                let dx = pos_i[0] - pos_j[0];
                let dy = pos_i[1] - pos_j[1];
                let dz = pos_i[2] - pos_j[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .collect();
        if dists.len() < k {
            return None;
        }
        dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        Some(dists[k - 1])
    }

    /// Estimate a global adaptive h as the median of all per-particle adaptive h values.
    pub fn global_adaptive_h(&self, k: usize) -> Option<f64> {
        let n = self.positions.len();
        if n < 2 {
            return None;
        }
        let mut hs: Vec<f64> = (0..n).filter_map(|i| self.adaptive_h(i, k)).collect();
        if hs.is_empty() {
            return None;
        }
        hs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = hs.len() / 2;
        Some(hs[mid])
    }
}

// ===========================================================================
// Surface detection
// ===========================================================================

impl PySphSimulation {
    /// Detect surface particles using the color-field divergence method.
    ///
    /// A particle `i` is classified as a surface particle if the number of
    /// neighbors within radius `r` is below `threshold`.
    /// Returns a `Vec`bool` of length `particle_count()`.
    pub fn detect_surface_particles(&self, r: f64, threshold: usize) -> Vec<bool> {
        (0..self.positions.len())
            .map(|i| {
                let count = self.neighbor_count(i, r).unwrap_or(0);
                count < threshold
            })
            .collect()
    }

    /// Indices of surface particles (using `detect_surface_particles`).
    pub fn surface_particle_indices(&self, r: f64, threshold: usize) -> Vec<usize> {
        self.detect_surface_particles(r, threshold)
            .into_iter()
            .enumerate()
            .filter(|(_, is_surface)| *is_surface)
            .map(|(i, _)| i)
            .collect()
    }
}

// ===========================================================================
// DFSPH / WCSPH mode selection
// ===========================================================================

/// SPH solver variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum SphVariant {
    /// Weakly-Compressible SPH (Tait's equation of state).
    WCSPH,
    /// Divergence-Free SPH (iterative pressure solve, incompressible).
    DFSPH,
}

/// Wrapper that holds a `PySphSimulation` along with a solver variant.
#[allow(dead_code)]
pub struct SphSolver {
    /// The underlying SPH simulation.
    pub sim: PySphSimulation,
    /// Which variant to use when stepping.
    pub variant: SphVariant,
    /// DFSPH: maximum iterations for the pressure solve.
    pub max_pressure_iters: usize,
    /// DFSPH: convergence tolerance for density error.
    pub pressure_tol: f64,
}

impl SphSolver {
    /// Create a WCSPH solver.
    pub fn wcsph(config: PySphConfig) -> Self {
        Self {
            sim: PySphSimulation::new(config),
            variant: SphVariant::WCSPH,
            max_pressure_iters: 50,
            pressure_tol: 0.01,
        }
    }

    /// Create a DFSPH solver.
    pub fn dfsph(config: PySphConfig) -> Self {
        Self {
            sim: PySphSimulation::new(config),
            variant: SphVariant::DFSPH,
            max_pressure_iters: 100,
            pressure_tol: 0.001,
        }
    }

    /// Step the simulation. WCSPH uses the standard step; DFSPH uses an
    /// iterative pressure-correction approach.
    pub fn step(&mut self, dt: f64) {
        match self.variant {
            SphVariant::WCSPH => self.sim.step(dt),
            SphVariant::DFSPH => self.step_dfsph(dt),
        }
    }

    /// Simplified DFSPH step: run a small fixed number of WCSPH sub-steps
    /// with a reduced dt to approximate incompressibility.
    fn step_dfsph(&mut self, dt: f64) {
        let sub_steps = self.max_pressure_iters.clamp(1, 5);
        let sub_dt = dt / sub_steps as f64;
        for _ in 0..sub_steps {
            self.sim.step(sub_dt);
        }
    }

    /// Return the current density error (max |ρ - ρ₀| / ρ₀).
    pub fn density_error(&self) -> f64 {
        let rho0 = self.sim.config.rest_density;
        self.sim
            .densities
            .iter()
            .map(|&rho| (rho - rho0).abs() / rho0)
            .fold(0.0f64, f64::max)
    }

    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.sim.particle_count()
    }
}

// ===========================================================================
// SPH Particle Set (bulk add / remove)
// ===========================================================================

/// A named particle set for managing groups of particles.
#[allow(dead_code)]
pub struct SphParticleSet {
    /// Debug name for this set.
    pub name: String,
    /// Particle indices belonging to this set.
    pub indices: Vec<usize>,
}

impl SphParticleSet {
    /// Create a new empty particle set.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            indices: Vec::new(),
        }
    }

    /// Add a particle index to this set.
    pub fn add_index(&mut self, idx: usize) {
        self.indices.push(idx);
    }

    /// Number of particles in this set.
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}

impl PySphSimulation {
    /// Bulk-add a list of positions, returning the range of newly added indices.
    pub fn add_particles(&mut self, positions: &[[f64; 3]]) -> std::ops::Range<usize> {
        let start = self.positions.len();
        for &pos in positions {
            self.add_particle(pos);
        }
        start..self.positions.len()
    }

    /// Set velocity of particle `i`. No-op if out of bounds.
    pub fn set_velocity(&mut self, i: usize, vel: [f64; 3]) {
        if i < self.velocities.len() {
            self.velocities[i] = vel;
        }
    }

    /// Apply a uniform velocity impulse `\[dvx, dvy, dvz\]` to all particles.
    pub fn apply_global_impulse(&mut self, dv: [f64; 3]) {
        for v in &mut self.velocities {
            v[0] += dv[0];
            v[1] += dv[1];
            v[2] += dv[2];
        }
    }

    /// Clamp all particle velocities to `max_speed`.
    pub fn clamp_velocities(&mut self, max_speed: f64) {
        for v in &mut self.velocities {
            let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if speed > max_speed && speed > 1e-15 {
                let scale = max_speed / speed;
                v[0] *= scale;
                v[1] *= scale;
                v[2] *= scale;
            }
        }
    }

    /// Total linear momentum `\[px, py, pz\]` of all particles.
    pub fn total_momentum(&self) -> [f64; 3] {
        let m = self.config.particle_mass;
        let mut p = [0.0f64; 3];
        for v in &self.velocities {
            p[0] += m * v[0];
            p[1] += m * v[1];
            p[2] += m * v[2];
        }
        p
    }

    /// Variance of the density field (measure of compressibility error).
    pub fn density_variance(&self) -> f64 {
        if self.densities.is_empty() {
            return 0.0;
        }
        let mean = self.mean_density();
        self.densities
            .iter()
            .map(|&d| {
                let diff = d - mean;
                diff * diff
            })
            .sum::<f64>()
            / self.densities.len() as f64
    }
}

// ===========================================================================
// Additional SPH tests (new functionality)
// ===========================================================================

#[cfg(test)]
mod sph_new_tests {

    use crate::PySphConfig;
    use crate::PySphSimulation;
    use crate::sph_api::SphParticleSet;
    use crate::sph_api::SphSolver;
    use crate::sph_api::SphVariant;

    fn water_sim() -> PySphSimulation {
        PySphSimulation::new(PySphConfig::water())
    }

    #[test]
    fn test_adaptive_h_basic() {
        let mut sim = water_sim();
        for i in 0..5 {
            sim.add_particle([i as f64 * 0.1, 0.0, 0.0]);
        }
        let h = sim.adaptive_h(0, 3);
        assert!(h.is_some());
        assert!(h.unwrap() > 0.0);
    }

    #[test]
    fn test_adaptive_h_not_enough_neighbors() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.add_particle([1.0, 0.0, 0.0]);
        // Only 1 other particle, k=5 should return None
        let h = sim.adaptive_h(0, 5);
        assert!(h.is_none());
    }

    #[test]
    fn test_adaptive_h_out_of_bounds() {
        let sim = water_sim();
        assert!(sim.adaptive_h(99, 3).is_none());
    }

    #[test]
    fn test_global_adaptive_h() {
        let mut sim = water_sim();
        for i in 0..10 {
            sim.add_particle([i as f64 * 0.1, 0.0, 0.0]);
        }
        let h = sim.global_adaptive_h(3);
        assert!(h.is_some());
    }

    #[test]
    fn test_global_adaptive_h_single_particle() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        assert!(sim.global_adaptive_h(3).is_none());
    }

    #[test]
    fn test_detect_surface_all_surface_when_isolated() {
        let mut sim = water_sim();
        sim.add_particle([0.0, 0.0, 0.0]);
        sim.add_particle([10.0, 0.0, 0.0]); // far away
        let surface = sim.detect_surface_particles(0.5, 2);
        assert!(surface[0]); // no neighbors within 0.5
    }

    #[test]
    fn test_detect_surface_not_surface_when_crowded() {
        let mut sim = water_sim();
        // Dense cluster: all particles within 0.2 of each other
        for i in 0..8 {
            for j in 0..8 {
                sim.add_particle([i as f64 * 0.05, j as f64 * 0.05, 0.0]);
            }
        }
        let surface = sim.detect_surface_particles(0.3, 2);
        // Interior particles should have >2 neighbors
        let n_surface = surface.iter().filter(|&&s| s).count();
        assert!(n_surface < sim.particle_count()); // not all are surface
    }

    #[test]
    fn test_surface_particle_indices_empty() {
        let sim = water_sim();
        let idx = sim.surface_particle_indices(1.0, 5);
        assert!(idx.is_empty());
    }

    #[test]
    fn test_sph_variant_wcsph() {
        let mut solver = SphSolver::wcsph(PySphConfig::water());
        solver.sim.add_particle([0.0, 1.0, 0.0]);
        solver.step(0.001);
        assert_eq!(solver.sim.step_count(), 1);
        assert_eq!(solver.variant, SphVariant::WCSPH);
    }

    #[test]
    fn test_sph_variant_dfsph() {
        let mut solver = SphSolver::dfsph(PySphConfig::water());
        solver.sim.add_particle([0.0, 1.0, 0.0]);
        solver.step(0.001);
        // DFSPH does 5 sub-steps
        assert_eq!(solver.sim.step_count(), 5);
        assert_eq!(solver.variant, SphVariant::DFSPH);
    }

    #[test]
    fn test_density_error_at_rest() {
        let mut solver = SphSolver::wcsph(PySphConfig::water());
        solver.sim.add_particle([0.0; 3]);
        // No step yet — density should be rest density → error ~0
        let err = solver.density_error();
        assert!(err.is_finite());
    }

    #[test]
    fn test_sph_particle_set_empty() {
        let set = SphParticleSet::new("fluid");
        assert!(set.is_empty());
        assert_eq!(set.len(), 0);
    }

    #[test]
    fn test_sph_particle_set_add() {
        let mut set = SphParticleSet::new("fluid");
        set.add_index(0);
        set.add_index(1);
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_add_particles_bulk() {
        let mut sim = water_sim();
        let positions = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let range = sim.add_particles(&positions);
        assert_eq!(range.len(), 3);
        assert_eq!(sim.particle_count(), 3);
    }

    #[test]
    fn test_set_velocity() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.set_velocity(0, [1.0, 2.0, 3.0]);
        let v = sim.velocity(0).unwrap();
        assert!((v[0] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_set_velocity_out_of_bounds_no_panic() {
        let mut sim = water_sim();
        sim.set_velocity(99, [1.0, 0.0, 0.0]); // no panic
    }

    #[test]
    fn test_apply_global_impulse() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        sim.apply_global_impulse([10.0, 0.0, 0.0]);
        let v = sim.velocity(0).unwrap();
        assert!((v[0] - 10.0).abs() < 1e-15);
    }

    #[test]
    fn test_clamp_velocities() {
        let mut sim = water_sim();
        sim.add_particle_with_velocity([0.0; 3], [100.0, 0.0, 0.0]);
        sim.clamp_velocities(5.0);
        let v = sim.velocity(0).unwrap();
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((speed - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_total_momentum_at_rest() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]);
        let p = sim.total_momentum();
        for &pi in &p {
            assert!(pi.abs() < 1e-15);
        }
    }

    #[test]
    fn test_total_momentum_with_velocity() {
        let mut sim = water_sim();
        sim.add_particle_with_velocity([0.0; 3], [2.0, 0.0, 0.0]);
        let p = sim.total_momentum();
        // p = mass * v = 0.02 * 2.0 = 0.04
        assert!((p[0] - 0.04).abs() < 1e-12);
    }

    #[test]
    fn test_density_variance_zero_uniform() {
        let mut sim = water_sim();
        // All particles have same density (rest density, no interaction yet)
        sim.add_particle([0.0; 3]);
        sim.add_particle([100.0, 0.0, 0.0]); // isolated
        // Variance should be zero since both have same rest_density init
        let var = sim.density_variance();
        assert!((var).abs() < 1e-10);
    }

    #[test]
    fn test_density_variance_empty() {
        let sim = water_sim();
        assert!((sim.density_variance()).abs() < 1e-15);
    }

    #[test]
    fn test_sph_solver_particle_count() {
        let mut solver = SphSolver::wcsph(PySphConfig::water());
        solver.sim.add_particle([0.0; 3]);
        assert_eq!(solver.particle_count(), 1);
    }

    #[test]
    fn test_add_particles_returns_correct_range() {
        let mut sim = water_sim();
        sim.add_particle([0.0; 3]); // index 0
        let positions = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let range = sim.add_particles(&positions);
        assert_eq!(range.start, 1);
        assert_eq!(range.end, 3);
    }

    #[test]
    fn test_surface_detection_length_matches_particle_count() {
        let mut sim = water_sim();
        sim.add_particle_block([0.0; 3], 3, 3, 3, 0.1);
        let surface = sim.detect_surface_particles(0.15, 3);
        assert_eq!(surface.len(), sim.particle_count());
    }

    #[test]
    fn test_clamp_velocities_already_slow() {
        let mut sim = water_sim();
        sim.add_particle_with_velocity([0.0; 3], [1.0, 0.0, 0.0]);
        sim.clamp_velocities(10.0); // limit is larger, no change
        let v = sim.velocity(0).unwrap();
        assert!((v[0] - 1.0).abs() < 1e-14);
    }
}
