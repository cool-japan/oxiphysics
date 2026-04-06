#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH (Smoothed Particle Hydrodynamics) Simulation.

#![allow(missing_docs)]

// ===========================================================================
// SPH (Smoothed Particle Hydrodynamics) Simulation
// ===========================================================================

/// Configuration for an SPH particle simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PySphConfig {
    /// Smoothing length (kernel radius).
    pub h: f64,
    /// Reference density (kg/m^3).
    pub rest_density: f64,
    /// Pressure stiffness constant k.
    pub stiffness: f64,
    /// Dynamic viscosity coefficient.
    pub viscosity: f64,
    /// Gravity `[gx, gy, gz]`.
    pub gravity: [f64; 3],
    /// Particle mass.
    pub particle_mass: f64,
}

impl PySphConfig {
    /// Create a default water-like SPH configuration.
    pub fn water() -> Self {
        Self {
            h: 0.1,
            rest_density: 1000.0,
            stiffness: 200.0,
            viscosity: 0.01,
            gravity: [0.0, -9.81, 0.0],
            particle_mass: 0.02,
        }
    }
}

impl Default for PySphConfig {
    fn default() -> Self {
        Self::water()
    }
}

/// SPH particle simulation (simplified 3-D, poly6/spiky kernels).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PySphSim {
    /// Particle positions `[x, y, z]` per particle.
    positions: Vec<[f64; 3]>,
    /// Particle velocities `[vx, vy, vz]` per particle.
    velocities: Vec<[f64; 3]>,
    /// Particle densities (one per particle).
    densities: Vec<f64>,
    /// Particle pressures.
    pressures: Vec<f64>,
    /// Configuration.
    config: PySphConfig,
    /// Accumulated simulation time.
    time: f64,
}

impl PySphSim {
    /// Create an empty SPH simulation.
    pub fn new(config: PySphConfig) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            config,
            time: 0.0,
        }
    }

    /// Add a particle at position `[x, y, z]` with zero velocity.
    pub fn add_particle(&mut self, pos: [f64; 3]) {
        self.positions.push(pos);
        self.velocities.push([0.0; 3]);
        self.densities.push(self.config.rest_density);
        self.pressures.push(0.0);
    }

    /// Add `n` particles arranged in a 3-D lattice starting at `origin` with
    /// spacing `dx`.
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

    /// Get the position of particle `i`, or `None` if out of bounds.
    pub fn position(&self, i: usize) -> Option<[f64; 3]> {
        self.positions.get(i).copied()
    }

    /// Get the velocity of particle `i`, or `None` if out of bounds.
    pub fn velocity(&self, i: usize) -> Option<[f64; 3]> {
        self.velocities.get(i).copied()
    }

    /// Get the density of particle `i`, or `None` if out of bounds.
    pub fn density(&self, i: usize) -> Option<f64> {
        self.densities.get(i).copied()
    }

    /// Get the pressure of particle `i`, or `None` if out of bounds.
    pub fn pressure(&self, i: usize) -> Option<f64> {
        self.pressures.get(i).copied()
    }

    /// Mutable access to the velocity slice (for initial-condition setup).
    pub fn velocities_mut(&mut self) -> &mut Vec<[f64; 3]> {
        &mut self.velocities
    }

    /// Smoothing length `h`.
    pub fn smoothing_length(&self) -> f64 {
        self.config.h
    }

    /// Rest density ρ₀.
    pub fn rest_density(&self) -> f64 {
        self.config.rest_density
    }

    /// Particle mass.
    pub fn particle_mass(&self) -> f64 {
        self.config.particle_mass
    }

    /// Return all positions as a flat `Vec`f64` of `\[x, y, z\]` triples.
    pub fn all_positions(&self) -> Vec<f64> {
        self.positions
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect()
    }

    /// Return all velocities as a flat `Vec`f64` of `[vx, vy, vz]` triples.
    pub fn all_velocities(&self) -> Vec<f64> {
        self.velocities
            .iter()
            .flat_map(|v| v.iter().copied())
            .collect()
    }

    /// Accumulated simulation time in seconds.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Advance the simulation by `dt` seconds.
    ///
    /// Uses a simple WCSPH scheme:
    /// 1. Compute densities and pressures (poly6 kernel).
    /// 2. Compute pressure + viscosity forces (spiky/laplacian kernels).
    /// 3. Integrate with leapfrog.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        if n == 0 {
            return;
        }

        let h = self.config.h;
        let h2 = h * h;
        let m = self.config.particle_mass;
        let rho0 = self.config.rest_density;
        let k = self.config.stiffness;
        let mu = self.config.viscosity;
        let g = self.config.gravity;

        // -- Density computation (poly6 kernel) --
        let poly6_coeff = 315.0 / (64.0 * std::f64::consts::PI * h.powi(9));
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
            self.pressures[i] = k * (self.densities[i] - rho0);
        }

        // -- Force computation --
        let spiky_coeff = -45.0 / (std::f64::consts::PI * h.powi(6));
        let visc_coeff = 45.0 / (std::f64::consts::PI * h.powi(6));
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
                // Pressure force (spiky gradient)
                let pf = -m * (self.pressures[i] + self.pressures[j]) / (2.0 * self.densities[j])
                    * spiky_coeff
                    * hr
                    * hr
                    / r;
                fp[0] += pf * dx;
                fp[1] += pf * dy;
                fp[2] += pf * dz;
                // Viscosity force (laplacian)
                let vf = mu * m / self.densities[j] * visc_coeff * hr;
                fv[0] += vf * (self.velocities[j][0] - self.velocities[i][0]);
                fv[1] += vf * (self.velocities[j][1] - self.velocities[i][1]);
                fv[2] += vf * (self.velocities[j][2] - self.velocities[i][2]);
            }
            // Gravity
            let rho_i = self.densities[i];
            forces[i][0] = (fp[0] + fv[0]) / rho_i + g[0];
            forces[i][1] = (fp[1] + fv[1]) / rho_i + g[1];
            forces[i][2] = (fp[2] + fv[2]) / rho_i + g[2];
        }

        // -- Integration (Euler) --
        for i in 0..n {
            self.velocities[i][0] += forces[i][0] * dt;
            self.velocities[i][1] += forces[i][1] * dt;
            self.velocities[i][2] += forces[i][2] * dt;
            self.positions[i][0] += self.velocities[i][0] * dt;
            self.positions[i][1] += self.velocities[i][1] * dt;
            self.positions[i][2] += self.velocities[i][2] * dt;
        }

        self.time += dt;
    }
}
