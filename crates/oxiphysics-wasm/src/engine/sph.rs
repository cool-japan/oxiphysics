// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmSphSim` — Smoothed Particle Hydrodynamics for the WASM boundary.

use wasm_bindgen::prelude::*;

use super::Float64View;

/// Compact SPH particle data for the WASM boundary.
#[wasm_bindgen]
#[derive(Debug)]
pub struct WasmSphSim {
    positions: Vec<[f64; 3]>,
    velocities: Vec<[f64; 3]>,
    densities: Vec<f64>,
    /// Smoothing length.
    h: f64,
    /// Rest density.
    rho0: f64,
    /// Stiffness constant.
    k: f64,
    /// Viscosity coefficient.
    mu: f64,
    /// Particle mass.
    mass: f64,
    /// Gravity.
    gravity: [f64; 3],
    /// Accumulated time.
    time: f64,
}

#[wasm_bindgen]
impl WasmSphSim {
    /// Create a new SPH simulation with default water-like parameters.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            h: 0.1,
            rho0: 1000.0,
            k: 200.0,
            mu: 0.01,
            mass: 0.02,
            gravity: [0.0, -9.81, 0.0],
            time: 0.0,
        }
    }

    /// Configure the simulation parameters.
    pub fn configure(&mut self, h: f64, rho0: f64, stiffness: f64, viscosity: f64, mass: f64) {
        self.h = h.max(1e-6);
        self.rho0 = rho0;
        self.k = stiffness;
        self.mu = viscosity;
        self.mass = mass;
    }

    /// Set gravity vector.
    pub fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
    }

    /// Add a particle at position `(x, y, z)`.
    pub fn add_particle(&mut self, x: f64, y: f64, z: f64) {
        self.positions.push([x, y, z]);
        self.velocities.push([0.0; 3]);
        self.densities.push(self.rho0);
    }

    /// Number of particles.
    pub fn particle_count(&self) -> u32 {
        self.positions.len() as u32
    }

    /// Get position of particle `i` as `Vec<f64>` of `[x, y, z]` (JS-compatible).
    pub fn get_position_js(&self, i: u32) -> Vec<f64> {
        self.positions
            .get(i as usize)
            .map(|p| p.to_vec())
            .unwrap_or_else(|| vec![0.0; 3])
    }

    /// Get velocity of particle `i` as `Vec<f64>` of `[vx, vy, vz]` (JS-compatible).
    pub fn get_velocity_js(&self, i: u32) -> Vec<f64> {
        self.velocities
            .get(i as usize)
            .map(|v| v.to_vec())
            .unwrap_or_else(|| vec![0.0; 3])
    }

    /// Return all positions as a flat `Float64View` of `[x, y, z]` triples.
    pub fn all_positions(&self) -> Float64View {
        Float64View::from_vec(
            self.positions
                .iter()
                .flat_map(|p| p.iter().copied())
                .collect(),
        )
    }

    /// Return all velocities as a flat `Float64View`.
    pub fn all_velocities(&self) -> Float64View {
        Float64View::from_vec(
            self.velocities
                .iter()
                .flat_map(|v| v.iter().copied())
                .collect(),
        )
    }

    /// Return all densities as a `Float64View`.
    pub fn all_densities(&self) -> Float64View {
        Float64View::from_vec(self.densities.clone())
    }

    /// Return all positions as a `js_sys::Float64Array` (JS-typed-array, zero-copy-friendly).
    pub fn all_positions_typed(&self) -> js_sys::Float64Array {
        self.all_positions().to_typed_array()
    }

    /// Return all velocities as a `js_sys::Float64Array`.
    pub fn all_velocities_typed(&self) -> js_sys::Float64Array {
        self.all_velocities().to_typed_array()
    }

    /// Return all densities as a `js_sys::Float64Array`.
    pub fn all_densities_typed(&self) -> js_sys::Float64Array {
        self.all_densities().to_typed_array()
    }

    /// Accumulated simulation time.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Advance the simulation by `dt` seconds using a simple WCSPH scheme.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        if n == 0 {
            return;
        }

        let h2 = self.h * self.h;
        let poly6 = 315.0 / (64.0 * std::f64::consts::PI * self.h.powi(9));
        let spiky = -45.0 / (std::f64::consts::PI * self.h.powi(6));
        let visc_k = 45.0 / (std::f64::consts::PI * self.h.powi(6));

        // Density
        let mut rho = vec![0.0f64; n];
        for (i, rho_i) in rho.iter_mut().enumerate().take(n) {
            for j in 0..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < h2 {
                    let d = h2 - r2;
                    *rho_i += self.mass * poly6 * d * d * d;
                }
            }
            self.densities[i] = rho_i.max(1e-3);
        }

        // Pressure = k * (rho - rho0)
        let pressure: Vec<f64> = self
            .densities
            .iter()
            .map(|&r| self.k * (r - self.rho0))
            .collect();

        // Forces
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
                let hr = self.h - r;
                let pf = -self.mass * (pressure[i] + pressure[j]) / (2.0 * self.densities[j])
                    * spiky
                    * hr
                    * hr
                    / r;
                fp[0] += pf * dx;
                fp[1] += pf * dy;
                fp[2] += pf * dz;
                let vf = self.mu * self.mass / self.densities[j] * visc_k * hr;
                fv[0] += vf * (self.velocities[j][0] - self.velocities[i][0]);
                fv[1] += vf * (self.velocities[j][1] - self.velocities[i][1]);
                fv[2] += vf * (self.velocities[j][2] - self.velocities[i][2]);
            }
            let ri = self.densities[i];
            forces[i][0] = (fp[0] + fv[0]) / ri + self.gravity[0];
            forces[i][1] = (fp[1] + fv[1]) / ri + self.gravity[1];
            forces[i][2] = (fp[2] + fv[2]) / ri + self.gravity[2];
        }

        // Euler integration
        for (i, force) in forces.iter().enumerate().take(n) {
            for (k, &f) in force.iter().enumerate().take(3) {
                self.velocities[i][k] += f * dt;
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
        }

        self.time += dt;
    }
}

impl WasmSphSim {
    /// Get position of particle `i` as `[x, y, z]` (Rust-only).
    pub fn get_position(&self, i: u32) -> [f64; 3] {
        self.positions.get(i as usize).copied().unwrap_or([0.0; 3])
    }

    /// Get velocity of particle `i` as `[vx, vy, vz]` (Rust-only).
    pub fn get_velocity(&self, i: u32) -> [f64; 3] {
        self.velocities.get(i as usize).copied().unwrap_or([0.0; 3])
    }
}

impl Default for WasmSphSim {
    fn default() -> Self {
        Self::new()
    }
}
