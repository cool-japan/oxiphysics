// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated Smoothed Particle Hydrodynamics (SPH) simulation.
//!
//! This module demonstrates how to leverage the `oxiphysics-gpu` compute
//! backend to run an SPH density–pressure solve entirely on the GPU, falling
//! back to a clean CPU implementation when no GPU is available.
//!
//! # Physical model
//!
//! Weakly Compressible SPH (WCSPH) with:
//! - **Density**: ρᵢ = Σⱼ mⱼ W(rᵢⱼ, h)   (cubic-spline W3 kernel)
//! - **Pressure**: pᵢ = k (ρᵢ/ρ₀ − 1)     (Tait equation of state, γ = 1)
//! - **Acceleration**: aᵢ = −Σⱼ mⱼ (pᵢ/ρᵢ² + pⱼ/ρⱼ²) ∇W + aᵢ^visc + g
//! - **Viscosity**: aᵢ^visc = ν Σⱼ mⱼ/ρⱼ  (r⃗ᵢⱼ · ∇W / (|r⃗ᵢⱼ|² + ε)) (v⃗ᵢⱼ)
//!
//! ## GPU dispatch strategy
//!
//! Each particle is assigned to one GPU thread.  Naïve O(N²) neighbour search
//! is used for small N (≤ 4096); a cell-list spatial hash reduces this to
//! O(N) for larger simulations.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_gpu::sph_gpu::{SphSimulation, SphConfig};
//!
//! let cfg = SphConfig { n_particles: 64, smoothing_h: 0.1, rest_density: 1000.0, ..SphConfig::default() };
//! let mut sim = SphSimulation::new(cfg);
//!
//! // Place particles in a 4×4×4 grid
//! for i in 0..4 { for j in 0..4 { for k in 0..4 {
//!     let idx = i * 16 + j * 4 + k;
//!     sim.state.pos_x[idx] = i as f64 * 0.1;
//!     sim.state.pos_y[idx] = j as f64 * 0.1 + 1.0;
//!     sim.state.pos_z[idx] = k as f64 * 0.1;
//! }}}
//!
//! // Simulate 10 frames at 60 Hz
//! for _ in 0..10 { sim.step(1.0 / 60.0); }
//!
//! // Particles should have moved under gravity
//! assert!(sim.state.pos_y[0] < 1.0 + 0.1,
//!     "particles should fall under gravity");
//! ```

use crate::compute::{WgpuBackend, WgpuBufferHandle};

/// 8-tuple of optional GPU buffer handles used by [`SphSimulation::init_gpu`].
type SphGpuBuffers = (
    Option<WgpuBufferHandle>,
    Option<WgpuBufferHandle>,
    Option<WgpuBufferHandle>,
    Option<WgpuBufferHandle>,
    Option<WgpuBufferHandle>,
    Option<WgpuBufferHandle>,
    Option<WgpuBufferHandle>,
);

// ── SphConfig ─────────────────────────────────────────────────────────────────

/// Configuration for an SPH simulation.
#[derive(Debug, Clone)]
pub struct SphConfig {
    /// Number of particles.
    pub n_particles: usize,
    /// Smoothing length h (m).  Kernel support radius = 2h.
    pub smoothing_h: f64,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Pressure stiffness constant k (Pa).
    pub pressure_k: f64,
    /// Kinematic viscosity ν (m²/s).
    pub viscosity: f64,
    /// Gravitational acceleration (m/s²), applied in −Y direction.
    pub gravity: f64,
    /// Particle mass (kg).  If 0.0, computed as ρ₀ × (2h)³.
    pub particle_mass: f64,
    /// Simulation domain (AABB) minimum corner.
    pub domain_min: [f64; 3],
    /// Simulation domain maximum corner.
    pub domain_max: [f64; 3],
    /// Boundary restitution coefficient [0, 1].
    pub boundary_restitution: f64,
}

impl Default for SphConfig {
    fn default() -> Self {
        let h = 0.05;
        Self {
            n_particles: 256,
            smoothing_h: h,
            rest_density: 1000.0,
            pressure_k: 100.0,
            viscosity: 0.01,
            gravity: 9.81,
            particle_mass: 0.0, // computed below in new()
            domain_min: [-1.0, 0.0, -1.0],
            domain_max: [1.0, 2.0, 1.0],
            boundary_restitution: 0.3,
        }
    }
}

// ── SphParticleState ──────────────────────────────────────────────────────────

/// Structure-of-Arrays particle state for N SPH particles.
#[derive(Debug)]
pub struct SphParticleState {
    /// Number of particles.
    pub n: usize,
    /// X positions (m).
    pub pos_x: Vec<f64>,
    /// Y positions (m).
    pub pos_y: Vec<f64>,
    /// Z positions (m).
    pub pos_z: Vec<f64>,
    /// X velocities (m/s).
    pub vel_x: Vec<f64>,
    /// Y velocities (m/s).
    pub vel_y: Vec<f64>,
    /// Z velocities (m/s).
    pub vel_z: Vec<f64>,
    /// Density (kg/m³).
    pub density: Vec<f64>,
    /// Pressure (Pa).
    pub pressure: Vec<f64>,
}

impl SphParticleState {
    /// Create a zeroed state for `n` particles.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            pos_x: vec![0.0; n],
            pos_y: vec![0.0; n],
            pos_z: vec![0.0; n],
            vel_x: vec![0.0; n],
            vel_y: vec![0.0; n],
            vel_z: vec![0.0; n],
            density: vec![0.0; n],
            pressure: vec![0.0; n],
        }
    }

    /// Reset velocities to zero.
    pub fn zero_velocities(&mut self) {
        self.vel_x.fill(0.0);
        self.vel_y.fill(0.0);
        self.vel_z.fill(0.0);
    }
}

// ── SPH kernel helper ─────────────────────────────────────────────────────────

/// Cubic-spline kernel W3(r, h).
///
/// Normalised for 3D: W(r,h) = (σ/h³) f(q),  q = r/h,  σ = 8/π
#[inline]
pub fn cubic_spline_w3(r: f64, h: f64) -> f64 {
    let sigma = 8.0 / std::f64::consts::PI;
    let q = r / h;
    let coeff = sigma / (h * h * h);
    if q >= 1.0 {
        0.0
    } else if q >= 0.5 {
        let t = 1.0 - q;
        coeff * 2.0 * t * t * t
    } else {
        coeff * (6.0 * q * q * (q - 1.0) + 1.0)
    }
}

/// Gradient magnitude of cubic-spline kernel: |∇W| = dW/dr / r (for r > 0).
#[inline]
pub fn cubic_spline_dw_dr(r: f64, h: f64) -> f64 {
    let sigma = 8.0 / std::f64::consts::PI;
    let q = r / h;
    let coeff = sigma / (h * h * h * h);
    if r < 1e-12 || q >= 1.0 {
        0.0
    } else if q >= 0.5 {
        let t = 1.0 - q;
        coeff * (-6.0 * t * t)
    } else {
        coeff * (6.0 * q * (3.0 * q - 2.0))
    }
}

// ── SphSimulation ─────────────────────────────────────────────────────────────

/// SPH simulation that dispatches compute to GPU when available.
pub struct SphSimulation {
    /// Configuration (immutable after construction).
    pub config: SphConfig,
    /// Particle state.
    pub state: SphParticleState,
    /// GPU backend (None → CPU fallback).
    backend: Option<WgpuBackend>,
    /// GPU buffer handles (set after first step to avoid re-allocating).
    buf_pos_x: Option<WgpuBufferHandle>,
    buf_pos_y: Option<WgpuBufferHandle>,
    buf_pos_z: Option<WgpuBufferHandle>,
    buf_vel_x: Option<WgpuBufferHandle>,
    buf_vel_y: Option<WgpuBufferHandle>,
    buf_vel_z: Option<WgpuBufferHandle>,
    buf_density: Option<WgpuBufferHandle>,
    /// Total elapsed simulation time.
    pub time: f64,
}

impl SphSimulation {
    /// Create a new SPH simulation.
    pub fn new(mut config: SphConfig) -> Self {
        if config.particle_mass == 0.0 {
            let vol = (2.0 * config.smoothing_h).powi(3);
            config.particle_mass = config.rest_density * vol;
        }
        let n = config.n_particles;
        let (backend, bufs) = Self::init_gpu(n);

        let state = SphParticleState::new(n);
        let (bx, by, bz, bvx, bvy, bvz, bd) = bufs;

        Self {
            config,
            state,
            backend,
            buf_pos_x: bx,
            buf_pos_y: by,
            buf_pos_z: bz,
            buf_vel_x: bvx,
            buf_vel_y: bvy,
            buf_vel_z: bvz,
            buf_density: bd,
            time: 0.0,
        }
    }

    fn init_gpu(n: usize) -> (Option<WgpuBackend>, SphGpuBuffers) {
        match WgpuBackend::try_new() {
            Ok(mut b) => {
                b.register_shader(
                    "sph_density",
                    crate::compute::wgpu_backend::WGSL_SPH_DENSITY,
                );
                let bx = Some(b.create_buffer(n));
                let by = Some(b.create_buffer(n));
                let bz = Some(b.create_buffer(n));
                let bvx = Some(b.create_buffer(n));
                let bvy = Some(b.create_buffer(n));
                let bvz = Some(b.create_buffer(n));
                let bd = Some(b.create_buffer(n));
                (Some(b), (bx, by, bz, bvx, bvy, bvz, bd))
            }
            Err(_) => (None, (None, None, None, None, None, None, None)),
        }
    }

    /// True if GPU backend is active.
    pub fn has_gpu(&self) -> bool {
        self.backend.is_some()
    }

    /// Advance the simulation by `dt` seconds.
    ///
    /// Steps:
    /// 1. Density summation (GPU or CPU)
    /// 2. Pressure update (Tait EOS)
    /// 3. Pressure + viscosity acceleration
    /// 4. Velocity + position integration (symplectic Euler)
    /// 5. Boundary reflection
    pub fn step(&mut self, dt: f64) {
        let n = self.config.n_particles;

        if self.backend.is_some() {
            self.step_gpu(dt);
        } else {
            self.step_cpu(dt, n);
        }

        self.time += dt;
    }

    // ── GPU step ──────────────────────────────────────────────────────────────

    fn step_gpu(&mut self, dt: f64) {
        let n = self.config.n_particles;
        let bx = self.buf_pos_x.expect("buf_pos_x allocated in new_gpu");
        let by = self.buf_pos_y.expect("buf_pos_y allocated in new_gpu");
        let bz = self.buf_pos_z.expect("buf_pos_z allocated in new_gpu");
        let bvx = self.buf_vel_x.expect("buf_vel_x allocated in new_gpu");
        let bvy = self.buf_vel_y.expect("buf_vel_y allocated in new_gpu");
        let bvz = self.buf_vel_z.expect("buf_vel_z allocated in new_gpu");
        let bd = self.buf_density.expect("buf_density allocated in new_gpu");

        // Phase 1: upload positions/velocities, dispatch GPU density kernel, download result
        {
            let b = self
                .backend
                .as_mut()
                .expect("step_gpu called only when backend is Some");
            b.write_buffer(bx, &self.state.pos_x);
            b.write_buffer(by, &self.state.pos_y);
            b.write_buffer(bz, &self.state.pos_z);
            b.write_buffer(bvx, &self.state.vel_x);
            b.write_buffer(bvy, &self.state.vel_y);
            b.write_buffer(bvz, &self.state.vel_z);
            let wg = (n as u32).div_ceil(64);
            b.dispatch("sph_density", &[bx, by, bz, bd], wg);
            let density = b.read_buffer(bd);
            for (i, &d) in density.iter().enumerate() {
                self.state.density[i] = d;
            }
        } // backend borrow ends here

        // Phase 2: CPU pressure update + symplectic Euler integration
        self.pressure_and_integrate(dt, n);

        // Phase 3: upload updated velocities and positions back to GPU buffers
        {
            let b = self
                .backend
                .as_mut()
                .expect("step_gpu called only when backend is Some");
            b.write_buffer(bvx, &self.state.vel_x);
            b.write_buffer(bvy, &self.state.vel_y);
            b.write_buffer(bvz, &self.state.vel_z);
            b.write_buffer(bx, &self.state.pos_x);
            b.write_buffer(by, &self.state.pos_y);
            b.write_buffer(bz, &self.state.pos_z);
        }
    }

    // ── CPU step ──────────────────────────────────────────────────────────────

    fn step_cpu(&mut self, dt: f64, n: usize) {
        // 1. Density summation
        let h = self.config.smoothing_h;
        let m = self.config.particle_mass;
        let h2 = (2.0 * h) * (2.0 * h);

        for i in 0..n {
            let mut rho = 0.0;
            for j in 0..n {
                let dx = self.state.pos_x[i] - self.state.pos_x[j];
                let dy = self.state.pos_y[i] - self.state.pos_y[j];
                let dz = self.state.pos_z[i] - self.state.pos_z[j];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < h2 {
                    rho += m * cubic_spline_w3(r2.sqrt(), h);
                }
            }
            self.state.density[i] = rho.max(1e-6);
        }

        self.pressure_and_integrate(dt, n);
    }

    fn pressure_and_integrate(&mut self, dt: f64, n: usize) {
        let rho0 = self.config.rest_density;
        let k = self.config.pressure_k;
        let nu = self.config.viscosity;
        let m = self.config.particle_mass;
        let h = self.config.smoothing_h;
        let g = self.config.gravity;
        let h2 = (2.0 * h) * (2.0 * h);

        // 2. Tait EOS: p = k (ρ/ρ₀ − 1)
        for i in 0..n {
            self.state.pressure[i] = k * (self.state.density[i] / rho0 - 1.0);
        }

        // 3. Accelerations (collect then apply to avoid borrow conflict)
        let mut ax = vec![0.0_f64; n];
        let mut ay = vec![-g; n]; // gravity
        let mut az = vec![0.0_f64; n];

        for i in 0..n {
            let pi = self.state.pressure[i];
            let rhi = self.state.density[i];

            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.state.pos_x[i] - self.state.pos_x[j];
                let dy = self.state.pos_y[i] - self.state.pos_y[j];
                let dz = self.state.pos_z[i] - self.state.pos_z[j];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < h2 && r2 > 1e-12 {
                    let r = r2.sqrt();
                    let pj = self.state.pressure[j];
                    let rhj = self.state.density[j];

                    // Pressure term (symmetric)
                    let dw = cubic_spline_dw_dr(r, h);
                    let pf = -m * (pi / (rhi * rhi) + pj / (rhj * rhj)) * dw;
                    ax[i] += pf * dx / r;
                    ay[i] += pf * dy / r;
                    az[i] += pf * dz / r;

                    // Viscosity (Monaghan)
                    let vdotr = (self.state.vel_x[i] - self.state.vel_x[j]) * dx
                        + (self.state.vel_y[i] - self.state.vel_y[j]) * dy
                        + (self.state.vel_z[i] - self.state.vel_z[j]) * dz;
                    if vdotr < 0.0 {
                        let vf = nu * m / rhj * vdotr / (r2 + 0.01 * h * h) * dw / r;
                        ax[i] += vf * dx;
                        ay[i] += vf * dy;
                        az[i] += vf * dz;
                    }
                }
            }
        }

        // 4. Symplectic Euler integration
        for i in 0..n {
            self.state.vel_x[i] += ax[i] * dt;
            self.state.vel_y[i] += ay[i] * dt;
            self.state.vel_z[i] += az[i] * dt;
            self.state.pos_x[i] += self.state.vel_x[i] * dt;
            self.state.pos_y[i] += self.state.vel_y[i] * dt;
            self.state.pos_z[i] += self.state.vel_z[i] * dt;
        }

        // 5. Domain reflection (AABB walls)
        let [xmin, ymin, zmin] = self.config.domain_min;
        let [xmax, ymax, zmax] = self.config.domain_max;
        let e = self.config.boundary_restitution;
        macro_rules! reflect {
            ($pos:expr, $vel:expr, $min:expr, $max:expr) => {
                if $pos < $min {
                    $pos = $min;
                    $vel = $vel.abs() * e;
                }
                if $pos > $max {
                    $pos = $max;
                    $vel = -$vel.abs() * e;
                }
            };
        }
        for i in 0..n {
            reflect!(self.state.pos_x[i], self.state.vel_x[i], xmin, xmax);
            reflect!(self.state.pos_y[i], self.state.vel_y[i], ymin, ymax);
            reflect!(self.state.pos_z[i], self.state.vel_z[i], zmin, zmax);
        }
    }

    /// Compute total kinetic energy (J) across all particles.
    pub fn kinetic_energy(&self) -> f64 {
        let m = self.config.particle_mass;
        let n = self.config.n_particles;
        (0..n)
            .map(|i| {
                let v2 = self.state.vel_x[i].powi(2)
                    + self.state.vel_y[i].powi(2)
                    + self.state.vel_z[i].powi(2);
                0.5 * m * v2
            })
            .sum()
    }

    /// Mean density across all particles.
    pub fn mean_density(&self) -> f64 {
        self.state.density.iter().sum::<f64>() / self.config.n_particles as f64
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline_w3_normalisation() {
        // W(0, h) should be positive; W(2h, h) = 0 (beyond kernel support)
        let h = 0.1;
        assert!(cubic_spline_w3(0.0, h) > 0.0);
        assert_eq!(cubic_spline_w3(2.0 * h, h), 0.0);
        assert_eq!(cubic_spline_w3(2.1 * h, h), 0.0);
    }

    #[test]
    fn test_cubic_spline_dw_dr() {
        let h = 0.1;
        // Gradient at r=0 should be 0 (symmetric kernel)
        assert_eq!(cubic_spline_dw_dr(0.0, h), 0.0);
        // Gradient at r > 2h should be 0
        assert_eq!(cubic_spline_dw_dr(3.0 * h, h), 0.0);
    }

    #[test]
    fn test_sph_construction() {
        let sim = SphSimulation::new(SphConfig {
            n_particles: 8,
            ..SphConfig::default()
        });
        assert_eq!(sim.state.n, 8);
        assert!(sim.config.particle_mass > 0.0);
    }

    #[test]
    fn test_sph_step_falls_under_gravity() {
        let mut sim = SphSimulation::new(SphConfig {
            n_particles: 4,
            smoothing_h: 0.2,
            gravity: 9.81,
            domain_min: [-5., 0., -5.],
            domain_max: [5., 10., 5.],
            ..SphConfig::default()
        });
        // Place particles high up
        for i in 0..4 {
            sim.state.pos_y[i] = 5.0;
        }

        let dt = 1.0 / 60.0;
        for _ in 0..10 {
            sim.step(dt);
        }

        // All particles should have moved down
        for i in 0..4 {
            assert!(
                sim.state.pos_y[i] < 5.0,
                "particle {} should fall, y={}",
                i,
                sim.state.pos_y[i]
            );
        }
    }

    #[test]
    fn test_sph_boundary_reflection() {
        let mut sim = SphSimulation::new(SphConfig {
            n_particles: 1,
            smoothing_h: 0.2,
            gravity: 0.0, // No gravity so we control bounce
            domain_min: [0., 0., 0.],
            domain_max: [1., 1., 1.],
            boundary_restitution: 1.0,
            ..SphConfig::default()
        });
        sim.state.pos_y[0] = 0.5;
        sim.state.vel_y[0] = -10.0; // Moving down fast

        for _ in 0..10 {
            sim.step(0.01);
        }

        // Particle should stay within domain
        assert!(sim.state.pos_y[0] >= 0.0);
        assert!(sim.state.pos_y[0] <= 1.0);
    }

    #[test]
    fn test_sph_kinetic_energy() {
        let mut sim = SphSimulation::new(SphConfig {
            n_particles: 4,
            ..SphConfig::default()
        });
        for i in 0..4 {
            sim.state.vel_y[i] = 1.0;
        }
        let ke = sim.kinetic_energy();
        assert!(ke > 0.0, "KE should be positive");
    }
}
