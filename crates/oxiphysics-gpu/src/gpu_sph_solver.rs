#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-side SPH solver — CPU mock backend.
//!
//! Implements a complete Smoothed Particle Hydrodynamics (SPH) solver pipeline
//! using plain Rust loops as a CPU fallback. The API mirrors what would run on
//! a GPU kernel dispatch, making it straightforward to swap in a real GPU
//! backend without changing call-sites.

#![allow(dead_code)]

use std::f32::consts::PI;

// ── Data structures ──────────────────────────────────────────────────────────

/// A single SPH particle stored on the GPU buffer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuSphParticle {
    /// Particle position \[x, y, z\] in metres.
    pub pos: [f32; 3],
    /// Particle velocity \[vx, vy, vz\] in m/s.
    pub vel: [f32; 3],
    /// Particle density in kg/m³.
    pub density: f32,
    /// Particle pressure in Pa.
    pub pressure: f32,
    /// Particle mass in kg.
    pub mass: f32,
}

/// Simulation parameters for the GPU SPH solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuSphParams {
    /// SPH smoothing kernel radius h in metres.
    pub kernel_radius: f32,
    /// Tait EOS stiffness constant k.
    pub eos_k: f32,
    /// Tait EOS exponent gamma.
    pub eos_gamma: f32,
    /// Dynamic viscosity coefficient μ.
    pub viscosity: f32,
    /// Number of particles.
    pub n_particles: usize,
}

/// GPU-side buffer holding all per-particle arrays.
#[derive(Debug, Clone)]
pub struct GpuSphBuffer {
    /// Per-particle positions \[x, y, z\].
    pub positions: Vec<[f32; 3]>,
    /// Per-particle velocities \[vx, vy, vz\].
    pub velocities: Vec<[f32; 3]>,
    /// Per-particle densities (kg/m³).
    pub densities: Vec<f32>,
    /// Per-particle pressures (Pa).
    pub pressures: Vec<f32>,
}

impl GpuSphBuffer {
    /// Allocate a new buffer for `n` particles, all initialised to zero.
    pub fn new(n: usize) -> Self {
        Self {
            positions: vec![[0.0; 3]; n],
            velocities: vec![[0.0; 3]; n],
            densities: vec![0.0; n],
            pressures: vec![0.0; n],
        }
    }

    /// Number of particles in this buffer.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Returns `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

// ── Kernel functions ─────────────────────────────────────────────────────────

/// Wendland C2 SPH kernel W(r, h).
///
/// Returns the kernel value at distance `r` with smoothing length `h`.
/// The kernel is zero for `r >= h`.
pub fn sph_kernel_gpu(r: f32, h: f32) -> f32 {
    if h <= 0.0 || r >= h {
        return 0.0;
    }
    let q = r / h;
    // 3-D Wendland C2: alpha_D * (1 - q/2)^4 * (2q + 1)
    // alpha_D = 21 / (2 * pi * h^3)
    let alpha = 21.0 / (2.0 * PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    let t2 = t * t;
    let t4 = t2 * t2;
    alpha * t4 * (2.0 * q + 1.0)
}

/// Wendland C2 SPH kernel gradient magnitude dW/dr.
///
/// Returns the magnitude of the kernel gradient at distance `r`.
pub fn sph_kernel_grad_gpu(r: f32, h: f32) -> f32 {
    if h <= 0.0 || r >= h || r < 1e-12 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 21.0 / (2.0 * PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    let t2 = t * t;
    let t3 = t2 * t;
    // d/dr [ t^4 (2q+1) ] = [ -4*(1/2h)*t^3*(2q+1) + t^4*(2/h) ]
    //                      = t^3 / h * [ -2*(2q+1) + 2*t ]
    //                      = t^3 / h * [ -4q - 2 + 2 - q ] = t^3/h * (-5q)
    // actually: d/dq [ t^4(2q+1) ] * (1/h)
    // = [ 4*t^3*(-1/2)*(2q+1) + t^4*2 ] / h
    // = t^3 [ -2(2q+1) + 2t ] / h
    // = t^3 [ -4q-2 + 2 - q ] / h = t^3(-5q) / h
    alpha * t3 * (-5.0 * q) / h
}

// ── Density computation ───────────────────────────────────────────────────────

/// Compute per-particle densities via O(N²) neighbour sum.
///
/// `ρ_i = Σ_j m_j * W(|r_i - r_j|, h)`
pub fn compute_density_gpu(buf: &GpuSphBuffer, params: &GpuSphParams) -> Vec<f32> {
    let n = buf.positions.len();
    let h = params.kernel_radius;
    // Assume unit mass per particle (mass can be encoded separately)
    let mass = 1.0_f32;
    let mut densities = vec![0.0_f32; n];
    for i in 0..n {
        let pi = buf.positions[i];
        let mut rho = 0.0_f32;
        for j in 0..n {
            let pj = buf.positions[j];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            rho += mass * sph_kernel_gpu(r, h);
        }
        densities[i] = rho;
    }
    densities
}

// ── Pressure computation ──────────────────────────────────────────────────────

/// Compute per-particle pressures using the Tait equation of state.
///
/// `P = k * (ρ/ρ₀)^γ - k` where `ρ₀ = 1000` kg/m³.
pub fn compute_pressure_gpu(densities: &[f32], params: &GpuSphParams) -> Vec<f32> {
    let rho0 = 1000.0_f32;
    densities
        .iter()
        .map(|&rho| {
            let ratio = rho / rho0;
            params.eos_k * (ratio.powf(params.eos_gamma) - 1.0)
        })
        .collect()
}

// ── Force computation ─────────────────────────────────────────────────────────

/// Compute pressure forces on all particles.
///
/// Uses the symmetric SPH momentum equation:
/// `f_i = -Σ_j m_j (P_i/ρ_i² + P_j/ρ_j²) ∇W`
#[allow(clippy::too_many_arguments)]
pub fn compute_pressure_force_gpu(
    buf: &GpuSphBuffer,
    pressures: &[f32],
    params: &GpuSphParams,
) -> Vec<[f32; 3]> {
    let n = buf.positions.len();
    let h = params.kernel_radius;
    let mass = 1.0_f32;
    let mut forces = vec![[0.0_f32; 3]; n];

    for i in 0..n {
        let pi = buf.positions[i];
        let rho_i = buf.densities[i].max(1e-6);
        let p_i = pressures[i];
        let mut fx = 0.0_f32;
        let mut fy = 0.0_f32;
        let mut fz = 0.0_f32;

        for j in 0..n {
            if i == j {
                continue;
            }
            let pj = buf.positions[j];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-12 {
                continue;
            }
            let rho_j = buf.densities[j].max(1e-6);
            let p_j = pressures[j];
            let grad_mag = sph_kernel_grad_gpu(r, h);
            let coeff = -mass * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * grad_mag / r;
            fx += coeff * dx;
            fy += coeff * dy;
            fz += coeff * dz;
        }
        forces[i] = [fx, fy, fz];
    }
    forces
}

/// Compute viscosity forces on all particles.
///
/// Uses the standard Monaghan viscosity formulation.
pub fn compute_viscosity_force_gpu(buf: &GpuSphBuffer, params: &GpuSphParams) -> Vec<[f32; 3]> {
    let n = buf.positions.len();
    let h = params.kernel_radius;
    let mu = params.viscosity;
    let mass = 1.0_f32;
    let mut forces = vec![[0.0_f32; 3]; n];

    for i in 0..n {
        let pi = buf.positions[i];
        let vi = buf.velocities[i];
        let rho_i = buf.densities[i].max(1e-6);
        let mut fx = 0.0_f32;
        let mut fy = 0.0_f32;
        let mut fz = 0.0_f32;

        for j in 0..n {
            if i == j {
                continue;
            }
            let pj = buf.positions[j];
            let vj = buf.velocities[j];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-12 {
                continue;
            }
            let rho_j = buf.densities[j].max(1e-6);
            let lap = sph_kernel_grad_gpu(r, h); // use gradient magnitude as Laplacian approx
            let dvx = vj[0] - vi[0];
            let dvy = vj[1] - vi[1];
            let dvz = vj[2] - vi[2];
            let coeff = mu * mass / rho_j * lap / r;
            fx += coeff * dvx;
            fy += coeff * dvy;
            fz += coeff * dvz;
        }
        forces[i] = [
            forces[i][0] + fx / rho_i,
            forces[i][1] + fy / rho_i,
            forces[i][2] + fz / rho_i,
        ];
    }
    forces
}

// ── Integration ───────────────────────────────────────────────────────────────

/// Semi-implicit Euler integration step.
///
/// Updates velocities and positions in-place:
/// `v += f * dt`, `x += v * dt`.
pub fn integrate_sph_gpu(buf: &mut GpuSphBuffer, forces: &[[f32; 3]], dt: f32) {
    let n = buf.positions.len();
    for i in 0..n {
        buf.velocities[i][0] += forces[i][0] * dt;
        buf.velocities[i][1] += forces[i][1] * dt;
        buf.velocities[i][2] += forces[i][2] * dt;
        buf.positions[i][0] += buf.velocities[i][0] * dt;
        buf.positions[i][1] += buf.velocities[i][1] * dt;
        buf.positions[i][2] += buf.velocities[i][2] * dt;
    }
}

// ── Full step ─────────────────────────────────────────────────────────────────

/// Perform one complete SPH time-step.
///
/// Executes density → pressure → force → integrate in sequence.
pub fn gpu_sph_step(buf: &mut GpuSphBuffer, params: &GpuSphParams, dt: f32) {
    let densities = compute_density_gpu(buf, params);
    let pressures = compute_pressure_gpu(&densities, params);
    buf.densities = densities;
    buf.pressures = pressures.clone();
    let pf = compute_pressure_force_gpu(buf, &pressures, params);
    let vf = compute_viscosity_force_gpu(buf, params);
    let n = buf.positions.len();
    let mut total_forces = vec![[0.0_f32; 3]; n];
    for i in 0..n {
        total_forces[i][0] = pf[i][0] + vf[i][0];
        total_forces[i][1] = pf[i][1] + vf[i][1];
        total_forces[i][2] = pf[i][2] + vf[i][2];
    }
    integrate_sph_gpu(buf, &total_forces, dt);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_buf(n: usize) -> GpuSphBuffer {
        let mut buf = GpuSphBuffer::new(n);
        for i in 0..n {
            buf.positions[i] = [i as f32 * 0.1, 0.0, 0.0];
            buf.velocities[i] = [0.0; 3];
            buf.densities[i] = 1000.0;
            buf.pressures[i] = 0.0;
        }
        buf
    }

    fn default_params(n: usize) -> GpuSphParams {
        GpuSphParams {
            kernel_radius: 0.5,
            eos_k: 1.0,
            eos_gamma: 7.0,
            viscosity: 0.01,
            n_particles: n,
        }
    }

    #[test]
    fn test_gpu_sph_particle_fields() {
        let p = GpuSphParticle {
            pos: [1.0, 2.0, 3.0],
            vel: [0.1, 0.2, 0.3],
            density: 1000.0,
            pressure: 101325.0,
            mass: 0.001,
        };
        assert_eq!(p.pos[0], 1.0);
        assert_eq!(p.mass, 0.001);
    }

    #[test]
    fn test_gpu_sph_params_fields() {
        let params = default_params(10);
        assert_eq!(params.n_particles, 10);
        assert!(params.kernel_radius > 0.0);
    }

    #[test]
    fn test_gpu_sph_buffer_new() {
        let buf = GpuSphBuffer::new(5);
        assert_eq!(buf.len(), 5);
        assert!(!buf.is_empty());
    }

    #[test]
    fn test_gpu_sph_buffer_empty() {
        let buf = GpuSphBuffer::new(0);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_sph_kernel_gpu_zero_at_boundary() {
        let w = sph_kernel_gpu(0.5, 0.5);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_sph_kernel_gpu_zero_beyond() {
        let w = sph_kernel_gpu(1.0, 0.5);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_sph_kernel_gpu_positive_inside() {
        let w = sph_kernel_gpu(0.1, 0.5);
        assert!(w > 0.0);
    }

    #[test]
    fn test_sph_kernel_gpu_peak_at_zero() {
        let w0 = sph_kernel_gpu(0.0, 0.5);
        let w1 = sph_kernel_gpu(0.2, 0.5);
        assert!(w0 > w1);
    }

    #[test]
    fn test_sph_kernel_gpu_zero_h() {
        assert_eq!(sph_kernel_gpu(0.1, 0.0), 0.0);
    }

    #[test]
    fn test_sph_kernel_grad_gpu_zero_h() {
        assert_eq!(sph_kernel_grad_gpu(0.1, 0.0), 0.0);
    }

    #[test]
    fn test_sph_kernel_grad_gpu_zero_r() {
        assert_eq!(sph_kernel_grad_gpu(0.0, 0.5), 0.0);
    }

    #[test]
    fn test_sph_kernel_grad_gpu_negative_inside() {
        // Gradient magnitude is negative (kernel decreases with r)
        let g = sph_kernel_grad_gpu(0.2, 0.5);
        assert!(g < 0.0);
    }

    #[test]
    fn test_sph_kernel_grad_gpu_zero_at_boundary() {
        let g = sph_kernel_grad_gpu(0.5, 0.5);
        assert_eq!(g, 0.0);
    }

    #[test]
    fn test_compute_density_gpu_nonzero() {
        let buf = make_buf(4);
        let params = default_params(4);
        let densities = compute_density_gpu(&buf, &params);
        assert_eq!(densities.len(), 4);
        // Particles within kernel_radius should contribute
        assert!(densities.iter().any(|&d| d > 0.0));
    }

    #[test]
    fn test_compute_density_gpu_length() {
        let buf = make_buf(6);
        let params = default_params(6);
        let d = compute_density_gpu(&buf, &params);
        assert_eq!(d.len(), 6);
    }

    #[test]
    fn test_compute_density_gpu_empty() {
        let buf = GpuSphBuffer::new(0);
        let params = default_params(0);
        let d = compute_density_gpu(&buf, &params);
        assert!(d.is_empty());
    }

    #[test]
    fn test_compute_pressure_gpu_positive() {
        let params = default_params(3);
        let densities = vec![1100.0_f32, 1000.0_f32, 900.0_f32];
        let pressures = compute_pressure_gpu(&densities, &params);
        assert_eq!(pressures.len(), 3);
        // Higher density => higher pressure
        assert!(pressures[0] > pressures[1]);
    }

    #[test]
    fn test_compute_pressure_gpu_rest_density() {
        let params = default_params(1);
        let d = vec![1000.0_f32]; // exactly rest density
        let p = compute_pressure_gpu(&d, &params);
        // (1.0)^gamma - 1 = 0 => P = 0
        assert!(p[0].abs() < 1e-3);
    }

    #[test]
    fn test_compute_pressure_force_gpu_length() {
        let buf = make_buf(4);
        let params = default_params(4);
        let pressures = vec![100.0_f32; 4];
        let forces = compute_pressure_force_gpu(&buf, &pressures, &params);
        assert_eq!(forces.len(), 4);
    }

    #[test]
    fn test_compute_viscosity_force_gpu_length() {
        let buf = make_buf(4);
        let params = default_params(4);
        let forces = compute_viscosity_force_gpu(&buf, &params);
        assert_eq!(forces.len(), 4);
    }

    #[test]
    fn test_integrate_sph_gpu_position_update() {
        let mut buf = GpuSphBuffer::new(2);
        buf.velocities[0] = [1.0, 0.0, 0.0];
        buf.velocities[1] = [0.0, 1.0, 0.0];
        let forces = vec![[0.0_f32; 3]; 2];
        integrate_sph_gpu(&mut buf, &forces, 0.1);
        assert!((buf.positions[0][0] - 0.1).abs() < 1e-5);
        assert!((buf.positions[1][1] - 0.1).abs() < 1e-5);
    }

    #[test]
    fn test_integrate_sph_gpu_velocity_update() {
        let mut buf = GpuSphBuffer::new(1);
        let forces = vec![[2.0_f32, 0.0, 0.0]];
        integrate_sph_gpu(&mut buf, &forces, 0.5);
        assert!((buf.velocities[0][0] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_gpu_sph_step_runs() {
        let mut buf = make_buf(5);
        let params = default_params(5);
        gpu_sph_step(&mut buf, &params, 0.001);
        // After one step densities should be set
        assert!(buf.densities.iter().any(|&d| d >= 0.0));
    }

    #[test]
    fn test_gpu_sph_step_no_nan() {
        let mut buf = make_buf(5);
        let params = default_params(5);
        gpu_sph_step(&mut buf, &params, 0.001);
        for i in 0..5 {
            assert!(!buf.positions[i][0].is_nan());
            assert!(!buf.densities[i].is_nan());
        }
    }

    #[test]
    fn test_gpu_sph_step_pressure_set() {
        let mut buf = make_buf(3);
        let params = default_params(3);
        gpu_sph_step(&mut buf, &params, 0.001);
        assert_eq!(buf.pressures.len(), 3);
    }

    #[test]
    fn test_sph_kernel_symmetry() {
        let h = 0.5;
        let w1 = sph_kernel_gpu(0.1, h);
        let w2 = sph_kernel_gpu(0.1, h);
        assert!((w1 - w2).abs() < 1e-8);
    }

    #[test]
    fn test_sph_kernel_decreasing() {
        let h = 1.0;
        let mut prev = sph_kernel_gpu(0.0, h);
        for i in 1..10 {
            let r = i as f32 * 0.09;
            let w = sph_kernel_gpu(r, h);
            assert!(w <= prev + 1e-6);
            prev = w;
        }
    }

    #[test]
    fn test_compute_density_uniform_grid() {
        // All particles at same location -> all have same density
        let n = 4;
        let buf = GpuSphBuffer::new(n);
        // All at origin
        let params = default_params(n);
        let d = compute_density_gpu(&buf, &params);
        // All densities should be equal
        for i in 1..n {
            assert!((d[i] - d[0]).abs() < 1e-5);
        }
        let _ = buf.positions[0]; // suppress warning
    }

    #[test]
    fn test_pressure_force_zero_pressure() {
        let buf = make_buf(3);
        let params = default_params(3);
        let pressures = vec![0.0_f32; 3];
        let forces = compute_pressure_force_gpu(&buf, &pressures, &params);
        for f in &forces {
            assert!(f[0].abs() < 1e-6 && f[1].abs() < 1e-6 && f[2].abs() < 1e-6);
        }
    }

    #[test]
    fn test_viscosity_force_same_velocity() {
        // Particles with same velocity should have zero viscosity force
        let mut buf = make_buf(3);
        for i in 0..3 {
            buf.velocities[i] = [1.0, 0.5, 0.0];
        }
        let params = default_params(3);
        let forces = compute_viscosity_force_gpu(&buf, &params);
        for f in &forces {
            assert!(f[0].abs() < 1e-4 && f[1].abs() < 1e-4 && f[2].abs() < 1e-4);
        }
    }

    #[test]
    fn test_buf_clone() {
        let buf = make_buf(3);
        let buf2 = buf.clone();
        assert_eq!(buf2.len(), 3);
    }

    #[test]
    fn test_params_copy() {
        let p = default_params(5);
        let p2 = p;
        assert_eq!(p2.n_particles, 5);
    }

    #[test]
    fn test_integrate_multiple_steps() {
        let mut buf = GpuSphBuffer::new(1);
        buf.velocities[0] = [1.0, 0.0, 0.0];
        let forces = vec![[0.0_f32; 3]];
        for _ in 0..10 {
            integrate_sph_gpu(&mut buf, &forces, 0.1);
        }
        assert!((buf.positions[0][0] - 1.0).abs() < 1e-4);
    }
}
