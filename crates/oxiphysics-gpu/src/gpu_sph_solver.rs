// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-side SPH solver — CPU mock backend.
//!
//! Implements a complete Smoothed Particle Hydrodynamics (SPH) solver pipeline
//! using plain Rust loops as a CPU fallback. The API mirrors what would run on
//! a GPU kernel dispatch, making it straightforward to swap in a real GPU
//! backend without changing call-sites.

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

/// Wendland C2 SPH kernel W(r, h) with compact support radius `h`.
///
/// In 3-D, with `q = r / h`:
///
/// `W(q) = 21 / (2π h³) · (1 − q)⁴ · (4q + 1)` for `0 ≤ q < 1`, and `0` otherwise.
///
/// The kernel is C² across the support boundary (`W`, `dW/dr` and `d²W/dr²` all
/// vanish at `r = h`) and normalised, `∫ W dV = 1`.
pub fn sph_kernel_gpu(r: f32, h: f32) -> f32 {
    if h <= 0.0 || r >= h {
        return 0.0;
    }
    let q = r / h;
    let alpha = 21.0 / (2.0 * PI * h * h * h);
    let t = 1.0 - q;
    let t2 = t * t;
    alpha * t2 * t2 * (4.0 * q + 1.0)
}

/// Radial derivative `dW/dr` of the Wendland C2 kernel (support radius `h`).
///
/// `dW/dr = 21 / (2π h³) · (−20 q (1 − q)³) / h` for `0 < r < h` with
/// `q = r / h`; the value is non-positive.  Returns `0` for `r ≥ h` and for
/// (near-)coincident points where the gradient direction is undefined.
pub fn sph_kernel_grad_gpu(r: f32, h: f32) -> f32 {
    if h <= 0.0 || r >= h || r < 1e-12 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 21.0 / (2.0 * PI * h * h * h);
    let t = 1.0 - q;
    alpha * (-20.0 * q * t * t * t) / h
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
    for (rho, &pi) in densities.iter_mut().zip(buf.positions.iter()) {
        let mut acc = 0.0_f32;
        for j in 0..n {
            let pj = buf.positions[j];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            acc += mass * sph_kernel_gpu(r, h);
        }
        *rho = acc;
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
pub fn compute_pressure_force_gpu(
    buf: &GpuSphBuffer,
    pressures: &[f32],
    params: &GpuSphParams,
) -> Vec<[f32; 3]> {
    let n = buf.positions.len();
    let h = params.kernel_radius;
    let mass = 1.0_f32;
    let mut forces = vec![[0.0_f32; 3]; n];

    for (i, (force, &pi)) in forces.iter_mut().zip(buf.positions.iter()).enumerate() {
        let rho_i = buf.densities[i].max(1e-6);
        let p_i = pressures[i];
        let mut fx = 0.0_f32;
        let mut fy = 0.0_f32;
        let mut fz = 0.0_f32;

        for (j, &pj) in buf.positions.iter().enumerate() {
            if i == j {
                continue;
            }
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
        *force = [fx, fy, fz];
    }
    forces
}

/// Compute viscous accelerations on all particles.
///
/// Uses the Morris et al. (1997) Laplacian formulation for a constant dynamic
/// viscosity `μ`:
///
/// `a_i = Σ_j 2 μ m_j / (ρ_i ρ_j) · (r_ij · ∇_i W_ij) / (|r_ij|² + η²) · (v_i − v_j)`
///
/// with `r_ij = r_i − r_j`, `r_ij · ∇_i W_ij = r · dW/dr` and `η² = 0.01 h²`.
/// Because `dW/dr ≤ 0` every pair term opposes the relative velocity, so the
/// viscosity damps shear instead of amplifying it.  The result is an
/// acceleration (it is integrated as `v += a·dt` by [`integrate_sph_gpu`]).
pub fn compute_viscosity_force_gpu(buf: &GpuSphBuffer, params: &GpuSphParams) -> Vec<[f32; 3]> {
    let n = buf.positions.len();
    let h = params.kernel_radius;
    let mu = params.viscosity;
    let mass = 1.0_f32;
    let eta2 = 0.01 * h * h;
    let mut forces = vec![[0.0_f32; 3]; n];

    for (i, (force, (&pi, &vi))) in forces
        .iter_mut()
        .zip(buf.positions.iter().zip(buf.velocities.iter()))
        .enumerate()
    {
        let rho_i = buf.densities[i].max(1e-6);
        let mut fx = 0.0_f32;
        let mut fy = 0.0_f32;
        let mut fz = 0.0_f32;

        for (j, (&pj, &vj)) in buf.positions.iter().zip(buf.velocities.iter()).enumerate() {
            if i == j {
                continue;
            }
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            if r < 1e-12 {
                continue;
            }
            let rho_j = buf.densities[j].max(1e-6);
            let dw_dr = sph_kernel_grad_gpu(r, h);
            let coeff = 2.0 * mu * mass * (r * dw_dr) / (rho_i * rho_j * (r2 + eta2));
            fx += coeff * (vi[0] - vj[0]);
            fy += coeff * (vi[1] - vj[1]);
            fz += coeff * (vi[2] - vj[2]);
        }
        *force = [fx, fy, fz];
    }
    forces
}

// ── Integration ───────────────────────────────────────────────────────────────

/// Semi-implicit Euler integration step.
///
/// Updates velocities and positions in-place:
/// `v += f * dt`, `x += v * dt`.
pub fn integrate_sph_gpu(buf: &mut GpuSphBuffer, forces: &[[f32; 3]], dt: f32) {
    for (vel, (pos, &f)) in buf
        .velocities
        .iter_mut()
        .zip(buf.positions.iter_mut().zip(forces.iter()))
    {
        vel[0] += f[0] * dt;
        vel[1] += f[1] * dt;
        vel[2] += f[2] * dt;
        pos[0] += vel[0] * dt;
        pos[1] += vel[1] * dt;
        pos[2] += vel[2] * dt;
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
    // ── Wendland C2 kernel with support radius h ─────────────────────────

    /// f64 reference of the 3-D Wendland C2 kernel (support radius `h`).
    fn wendland_reference(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 1.0 {
            return 0.0;
        }
        21.0 / (2.0 * std::f64::consts::PI * h.powi(3)) * (1.0 - q).powi(4) * (4.0 * q + 1.0)
    }

    #[test]
    fn test_sph_kernel_gpu_is_normalised_over_its_support() {
        // Simpson quadrature of the volume integral ∫₀ʰ 4π r² W(r) dr must be 1.
        let h = 0.5_f32;
        let n = 2000usize; // even number of sub-intervals
        let dr = h as f64 / n as f64;
        let integrand =
            |r: f64| 4.0 * std::f64::consts::PI * r * r * sph_kernel_gpu(r as f32, h) as f64;
        let mut acc = integrand(0.0) + integrand(h as f64);
        for i in 1..n {
            let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
            acc += weight * integrand(i as f64 * dr);
        }
        let integral = acc * dr / 3.0;
        assert!(
            (integral - 1.0).abs() < 1e-3,
            "kernel must integrate to 1 over its support, got {integral}"
        );
    }

    #[test]
    fn test_sph_kernel_gpu_vanishes_continuously_at_support_boundary() {
        // W must vanish continuously at the support boundary r = h (no jump).
        let h = 0.5_f32;
        let w0 = sph_kernel_gpu(0.0, h);
        let w_edge = sph_kernel_gpu(0.999 * h, h);
        assert!(w0 > 0.0);
        assert!(
            w_edge < 1e-6 * w0,
            "W(0.999h) = {w_edge} must be negligible next to W(0) = {w0}"
        );
        let g_edge = sph_kernel_grad_gpu(0.999 * h, h).abs();
        let g_peak = sph_kernel_grad_gpu(0.25 * h, h).abs();
        assert!(
            g_edge < 1e-6 * g_peak,
            "dW/dr must vanish at the support edge"
        );
    }

    #[test]
    fn test_sph_kernel_gpu_matches_reference_formula() {
        let h = 0.37_f32;
        for i in 0..100 {
            let r = i as f32 * 0.01 * h;
            let expected = wendland_reference(r as f64, h as f64);
            let got = sph_kernel_gpu(r, h) as f64;
            assert!(
                (got - expected).abs() <= 1e-5 * expected.abs().max(1e-3),
                "W({r}, {h}) = {got}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_sph_kernel_grad_gpu_matches_finite_difference_of_kernel() {
        let h = 0.4_f64;
        for &q in &[0.1_f64, 0.25, 0.4, 0.55, 0.7, 0.85] {
            let r = q * h;
            let eps = 1e-6;
            let numeric =
                (wendland_reference(r + eps, h) - wendland_reference(r - eps, h)) / (2.0 * eps);
            let analytic = sph_kernel_grad_gpu(r as f32, h as f32) as f64;
            assert!(
                analytic < 0.0,
                "dW/dr must be negative inside the support (q = {q})"
            );
            assert!(
                (analytic - numeric).abs() <= 2e-4 * numeric.abs(),
                "dW/dr at q = {q}: analytic {analytic} vs finite difference {numeric}"
            );
        }
    }

    // ── Morris viscosity ─────────────────────────────────────────────────

    /// Two particles `0.5 h` apart along x with opposite velocities along y.
    fn shear_pair(h: f32, viscosity: f32) -> (GpuSphBuffer, GpuSphParams) {
        let mut buf = GpuSphBuffer::new(2);
        buf.positions[1] = [0.5 * h, 0.0, 0.0];
        buf.velocities[0] = [0.0, 1.0, 0.0];
        buf.velocities[1] = [0.0, -1.0, 0.0];
        let params = GpuSphParams {
            kernel_radius: h,
            eos_k: 1.0,
            eos_gamma: 7.0,
            viscosity,
            n_particles: 2,
        };
        buf.densities = compute_density_gpu(&buf, &params);
        (buf, params)
    }

    #[test]
    fn test_viscosity_acceleration_matches_morris_formula() {
        let h = 0.2_f32;
        let mu = 0.5_f32;
        let (buf, params) = shear_pair(h, mu);
        let acc = compute_viscosity_force_gpu(&buf, &params);

        let (h64, mu64) = (h as f64, mu as f64);
        let r = 0.5 * h64;
        let rho_i = buf.densities[0] as f64;
        let rho_j = buf.densities[1] as f64;
        let dw_dr = {
            let eps = 1e-7;
            (wendland_reference(r + eps, h64) - wendland_reference(r - eps, h64)) / (2.0 * eps)
        };
        let coeff = 2.0 * mu64 * 1.0 / (rho_i * rho_j) * (r * dw_dr) / (r * r + 0.01 * h64 * h64);
        assert!(coeff < 0.0, "Morris coefficient must be negative (damping)");
        // a_i = coeff (v_i - v_j) = coeff * (0, 2, 0); a_j = -a_i.
        let expected_y = coeff * 2.0;
        assert!(
            (acc[0][1] as f64 - expected_y).abs() <= 1e-3 * expected_y.abs(),
            "a_i.y = {} but the Morris formula gives {expected_y}",
            acc[0][1]
        );
        assert!((acc[1][1] as f64 + expected_y).abs() <= 1e-3 * expected_y.abs());
        assert!(acc[0][0].abs() < 1e-12 && acc[0][2].abs() < 1e-12);
        // Equal masses and densities: momentum is conserved pairwise.
        assert!((acc[0][1] + acc[1][1]).abs() < 1e-6 * acc[0][1].abs());
    }

    #[test]
    fn test_viscosity_damps_relative_velocity_without_overshoot() {
        let h = 0.2_f32;
        let (mut buf, params) = shear_pair(h, 0.5);
        let initial_positions = buf.positions.clone();
        let mut previous = (buf.velocities[0][1] - buf.velocities[1][1]).abs();
        assert!((previous - 2.0).abs() < 1e-6);
        for step in 0..10 {
            let acc = compute_viscosity_force_gpu(&buf, &params);
            integrate_sph_gpu(&mut buf, &acc, 0.1);
            // Isolate the viscous effect from advection.
            buf.positions = initial_positions.clone();
            let relative = (buf.velocities[0][1] - buf.velocities[1][1]).abs();
            assert!(
                relative < previous,
                "step {step}: |v_i - v_j| must strictly decrease ({previous} -> {relative})"
            );
            assert!(
                buf.velocities[0][1] > 0.0 && buf.velocities[1][1] < 0.0,
                "step {step}: viscosity must not reverse the shear direction"
            );
            previous = relative;
        }
        assert!(
            previous < 1.9,
            "ten damping steps should visibly reduce shear"
        );
    }

    #[test]
    fn test_viscosity_never_accelerates_away_from_neighbour_velocity() {
        // A 1-D line of particles with a linear velocity ramp: viscosity must pull
        // every particle towards its neighbours' velocities, i.e. the slowest end
        // is accelerated upwards and the fastest end downwards.
        let n = 6;
        let h = 0.5_f32;
        let mut buf = GpuSphBuffer::new(n);
        for i in 0..n {
            buf.positions[i] = [0.1 * i as f32, 0.0, 0.0];
            buf.velocities[i] = [0.0, (i as f32 - 2.5) * 0.3, 0.0];
        }
        let params = GpuSphParams {
            kernel_radius: h,
            eos_k: 1.0,
            eos_gamma: 7.0,
            viscosity: 0.2,
            n_particles: n,
        };
        buf.densities = compute_density_gpu(&buf, &params);
        let acc = compute_viscosity_force_gpu(&buf, &params);
        // The two end particles sit at the extremes of the ramp: the lowest-speed
        // end must be accelerated upwards and the highest downwards.
        assert!(
            acc[0][1] > 0.0,
            "slowest particle must be pulled up, got {}",
            acc[0][1]
        );
        assert!(
            acc[n - 1][1] < 0.0,
            "fastest particle must be pulled down, got {}",
            acc[n - 1][1]
        );
    }
}
