// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-side molecular dynamics (MD) solver — CPU mock backend.
//!
//! Implements a Lennard-Jones MD solver pipeline using plain Rust loops as
//! a CPU fallback. Periodic boundary conditions (minimum image convention)
//! are applied. The API mirrors a GPU kernel dispatch for easy substitution.

// ── Data structures ──────────────────────────────────────────────────────────

/// A single MD atom stored in the GPU buffer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuMdAtom {
    /// Atom position \[x, y, z\] in Angstroms.
    pub pos: [f32; 3],
    /// Atom velocity \[vx, vy, vz\] in Å/ps.
    pub vel: [f32; 3],
    /// Current force on atom \[fx, fy, fz\] in kJ/(mol·Å).
    pub force: [f32; 3],
    /// Atom mass in atomic mass units (amu).
    pub mass: f32,
    /// Partial charge in elementary charge units.
    pub charge: f32,
}

/// Simulation parameters for the GPU MD solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuMdParams {
    /// Lennard-Jones well depth ε (kJ/mol).
    pub epsilon: f32,
    /// Lennard-Jones radius σ (Å).
    pub sigma: f32,
    /// Interaction cutoff distance (Å).
    pub cutoff: f32,
    /// Periodic simulation box dimensions \[Lx, Ly, Lz\] in Å.
    pub box_size: [f32; 3],
    /// Number of atoms.
    pub n_atoms: usize,
}

/// GPU-side buffer holding the full MD system state.
#[derive(Debug, Clone)]
pub struct GpuMdBuffer {
    /// Per-atom data.
    pub atoms: Vec<GpuMdAtom>,
    /// Current simulation step counter.
    pub step: usize,
}

impl GpuMdBuffer {
    /// Allocate a new buffer for `n` atoms, all at origin with zero velocity.
    pub fn new(n: usize) -> Self {
        Self {
            atoms: vec![
                GpuMdAtom {
                    pos: [0.0; 3],
                    vel: [0.0; 3],
                    force: [0.0; 3],
                    mass: 1.0,
                    charge: 0.0,
                };
                n
            ],
            step: 0,
        }
    }

    /// Number of atoms in this buffer.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Returns `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
}

// ── LJ kernel functions ───────────────────────────────────────────────────────

/// Lennard-Jones 12-6 force magnitude dU/dr (negative of force along r̂).
///
/// Returns `dU/dr = 4ε [ -12σ¹²/r¹³ + 6σ⁶/r⁷ ]`.
/// Returns zero for `r >= cutoff` or `r < 1e-10`.
pub fn lj_force_gpu(r: f32, params: &GpuMdParams) -> f32 {
    if r >= params.cutoff || r < 1e-10 {
        return 0.0;
    }
    let sr = params.sigma / r;
    let sr6 = sr * sr * sr * sr * sr * sr;
    let sr12 = sr6 * sr6;
    // Return the force F = -dU/dr = 4ε(12σ¹²/r¹³ - 6σ⁶/r⁷)
    4.0 * params.epsilon * (12.0 * sr12 - 6.0 * sr6) / r
}

/// Lennard-Jones 12-6 pair potential.
///
/// Returns `U(r) = 4ε [ (σ/r)¹² − (σ/r)⁶ ]`.
/// Returns zero for `r >= cutoff`.
pub fn lj_potential_gpu(r: f32, params: &GpuMdParams) -> f32 {
    if r >= params.cutoff || r < 1e-10 {
        return 0.0;
    }
    let sr = params.sigma / r;
    let sr6 = sr * sr * sr * sr * sr * sr;
    let sr12 = sr6 * sr6;
    4.0 * params.epsilon * (sr12 - sr6)
}

// ── Periodic boundary conditions ──────────────────────────────────────────────

/// Minimum-image distance between two atoms under periodic boundary conditions.
///
/// Applies the minimum image convention along each box axis.
pub fn pbc_distance_gpu(a: [f32; 3], b: [f32; 3], box_size: [f32; 3]) -> f32 {
    let mut r2 = 0.0_f32;
    for k in 0..3 {
        let mut d = a[k] - b[k];
        let l = box_size[k];
        if l > 0.0 {
            d -= (d / l).round() * l;
        }
        r2 += d * d;
    }
    r2.sqrt()
}

/// Minimum-image displacement vector from atom `b` to atom `a`.
fn pbc_displacement_gpu(a: [f32; 3], b: [f32; 3], box_size: [f32; 3]) -> [f32; 3] {
    let mut disp = [0.0_f32; 3];
    for k in 0..3 {
        let mut d = a[k] - b[k];
        let l = box_size[k];
        if l > 0.0 {
            d -= (d / l).round() * l;
        }
        disp[k] = d;
    }
    disp
}

// ── Force computation ─────────────────────────────────────────────────────────

/// Compute all-pairs LJ forces on every atom and store in `atom.force`.
///
/// Clears previous forces before accumulating new ones.
pub fn compute_forces_gpu(buf: &mut GpuMdBuffer, params: &GpuMdParams) {
    let n = buf.atoms.len();
    // Zero forces
    for atom in buf.atoms.iter_mut() {
        atom.force = [0.0; 3];
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let pi = buf.atoms[i].pos;
            let pj = buf.atoms[j].pos;
            let disp = pbc_displacement_gpu(pi, pj, params.box_size);
            let r2 = disp[0] * disp[0] + disp[1] * disp[1] + disp[2] * disp[2];
            let r = r2.sqrt();
            if r < 1e-10 || r >= params.cutoff {
                continue;
            }
            let f_mag = lj_force_gpu(r, params);
            // Force on i from j: F = f_mag * (disp / r) — but lj_force is dU/dr
            // F_i = -(dU/dr) * r̂ = -(dU/dr) * disp/r
            let scale = -f_mag / r;
            buf.atoms[i].force[0] += scale * disp[0];
            buf.atoms[i].force[1] += scale * disp[1];
            buf.atoms[i].force[2] += scale * disp[2];
            buf.atoms[j].force[0] -= scale * disp[0];
            buf.atoms[j].force[1] -= scale * disp[1];
            buf.atoms[j].force[2] -= scale * disp[2];
        }
    }
}

// ── Integration ───────────────────────────────────────────────────────────────

/// Velocity-Verlet integration step.
///
/// Updates positions and velocities using the current forces.
/// `v += (F/m) * dt`, `x += v * dt`.
/// Increments the step counter.
pub fn verlet_integrate_gpu(buf: &mut GpuMdBuffer, dt: f32) {
    for atom in buf.atoms.iter_mut() {
        let inv_mass = if atom.mass > 1e-10 {
            1.0 / atom.mass
        } else {
            0.0
        };
        atom.vel[0] += atom.force[0] * inv_mass * dt;
        atom.vel[1] += atom.force[1] * inv_mass * dt;
        atom.vel[2] += atom.force[2] * inv_mass * dt;
        atom.pos[0] += atom.vel[0] * dt;
        atom.pos[1] += atom.vel[1] * dt;
        atom.pos[2] += atom.vel[2] * dt;
    }
    buf.step += 1;
}

// ── Thermodynamic observables ─────────────────────────────────────────────────

/// Compute total kinetic energy: `KE = Σ 0.5 * m_i * |v_i|²`.
pub fn kinetic_energy_gpu(buf: &GpuMdBuffer) -> f32 {
    buf.atoms
        .iter()
        .map(|a| 0.5 * a.mass * (a.vel[0] * a.vel[0] + a.vel[1] * a.vel[1] + a.vel[2] * a.vel[2]))
        .sum()
}

/// Compute total LJ potential energy.
pub fn potential_energy_gpu(buf: &GpuMdBuffer, params: &GpuMdParams) -> f32 {
    let n = buf.atoms.len();
    let mut pe = 0.0_f32;
    for i in 0..n {
        for j in (i + 1)..n {
            let r = pbc_distance_gpu(buf.atoms[i].pos, buf.atoms[j].pos, params.box_size);
            pe += lj_potential_gpu(r, params);
        }
    }
    pe
}

/// Estimate instantaneous temperature from kinetic energy.
///
/// Uses the equipartition theorem: `T = 2 * KE / (3 * N * kB)`.
/// Boltzmann constant `kB = 8.314e-3` kJ/(mol·K).
pub fn temperature_gpu(buf: &GpuMdBuffer) -> f32 {
    let n = buf.atoms.len();
    if n == 0 {
        return 0.0;
    }
    let kb = 8.314e-3_f32; // kJ/(mol·K)
    let ke = kinetic_energy_gpu(buf);
    2.0 * ke / (3.0 * n as f32 * kb)
}

// ── Thermostat ────────────────────────────────────────────────────────────────

/// Rescale all atom velocities to achieve `target_temp` (velocity scaling).
///
/// No-op when the current temperature is zero.
pub fn rescale_velocities_gpu(buf: &mut GpuMdBuffer, target_temp: f32) {
    let t_curr = temperature_gpu(buf);
    if t_curr < 1e-10 || target_temp < 0.0 {
        return;
    }
    let scale = (target_temp / t_curr).sqrt();
    for atom in buf.atoms.iter_mut() {
        atom.vel[0] *= scale;
        atom.vel[1] *= scale;
        atom.vel[2] *= scale;
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> GpuMdParams {
        GpuMdParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 3.5,
            box_size: [10.0; 3],
            n_atoms: 4,
        }
    }

    fn make_buf_grid(n: usize) -> GpuMdBuffer {
        let mut buf = GpuMdBuffer::new(n);
        for i in 0..n {
            buf.atoms[i].pos = [i as f32 * 1.5, 0.0, 0.0];
            buf.atoms[i].mass = 1.0;
        }
        buf
    }

    #[test]
    fn test_gpu_md_atom_fields() {
        let a = GpuMdAtom {
            pos: [1.0, 2.0, 3.0],
            vel: [0.1, 0.2, 0.3],
            force: [0.0; 3],
            mass: 12.0,
            charge: -0.5,
        };
        assert_eq!(a.mass, 12.0);
        assert_eq!(a.charge, -0.5);
    }

    #[test]
    fn test_gpu_md_params_fields() {
        let p = default_params();
        assert_eq!(p.n_atoms, 4);
        assert!(p.cutoff > p.sigma);
    }

    #[test]
    fn test_gpu_md_buffer_new() {
        let buf = GpuMdBuffer::new(5);
        assert_eq!(buf.len(), 5);
        assert!(!buf.is_empty());
        assert_eq!(buf.step, 0);
    }

    #[test]
    fn test_gpu_md_buffer_empty() {
        let buf = GpuMdBuffer::new(0);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_lj_potential_minimum() {
        let params = default_params();
        // LJ minimum at r = 2^(1/6) * sigma
        let r_min = 2.0_f32.powf(1.0 / 6.0) * params.sigma;
        let u = lj_potential_gpu(r_min, &params);
        // At minimum: U = -epsilon
        assert!((u - (-params.epsilon)).abs() < 0.01);
    }

    #[test]
    fn test_lj_potential_zero_beyond_cutoff() {
        let params = default_params();
        assert_eq!(lj_potential_gpu(params.cutoff + 0.1, &params), 0.0);
    }

    #[test]
    fn test_lj_potential_zero_near_zero_r() {
        let params = default_params();
        assert_eq!(lj_potential_gpu(0.0, &params), 0.0);
    }

    #[test]
    fn test_lj_force_repulsive_close() {
        let params = default_params();
        // At r < sigma the force is repulsive (dU/dr > 0)
        let f = lj_force_gpu(0.8, &params);
        assert!(f > 0.0);
    }

    #[test]
    fn test_lj_force_attractive_far() {
        let params = default_params();
        // At r between sigma and 2^(1/6)*sigma force is attractive (dU/dr < 0)
        let f = lj_force_gpu(1.2, &params);
        assert!(f < 0.0);
    }

    #[test]
    fn test_lj_force_zero_beyond_cutoff() {
        let params = default_params();
        assert_eq!(lj_force_gpu(params.cutoff + 1.0, &params), 0.0);
    }

    #[test]
    fn test_pbc_distance_no_wrap() {
        let params = default_params();
        let a = [1.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let d = pbc_distance_gpu(a, b, params.box_size);
        assert!((d - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_pbc_distance_wrap() {
        let box_size = [10.0_f32; 3];
        let a = [0.5, 0.0, 0.0];
        let b = [9.5, 0.0, 0.0];
        // Minimum image: 1.0, not 9.0
        let d = pbc_distance_gpu(a, b, box_size);
        assert!((d - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_pbc_distance_self() {
        let box_size = [10.0_f32; 3];
        let a = [3.0, 4.0, 5.0];
        let d = pbc_distance_gpu(a, a, box_size);
        assert!(d < 1e-5);
    }

    #[test]
    fn test_compute_forces_gpu_newton3() {
        let mut buf = GpuMdBuffer::new(2);
        buf.atoms[0].pos = [0.0; 3];
        buf.atoms[1].pos = [1.2, 0.0, 0.0];
        buf.atoms[0].mass = 1.0;
        buf.atoms[1].mass = 1.0;
        let params = default_params();
        compute_forces_gpu(&mut buf, &params);
        // Newton's third law
        assert!((buf.atoms[0].force[0] + buf.atoms[1].force[0]).abs() < 1e-5);
    }

    #[test]
    fn test_compute_forces_gpu_zero_beyond_cutoff() {
        let mut buf = GpuMdBuffer::new(2);
        buf.atoms[0].pos = [0.0; 3];
        buf.atoms[1].pos = [5.0, 0.0, 0.0]; // beyond cutoff=3.5
        let params = default_params();
        compute_forces_gpu(&mut buf, &params);
        assert!(buf.atoms[0].force[0].abs() < 1e-8);
    }

    #[test]
    fn test_verlet_integrate_gpu_position() {
        let mut buf = GpuMdBuffer::new(1);
        buf.atoms[0].vel = [1.0, 0.0, 0.0];
        buf.atoms[0].force = [0.0; 3];
        verlet_integrate_gpu(&mut buf, 0.1);
        assert!((buf.atoms[0].pos[0] - 0.1).abs() < 1e-5);
    }

    #[test]
    fn test_verlet_integrate_gpu_step_counter() {
        let mut buf = GpuMdBuffer::new(1);
        verlet_integrate_gpu(&mut buf, 0.01);
        assert_eq!(buf.step, 1);
    }

    #[test]
    fn test_kinetic_energy_gpu_zero_vel() {
        let buf = make_buf_grid(4);
        let ke = kinetic_energy_gpu(&buf);
        assert!(ke.abs() < 1e-8);
    }

    #[test]
    fn test_kinetic_energy_gpu_nonzero() {
        let mut buf = GpuMdBuffer::new(2);
        buf.atoms[0].vel = [1.0, 0.0, 0.0];
        buf.atoms[0].mass = 2.0;
        buf.atoms[1].vel = [0.0, 0.0, 0.0];
        buf.atoms[1].mass = 1.0;
        let ke = kinetic_energy_gpu(&buf);
        // 0.5 * 2 * 1^2 = 1.0
        assert!((ke - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_potential_energy_gpu_empty() {
        let buf = GpuMdBuffer::new(0);
        let params = default_params();
        assert_eq!(potential_energy_gpu(&buf, &params), 0.0);
    }

    #[test]
    fn test_potential_energy_gpu_single() {
        let buf = GpuMdBuffer::new(1);
        let params = default_params();
        // Only one atom -> no pairs -> PE = 0
        assert_eq!(potential_energy_gpu(&buf, &params), 0.0);
    }

    #[test]
    fn test_temperature_gpu_zero_vel() {
        let buf = make_buf_grid(4);
        let t = temperature_gpu(&buf);
        assert!(t < 1e-6);
    }

    #[test]
    fn test_temperature_gpu_empty() {
        let buf = GpuMdBuffer::new(0);
        assert_eq!(temperature_gpu(&buf), 0.0);
    }

    #[test]
    fn test_temperature_gpu_nonzero() {
        let mut buf = GpuMdBuffer::new(3);
        for a in buf.atoms.iter_mut() {
            a.vel = [1.0, 1.0, 1.0];
            a.mass = 1.0;
        }
        let t = temperature_gpu(&buf);
        assert!(t > 0.0);
    }

    #[test]
    fn test_rescale_velocities_gpu() {
        let mut buf = GpuMdBuffer::new(4);
        for a in buf.atoms.iter_mut() {
            a.vel = [1.0, 0.5, 0.2];
            a.mass = 1.0;
        }
        let target = 300.0;
        rescale_velocities_gpu(&mut buf, target);
        let t_after = temperature_gpu(&buf);
        assert!((t_after - target).abs() < 1.0);
    }

    #[test]
    fn test_rescale_velocities_gpu_zero_vel_noop() {
        let mut buf = GpuMdBuffer::new(2);
        // All velocities zero -> temperature is zero -> no rescaling
        rescale_velocities_gpu(&mut buf, 300.0);
        for a in &buf.atoms {
            assert!(a.vel[0].abs() < 1e-8);
        }
    }

    #[test]
    fn test_buf_clone() {
        let buf = make_buf_grid(3);
        let buf2 = buf.clone();
        assert_eq!(buf2.len(), 3);
    }

    #[test]
    fn test_compute_forces_accumulate_many() {
        let mut buf = make_buf_grid(4);
        let params = default_params();
        compute_forces_gpu(&mut buf, &params);
        // Total force should be ~zero (Newton's 3rd law globally)
        let fx_total: f32 = buf.atoms.iter().map(|a| a.force[0]).sum();
        assert!(fx_total.abs() < 1e-4);
    }

    #[test]
    fn test_lj_potential_positive_repulsive() {
        let params = default_params();
        // Very short r -> strongly repulsive -> U > 0
        let u = lj_potential_gpu(0.5, &params);
        assert!(u > 0.0);
    }

    #[test]
    fn test_verlet_integrate_gpu_velocity_from_force() {
        let mut buf = GpuMdBuffer::new(1);
        buf.atoms[0].force = [2.0, 0.0, 0.0];
        buf.atoms[0].mass = 1.0;
        verlet_integrate_gpu(&mut buf, 0.5);
        // v = 0 + (2/1)*0.5 = 1.0
        assert!((buf.atoms[0].vel[0] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_pbc_distance_3d_wrap() {
        let box_size = [5.0_f32; 3];
        let a = [0.1, 0.1, 0.1];
        let b = [4.9, 4.9, 4.9];
        let d = pbc_distance_gpu(a, b, box_size);
        // Minimum image: delta = (-0.2, -0.2, -0.2) => r = 0.2*sqrt(3)
        let expected = 0.2 * 3.0_f32.sqrt();
        assert!((d - expected).abs() < 1e-4);
    }

    #[test]
    fn test_total_energy_two_atoms() {
        let mut buf = GpuMdBuffer::new(2);
        buf.atoms[0].pos = [0.0; 3];
        buf.atoms[0].mass = 1.0;
        buf.atoms[1].pos = [1.1, 0.0, 0.0]; // near LJ minimum
        buf.atoms[1].mass = 1.0;
        let params = default_params();
        let pe = potential_energy_gpu(&buf, &params);
        // Should be negative near minimum
        assert!(pe < 0.0);
    }
}
