// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Magnetohydrodynamics (MHD) SPH.
//!
//! Implements ideal and resistive MHD in SPH formulation:
//!
//! - [`MhdParticle`]: SPH particle carrying magnetic field B
//! - [`alfven_velocity`]: Alfvén wave speed v_A = |B|/√(μ₀ρ)
//! - [`magnetosonic_speed`]: Fast magnetosonic speed √(cs² + vA²)
//! - [`magnetic_pressure`]: Magnetic pressure p_mag = B²/(2μ₀)
//! - [`total_pressure`]: Total pressure including magnetic contribution
//! - [`LorentzForce`]: Lorentz force density **J** × **B**
//! - [`curl_b`]: SPH estimate of ∇ × **B**
//! - [`div_b`]: SPH estimate of ∇·**B** (divergence cleaning monitor)
//! - [`induction_equation_rhs`]: Ideal MHD induction: d**B**/dt = ∇×(**v**×**B**)
//! - [`dedner_cleaning`]: Dedner hyperbolic/parabolic (GLM) divergence cleaning
//!   with the per-particle potential [`MhdParticle::psi`]
//! - [`MhdSolver`]: Resistive MHD solver accumulating particles
//! - [`lundquist_number`]: Magnetic Reynolds number S = vA·L/η
//! - [`plasma_beta`]: Ratio of thermal to magnetic pressure

use std::f64::consts::PI;

use crate::kernel::SphKernel;

// ============================================================================
// Physical constants
// ============================================================================

/// Permeability of free space μ₀ (H/m).
const MU_0: f64 = 4.0 * PI * 1e-7;

// ============================================================================
// Math helpers
// ============================================================================

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

// ============================================================================
// MhdParticle
// ============================================================================

/// SPH particle carrying a magnetic field vector for MHD simulation.
#[derive(Debug, Clone)]
pub struct MhdParticle {
    /// Position **x** (m).
    pub pos: [f64; 3],
    /// Velocity **v** (m/s).
    pub vel: [f64; 3],
    /// Magnetic flux density **B** (T).
    pub b: [f64; 3],
    /// Mass density ρ (kg/m³).
    pub rho: f64,
    /// Thermal pressure p (Pa).
    pub p: f64,
    /// Particle mass m (kg).
    pub mass: f64,
    /// Smoothing length h (m).
    pub h: f64,
    /// Dedner divergence-cleaning potential ψ (T·m/s).
    ///
    /// Advanced by [`dedner_cleaning`]; zero for a freshly created particle.
    pub psi: f64,
}

impl MhdParticle {
    /// Create a new `MhdParticle` with all fields initialised.
    pub fn new(
        pos: [f64; 3],
        vel: [f64; 3],
        b: [f64; 3],
        rho: f64,
        p: f64,
        mass: f64,
        h: f64,
    ) -> Self {
        Self {
            pos,
            vel,
            b,
            rho,
            p,
            mass,
            h,
            psi: 0.0,
        }
    }

    /// Magnetic energy of this particle: E_mag = |B|² / (2·μ₀) · V,
    /// where V = mass / rho is the particle volume.
    pub fn magnetic_energy(&self) -> f64 {
        if self.rho.abs() < 1e-300 {
            return 0.0;
        }
        let vol = self.mass / self.rho;
        magnetic_pressure(self.b) * vol
    }

    /// Kinetic energy of this particle: E_kin = ½·m·|v|².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }
}

// ============================================================================
// Free functions
// ============================================================================

/// Alfvén wave speed: v_A = |B| / √(μ₀·ρ).
///
/// # Arguments
/// * `b`   – magnetic flux density (T)
/// * `rho` – mass density (kg/m³)
pub fn alfven_velocity(b: [f64; 3], rho: f64) -> f64 {
    let denom = MU_0 * rho;
    if denom <= 0.0 {
        return 0.0;
    }
    len3(b) / denom.sqrt()
}

/// Fast magnetosonic speed: c_ms = √(cs² + vA²).
///
/// # Arguments
/// * `cs` – sound speed (m/s)
/// * `va` – Alfvén speed (m/s)
pub fn magnetosonic_speed(cs: f64, va: f64) -> f64 {
    (cs * cs + va * va).sqrt()
}

/// Magnetic pressure: p_mag = |B|² / (2·μ₀)  (Pa).
///
/// # Arguments
/// * `b` – magnetic flux density (T)
pub fn magnetic_pressure(b: [f64; 3]) -> f64 {
    dot3(b, b) / (2.0 * MU_0)
}

/// Total pressure: p_tot = p + p_mag.
///
/// # Arguments
/// * `p` – thermal pressure (Pa)
/// * `b` – magnetic flux density (T)
pub fn total_pressure(p: f64, b: [f64; 3]) -> f64 {
    p + magnetic_pressure(b)
}

// ============================================================================
// LorentzForce
// ============================================================================

/// Lorentz force density **f** = **J** × **B** (N/m³).
pub struct LorentzForce {
    /// Current density **J** (A/m²).
    pub j: [f64; 3],
    /// Magnetic flux density **B** (T).
    pub b: [f64; 3],
}

impl LorentzForce {
    /// Create a `LorentzForce`.
    pub fn new(j: [f64; 3], b: [f64; 3]) -> Self {
        Self { j, b }
    }

    /// Compute force density **f** = **J** × **B** (N/m³).
    pub fn force_density(&self) -> [f64; 3] {
        cross3(self.j, self.b)
    }
}

// ============================================================================
// SPH operators
// ============================================================================

/// SPH estimate of ∇ × **B** at particle `i`, yielding current density proxy.
///
/// Uses the standard SPH curl approximation:
/// (∇ × **B**)_i ≈ Σ_j (m_j/ρ_j) **B**_j × ∇W(r_{ij}, h_i)
///
/// The kernel gradient is evaluated along the unit vector r̂_{ij}.
///
/// # Arguments
/// * `particles` – slice of MHD particles
/// * `i`         – index of the particle at which to evaluate
/// * `kernel`    – SPH kernel function `W(r, h) -> f64`
pub fn curl_b(particles: &[MhdParticle], i: usize, kernel: impl Fn(f64, f64) -> f64) -> [f64; 3] {
    let pi = &particles[i];
    let mut result = [0.0f64; 3];
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij = sub3(pi.pos, pj.pos);
        let r = len3(rij);
        if r < 1e-300 {
            continue;
        }
        // Finite-difference kernel gradient: dW/dr · r̂/r
        let dr = 1e-6 * pi.h.max(1e-10);
        let dw_dr = (kernel(r + dr, pi.h) - kernel(r - dr, pi.h)) / (2.0 * dr);
        let grad_w = scale3(rij, dw_dr / r);
        // contribution: (m_j / rho_j) * B_j × grad_W
        let factor = pj.mass / (pj.rho + 1e-300);
        let bj_cross_gw = cross3(pj.b, grad_w);
        result = add3(result, scale3(bj_cross_gw, factor));
    }
    result
}

/// SPH estimate of ∇·**B** at particle `i`.
///
/// ∇·**B** ≈ Σ_j (m_j/ρ_j) (**B**_j - **B**_i) · ∇W(r_{ij}, h_i)
///
/// Should be ~0 for a divergence-free field; used as a quality monitor.
///
/// # Arguments
/// * `particles` – slice of MHD particles
/// * `i`         – index of the particle at which to evaluate
/// * `kernel`    – SPH kernel function `W(r, h) -> f64`
pub fn div_b(particles: &[MhdParticle], i: usize, kernel: impl Fn(f64, f64) -> f64) -> f64 {
    let pi = &particles[i];
    let mut result = 0.0f64;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij = sub3(pi.pos, pj.pos);
        let r = len3(rij);
        if r < 1e-300 {
            continue;
        }
        let dr = 1e-6 * pi.h.max(1e-10);
        let dw_dr = (kernel(r + dr, pi.h) - kernel(r - dr, pi.h)) / (2.0 * dr);
        let grad_w = scale3(rij, dw_dr / r);
        let db = sub3(pj.b, pi.b);
        let factor = pj.mass / (pj.rho + 1e-300);
        result += factor * dot3(db, grad_w);
    }
    result
}

/// Right-hand side of the ideal MHD induction equation: d**B**/dt = ∇×(**v**×**B**).
///
/// Uses the identity ∇×(**v**×**B**) = (**B**·∇)**v** - (**v**·∇)**B** and an SPH
/// approximation analogous to the curl operator above.
///
/// # Arguments
/// * `particles` – slice of MHD particles
/// * `i`         – index of the particle
/// * `kernel`    – SPH kernel function
pub fn induction_equation_rhs(
    particles: &[MhdParticle],
    i: usize,
    kernel: impl Fn(f64, f64) -> f64,
) -> [f64; 3] {
    let pi = &particles[i];
    let mut result = [0.0f64; 3];
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij = sub3(pi.pos, pj.pos);
        let r = len3(rij);
        if r < 1e-300 {
            continue;
        }
        let dr = 1e-6 * pi.h.max(1e-10);
        let dw_dr = (kernel(r + dr, pi.h) - kernel(r - dr, pi.h)) / (2.0 * dr);
        let grad_w = scale3(rij, dw_dr / r);
        // vxB at particle j
        let vxb_j = cross3(pj.vel, pj.b);
        let factor = pj.mass / (pj.rho + 1e-300);
        // SPH curl contribution
        let contrib = cross3(grad_w, vxb_j);
        result = add3(result, scale3(contrib, -factor));
    }
    result
}

// ============================================================================
// Dedner cleaning
// ============================================================================

/// Analytic kernel gradient ∇_i W(**x**_i − **x**_j, h_i) of the pair
/// `(pi, pj)`, or `None` when the two particles coincide.
fn pair_kernel_gradient(
    pi: &MhdParticle,
    pj: &MhdParticle,
    kernel: &dyn SphKernel,
) -> Option<[f64; 3]> {
    let rij = sub3(pi.pos, pj.pos);
    let r = len3(rij);
    if r < 1e-300 {
        return None;
    }
    Some(scale3(rij, kernel.grad_w(r, pi.h) / r))
}

/// SPH divergence (∇·**B**)_i = Σ_j (m_j/ρ_j) (**B**_j − **B**_i) · ∇_i W_ij
/// with the analytic kernel gradient.
fn divergence_b_analytic(particles: &[MhdParticle], i: usize, kernel: &dyn SphKernel) -> f64 {
    let pi = &particles[i];
    let mut result = 0.0f64;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let Some(grad_w) = pair_kernel_gradient(pi, pj, kernel) else {
            continue;
        };
        let volume_j = pj.mass / (pj.rho + 1e-300);
        result += volume_j * dot3(sub3(pj.b, pi.b), grad_w);
    }
    result
}

/// SPH gradient (∇ψ)_i = Σ_j (m_j/ρ_j) (ψ_j − ψ_i) ∇_i W_ij of the stored
/// cleaning potentials.
fn gradient_psi(particles: &[MhdParticle], i: usize, kernel: &dyn SphKernel) -> [f64; 3] {
    let pi = &particles[i];
    let mut result = [0.0f64; 3];
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let Some(grad_w) = pair_kernel_gradient(pi, pj, kernel) else {
            continue;
        };
        let volume_j = pj.mass / (pj.rho + 1e-300);
        result = add3(result, scale3(grad_w, volume_j * (pj.psi - pi.psi)));
    }
    result
}

/// Dedner hyperbolic/parabolic (GLM) divergence cleaning step for **B**.
///
/// Each particle carries a cleaning potential ψ ([`MhdParticle::psi`]) that is
/// sourced by the local divergence error and feeds back on the field:
///
/// ```text
/// dψ/dt = −ch² ∇·B − ψ/τ,      τ_i = h_i / (σ ch)
/// dB/dt = −∇ψ
/// ```
///
/// so that ∇·**B** errors propagate away as waves of speed `ch` and are damped
/// on the time scale `τ` (Dedner et al. 2002; Tricco & Price 2012).  A
/// divergence-free field gives ψ = 0 and is left exactly unchanged; in
/// particular the magnetic energy of a solenoidal field is not drained.
///
/// One call advances the system by `dt` in four snapshot passes:
///
/// 1. ∇·**B** is evaluated for every particle from the current field;
/// 2. ψ is advanced with the damping term integrated exactly,
///    ψ ← ψ e^{−dt/τ} − ch² (∇·**B**) τ (1 − e^{−dt/τ}), which stays stable for
///    any `σ` and reduces to ψ ← ψ − ch² (∇·**B**) dt when `σ = 0`;
/// 3. ∇ψ is evaluated from the updated potentials;
/// 4. **B** ← **B** − dt ∇ψ.
///
/// Both operators use the analytic kernel gradient ∇_i W_ij with the smoothing
/// length h_i of the receiving particle.  The explicit scheme is stable only
/// for `ch · dt ≲ 0.3 · h`.  The cost is O(N²) per call.
///
/// # Arguments
/// * `particles` – mutable particle slice
/// * `ch`        – hyperbolic cleaning speed (m/s); `0` disables the cleaning
/// * `sigma`     – dimensionless damping factor σ; values ≤ 0 disable the damping
/// * `dt`        – time step (s)
/// * `kernel`    – SPH smoothing kernel
pub fn dedner_cleaning(
    particles: &mut [MhdParticle],
    ch: f64,
    sigma: f64,
    dt: f64,
    kernel: &dyn SphKernel,
) {
    if particles.is_empty() || !ch.is_finite() || !dt.is_finite() || ch.abs() < 1e-300 {
        return;
    }
    let n = particles.len();
    let ch_sq = ch * ch;

    let divergence: Vec<f64> = (0..n)
        .map(|i| divergence_b_analytic(particles, i, kernel))
        .collect();

    for (p, &div_i) in particles.iter_mut().zip(divergence.iter()) {
        let damping_rate = if sigma > 0.0 && p.h > 0.0 {
            sigma * ch.abs() / p.h
        } else {
            0.0
        };
        let x = damping_rate * dt;
        let (decay, source_time) = if x == 0.0 {
            (1.0, dt)
        } else {
            ((-x).exp(), -(-x).exp_m1() / x * dt)
        };
        p.psi = p.psi * decay - ch_sq * div_i * source_time;
    }

    let gradients: Vec<[f64; 3]> = (0..n).map(|i| gradient_psi(particles, i, kernel)).collect();
    for (p, grad) in particles.iter_mut().zip(gradients.iter()) {
        p.b = sub3(p.b, scale3(*grad, dt));
    }
}

// ============================================================================
// MhdSolver
// ============================================================================

/// Resistive MHD SPH solver accumulating particles.
pub struct MhdSolver {
    /// Collection of MHD particles.
    pub particles: Vec<MhdParticle>,
    /// Magnetic resistivity η (m²/s = Ω·m / μ₀).
    pub eta: f64,
}

impl MhdSolver {
    /// Create an empty `MhdSolver`.
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
            eta: 0.0,
        }
    }

    /// Create an `MhdSolver` with a specified resistivity η.
    pub fn with_resistivity(eta: f64) -> Self {
        Self {
            particles: Vec::new(),
            eta,
        }
    }

    /// Add a particle to the solver.
    pub fn add_particle(&mut self, p: MhdParticle) {
        self.particles.push(p);
    }

    /// Number of particles currently in the solver.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }

    /// Total magnetic energy: E_mag = Σ_i |B_i|² / (2μ₀) · V_i.
    pub fn total_magnetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.magnetic_energy()).sum()
    }

    /// Total kinetic energy: E_kin = Σ_i ½·m_i·|v_i|².
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
}

impl Default for MhdSolver {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Dimensionless parameters
// ============================================================================

/// Lundquist number (magnetic Reynolds number for Alfvén waves):
/// S = v_A · L / η.
///
/// # Arguments
/// * `va`  – Alfvén speed (m/s)
/// * `l`   – characteristic length L (m)
/// * `eta` – magnetic diffusivity η (m²/s)
pub fn lundquist_number(va: f64, l: f64, eta: f64) -> f64 {
    if eta.abs() < 1e-300 {
        return f64::INFINITY;
    }
    va * l / eta
}

/// Plasma beta: ratio of thermal to magnetic pressure β = p / p_mag.
///
/// β ≫ 1 → thermally dominated; β ≪ 1 → magnetically dominated.
///
/// # Arguments
/// * `p` – thermal pressure (Pa)
/// * `b` – magnetic flux density (T)
pub fn plasma_beta(p: f64, b: [f64; 3]) -> f64 {
    let pmag = magnetic_pressure(b);
    if pmag.abs() < 1e-300 {
        return f64::INFINITY;
    }
    p / pmag
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::CubicSplineKernel;

    // ---- MhdParticle --------------------------------------------------------

    #[test]
    fn test_mhd_particle_kinetic_energy() {
        let p = MhdParticle::new([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, 2.0, 0.05);
        assert!((p.kinetic_energy() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_mhd_particle_magnetic_energy_positive() {
        let b = [0.01, 0.0, 0.0];
        let p = MhdParticle::new([0.0; 3], [0.0; 3], b, 1000.0, 1e5, 0.001, 0.05);
        assert!(p.magnetic_energy() >= 0.0);
    }

    #[test]
    fn test_mhd_particle_zero_velocity_zero_kinetic() {
        let p = MhdParticle::new([0.0; 3], [0.0; 3], [0.0; 3], 1.0, 0.0, 1.0, 0.1);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    // ---- alfven_velocity ----------------------------------------------------

    #[test]
    fn test_alfven_velocity_positive() {
        let b = [0.01, 0.0, 0.0];
        let va = alfven_velocity(b, 1000.0);
        assert!(va > 0.0);
    }

    #[test]
    fn test_alfven_velocity_zero_rho_returns_zero() {
        let b = [0.01, 0.0, 0.0];
        assert_eq!(alfven_velocity(b, 0.0), 0.0);
    }

    #[test]
    fn test_alfven_velocity_zero_b_returns_zero() {
        assert_eq!(alfven_velocity([0.0; 3], 1000.0), 0.0);
    }

    #[test]
    fn test_alfven_velocity_scales_with_b() {
        let va1 = alfven_velocity([0.01, 0.0, 0.0], 1000.0);
        let va2 = alfven_velocity([0.02, 0.0, 0.0], 1000.0);
        assert!((va2 - 2.0 * va1).abs() < 1e-15);
    }

    #[test]
    fn test_alfven_velocity_decreases_with_density() {
        let va1 = alfven_velocity([0.01, 0.0, 0.0], 100.0);
        let va2 = alfven_velocity([0.01, 0.0, 0.0], 1000.0);
        assert!(va1 > va2);
    }

    // ---- magnetosonic_speed -------------------------------------------------

    #[test]
    fn test_magnetosonic_speed_greater_than_alfven() {
        let cs = 300.0;
        let va = 200.0;
        assert!(magnetosonic_speed(cs, va) > va);
    }

    #[test]
    fn test_magnetosonic_speed_greater_than_sound_speed() {
        let cs = 300.0;
        let va = 200.0;
        assert!(magnetosonic_speed(cs, va) > cs);
    }

    #[test]
    fn test_magnetosonic_speed_zero_va_equals_cs() {
        let cs = 300.0;
        assert!((magnetosonic_speed(cs, 0.0) - cs).abs() < 1e-10);
    }

    // ---- magnetic_pressure --------------------------------------------------

    #[test]
    fn test_magnetic_pressure_non_negative() {
        assert!(magnetic_pressure([0.01, 0.02, 0.03]) >= 0.0);
    }

    #[test]
    fn test_magnetic_pressure_zero_field() {
        assert_eq!(magnetic_pressure([0.0; 3]), 0.0);
    }

    #[test]
    fn test_magnetic_pressure_scales_with_b_squared() {
        let pm1 = magnetic_pressure([0.01, 0.0, 0.0]);
        let pm2 = magnetic_pressure([0.02, 0.0, 0.0]);
        assert!((pm2 - 4.0 * pm1).abs() < 1e-30);
    }

    // ---- total_pressure -----------------------------------------------------

    #[test]
    fn test_total_pressure_greater_than_thermal() {
        let pt = total_pressure(1e5, [0.01, 0.0, 0.0]);
        assert!(pt > 1e5);
    }

    #[test]
    fn test_total_pressure_zero_b_equals_thermal() {
        let pt = total_pressure(1e5, [0.0; 3]);
        assert!((pt - 1e5).abs() < 1e-10);
    }

    // ---- LorentzForce -------------------------------------------------------

    #[test]
    fn test_lorentz_force_j_cross_b() {
        // J along x, B along y → F along z
        let lf = LorentzForce::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let f = lf.force_density();
        assert!((f[2] - 1.0).abs() < 1e-15);
        assert!(f[0].abs() < 1e-15);
        assert!(f[1].abs() < 1e-15);
    }

    #[test]
    fn test_lorentz_force_parallel_zero() {
        let lf = LorentzForce::new([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let f = lf.force_density();
        for &fi in &f {
            assert!(fi.abs() < 1e-15);
        }
    }

    // ---- dedner_cleaning ----------------------------------------------------

    /// Cubic lattice of `n`³ particles (spacing `s`, smoothing length `h = s`)
    /// centred on the origin, carrying the field `field(pos)`.  With
    /// `rho = 1000` and `m = rho s³` every particle owns the volume `s³`.
    fn lattice(n: usize, s: f64, field: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<MhdParticle> {
        let rho = 1000.0;
        let mass = rho * s * s * s;
        let half = (n as f64 - 1.0) / 2.0;
        let mut particles = Vec::with_capacity(n * n * n);
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let pos = [
                        (ix as f64 - half) * s,
                        (iy as f64 - half) * s,
                        (iz as f64 - half) * s,
                    ];
                    particles.push(MhdParticle::new(
                        pos,
                        [0.0; 3],
                        field(pos),
                        rho,
                        1e5,
                        mass,
                        s,
                    ));
                }
            }
        }
        particles
    }

    /// Particles whose complete support (2 h) lies inside a 7³ lattice.
    fn is_interior(p: &MhdParticle, s: f64) -> bool {
        p.pos.iter().all(|c| c.abs() <= s + 1e-12)
    }

    fn rms_div_interior(particles: &[MhdParticle], s: f64) -> f64 {
        let kernel = CubicSplineKernel;
        let values: Vec<f64> = (0..particles.len())
            .filter(|&i| is_interior(&particles[i], s))
            .map(|i| div_b(particles, i, |r, h| kernel.w(r, h)))
            .collect();
        assert!(!values.is_empty());
        (values.iter().map(|d| d * d).sum::<f64>() / values.len() as f64).sqrt()
    }

    fn max_b_norm(particles: &[MhdParticle]) -> f64 {
        particles.iter().map(|p| len3(p.b)).fold(0.0, f64::max)
    }

    #[test]
    fn test_dedner_uniform_field_unchanged_after_many_calls() {
        // A divergence-free field must be left alone (the former implementation
        // multiplied B by a decay factor on every call).
        let b0 = [0.01, -0.02, 0.03];
        let s = 0.1;
        let mut particles = lattice(5, s, |_| b0);
        let ch = 5.0;
        let dt = 0.25 * s / ch;
        for _ in 0..100 {
            dedner_cleaning(&mut particles, ch, 1.0, dt, &CubicSplineKernel);
        }
        for p in &particles {
            for (component, expected) in p.b.iter().zip(b0.iter()) {
                assert!(
                    (component - expected).abs() < 1e-12,
                    "uniform B changed: {:?}",
                    p.b
                );
            }
            assert!(p.psi.abs() < 1e-12, "psi must stay zero, got {}", p.psi);
        }
    }

    #[test]
    fn test_dedner_divergence_error_sources_psi_with_correct_sign() {
        // B = (a x, 0, 0) has ∇·B = a > 0; with σ = 0 one step gives ψ = −ch² dt ∇·B.
        let s = 0.1;
        let a = 0.5;
        let mut particles = lattice(7, s, |pos| [a * pos[0], 0.0, 0.0]);
        let kernel = CubicSplineKernel;
        let div_before: Vec<f64> = (0..particles.len())
            .map(|i| div_b(&particles, i, |r, h| kernel.w(r, h)))
            .collect();
        let (ch, dt) = (4.0, 1e-4);
        dedner_cleaning(&mut particles, ch, 0.0, dt, &kernel);
        let mut checked = 0;
        for (i, p) in particles.iter().enumerate() {
            let expected = -ch * ch * dt * div_before[i];
            assert!(
                (p.psi - expected).abs() <= 1e-6 * expected.abs() + 1e-18,
                "particle {i}: psi {} vs −ch² dt div {}",
                p.psi,
                expected
            );
            if is_interior(p, s) {
                checked += 1;
                assert!(p.psi < 0.0, "positive divergence must give negative psi");
                // the lattice sum reproduces the continuum divergence a
                assert!(
                    (div_before[i] - a).abs() < 0.05 * a,
                    "SPH divergence {} should be close to {a}",
                    div_before[i]
                );
            }
        }
        assert_eq!(checked, 27);
    }

    #[test]
    fn test_dedner_psi_gradient_drives_b_against_gradient() {
        // ψ = g x with a uniform B (∇·B = 0): B ← B − dt ∇ψ lowers B_x.
        let s = 0.1;
        let g = 3.0;
        let b0 = [0.2, 0.0, 0.0];
        let mut particles = lattice(7, s, |_| b0);
        for p in particles.iter_mut() {
            p.psi = g * p.pos[0];
        }
        let before = particles.clone();
        let dt = 1e-3;
        let kernel = CubicSplineKernel;
        dedner_cleaning(&mut particles, 2.0, 0.0, dt, &kernel);

        for (i, p) in particles.iter().enumerate() {
            // independent evaluation of Σ_j V_j (ψ_j − ψ_i) ∇_i W_ij
            let mut grad = [0.0f64; 3];
            for (j, q) in before.iter().enumerate() {
                if j == i {
                    continue;
                }
                let d = sub3(before[i].pos, q.pos);
                let r = len3(d);
                let w = kernel.grad_w(r, before[i].h) / r;
                let volume = q.mass / q.rho;
                for k in 0..3 {
                    grad[k] += volume * (q.psi - before[i].psi) * w * d[k];
                }
            }
            for (k, g_k) in grad.iter().enumerate() {
                let expected = before[i].b[k] - dt * g_k;
                assert!(
                    (p.b[k] - expected).abs() < 1e-12,
                    "particle {i} component {k}: {} vs {expected}",
                    p.b[k]
                );
            }
            assert_eq!(p.psi, before[i].psi, "∇·B = 0 leaves ψ untouched at σ = 0");
            if is_interior(p, s) {
                let delta_bx = p.b[0] - b0[0];
                assert!(
                    (delta_bx + dt * g).abs() < 0.05 * dt * g,
                    "B_x must change by −dt g = {}, got {delta_bx}",
                    -dt * g
                );
            }
        }
    }

    #[test]
    fn test_dedner_damping_decays_psi_exactly() {
        // Uniform B and uniform ψ: no source, no gradient, only the damping
        // ψ ← ψ exp(−σ ch dt / h) acts.
        let s = 0.1;
        let (ch, sigma, dt) = (5.0, 0.8, 0.01);
        let b0 = [0.01, 0.0, 0.0];
        let mut particles = lattice(3, s, |_| b0);
        for p in particles.iter_mut() {
            p.psi = 2.0;
        }
        let factor = (-sigma * ch * dt / s).exp();
        for step in 1..=5 {
            dedner_cleaning(&mut particles, ch, sigma, dt, &CubicSplineKernel);
            let expected = 2.0 * factor.powi(step);
            for p in &particles {
                assert!(
                    (p.psi - expected).abs() < 1e-13 * expected,
                    "step {step}: psi {} vs {expected}",
                    p.psi
                );
                assert_eq!(p.b, b0);
            }
        }
    }

    #[test]
    fn test_dedner_damping_is_unconditionally_stable() {
        // σ ch dt / h = 50 would flip the sign of ψ with a forward-Euler
        // damping term; the exponential integrator keeps it positive and tiny.
        let s = 0.1;
        let mut particles = lattice(3, s, |_| [0.0, 0.0, 0.01]);
        for p in particles.iter_mut() {
            p.psi = 1.0;
        }
        dedner_cleaning(&mut particles, 5.0, 1.0, 1.0, &CubicSplineKernel);
        for p in &particles {
            assert!(p.psi > 0.0 && p.psi < 1e-20, "psi = {}", p.psi);
        }
    }

    #[test]
    fn test_dedner_solenoidal_rotation_field_untouched_in_interior() {
        // B = g (−y, x, 0) is divergence free; on interior particles the
        // lattice sum is exact, so one step leaves B and ψ unchanged.
        let s = 0.1;
        let g = 0.4;
        let mut particles = lattice(7, s, |pos| [-g * pos[1], g * pos[0], 0.0]);
        let before: Vec<[f64; 3]> = particles.iter().map(|p| p.b).collect();
        dedner_cleaning(&mut particles, 4.0, 1.0, 1e-3, &CubicSplineKernel);
        let mut checked = 0;
        for (p, b_before) in particles.iter().zip(before.iter()) {
            if !is_interior(p, s) {
                continue;
            }
            checked += 1;
            assert!(p.psi.abs() < 1e-9, "psi = {}", p.psi);
            for (component, expected) in p.b.iter().zip(b_before.iter()) {
                assert!((component - expected).abs() < 1e-12);
            }
        }
        assert_eq!(checked, 27);
    }

    #[test]
    fn test_dedner_cleans_monopole_perturbation() {
        // Smooth monopole bump B = B0 + b_a (r/s) exp(−r²/w²): ∇·B > 0 around
        // the origin.  Cleaning must remove the divergence error without
        // blowing up the field.
        let s = 0.1;
        let w = 1.5 * s;
        let b0 = [0.003, -0.002, 0.001];
        let b_a = 0.01;
        let mut particles = lattice(7, s, |pos| {
            let envelope = (-(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]) / (w * w)).exp();
            [
                b0[0] + b_a * pos[0] / s * envelope,
                b0[1] + b_a * pos[1] / s * envelope,
                b0[2] + b_a * pos[2] / s * envelope,
            ]
        });
        let ch = 5.0;
        let dt = 0.25 * s / ch;
        let rms_initial = rms_div_interior(&particles, s);
        let b_max_initial = max_b_norm(&particles);
        assert!(
            rms_initial > 0.05,
            "seeded divergence too small: {rms_initial}"
        );

        for _ in 0..10 {
            dedner_cleaning(&mut particles, ch, 1.0, dt, &CubicSplineKernel);
        }
        let rms_mid = rms_div_interior(&particles, s);
        assert!(
            rms_mid < 0.5 * rms_initial,
            "rms ∇·B after 10 steps {rms_mid} vs initial {rms_initial}"
        );

        for _ in 0..40 {
            dedner_cleaning(&mut particles, ch, 1.0, dt, &CubicSplineKernel);
        }
        let rms_final = rms_div_interior(&particles, s);
        assert!(
            rms_final < 0.05 * rms_initial,
            "rms ∇·B after 50 steps {rms_final} vs initial {rms_initial}"
        );
        assert!(
            max_b_norm(&particles) < 2.0 * b_max_initial,
            "cleaning must not amplify the field"
        );
        assert!(particles.iter().all(|p| p.b.iter().all(|c| c.is_finite())));
    }

    #[test]
    fn test_dedner_disabled_inputs_are_no_ops() {
        let s = 0.1;
        let field = |pos: [f64; 3]| [0.5 * pos[0], 0.0, 0.0];
        let reference = lattice(3, s, field);
        for (ch, dt) in [(0.0, 1e-3), (5.0, 0.0), (f64::NAN, 1e-3), (5.0, f64::NAN)] {
            let mut particles = lattice(3, s, field);
            dedner_cleaning(&mut particles, ch, 1.0, dt, &CubicSplineKernel);
            for (p, r) in particles.iter().zip(reference.iter()) {
                assert_eq!(p.b, r.b, "ch = {ch}, dt = {dt}");
                assert_eq!(p.psi, 0.0, "ch = {ch}, dt = {dt}");
            }
        }
        let mut empty: Vec<MhdParticle> = Vec::new();
        dedner_cleaning(&mut empty, 5.0, 1.0, 1e-3, &CubicSplineKernel);
        assert!(empty.is_empty());
    }

    #[test]
    fn test_dedner_nonpositive_sigma_disables_damping() {
        let s = 0.1;
        let field = |pos: [f64; 3]| [0.5 * pos[0], 0.0, 0.0];
        let mut zero = lattice(3, s, field);
        let mut negative = lattice(3, s, field);
        dedner_cleaning(&mut zero, 5.0, 0.0, 1e-3, &CubicSplineKernel);
        dedner_cleaning(&mut negative, 5.0, -2.0, 1e-3, &CubicSplineKernel);
        for (a, b) in zero.iter().zip(negative.iter()) {
            assert_eq!(a.psi, b.psi);
            assert_eq!(a.b, b.b);
        }
    }

    #[test]
    fn test_mhd_particle_new_starts_with_zero_psi() {
        let p = MhdParticle::new([0.0; 3], [0.0; 3], [0.01, 0.0, 0.0], 1.0, 0.0, 1.0, 0.1);
        assert_eq!(p.psi, 0.0);
    }

    // ---- MhdSolver ----------------------------------------------------------

    #[test]
    fn test_mhd_solver_empty() {
        let solver = MhdSolver::new();
        assert_eq!(solver.particle_count(), 0);
        assert_eq!(solver.total_magnetic_energy(), 0.0);
        assert_eq!(solver.total_kinetic_energy(), 0.0);
    }

    #[test]
    fn test_mhd_solver_add_particle() {
        let mut solver = MhdSolver::new();
        solver.add_particle(MhdParticle::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            0.0,
            1.0,
            0.1,
        ));
        assert_eq!(solver.particle_count(), 1);
    }

    #[test]
    fn test_mhd_solver_total_energy_non_negative() {
        let mut solver = MhdSolver::new();
        solver.add_particle(MhdParticle::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.01, 0.0, 0.0],
            1000.0,
            1e5,
            0.001,
            0.05,
        ));
        assert!(solver.total_magnetic_energy() >= 0.0);
        assert!(solver.total_kinetic_energy() >= 0.0);
    }

    // ---- lundquist_number ---------------------------------------------------

    #[test]
    fn test_lundquist_number_positive() {
        let s = lundquist_number(1e4, 1.0, 1.0);
        assert!(s > 0.0);
    }

    #[test]
    fn test_lundquist_number_zero_eta_infinite() {
        let s = lundquist_number(1e4, 1.0, 0.0);
        assert!(s.is_infinite());
    }

    // ---- plasma_beta --------------------------------------------------------

    #[test]
    fn test_plasma_beta_large_when_b_small() {
        let beta = plasma_beta(1e5, [1e-10, 0.0, 0.0]);
        assert!(beta > 1.0, "Beta must be large when B is small, got {beta}");
    }

    #[test]
    fn test_plasma_beta_small_when_b_large() {
        let beta = plasma_beta(1.0, [10.0, 0.0, 0.0]);
        assert!(beta < 1.0, "Beta must be small when B is large, got {beta}");
    }

    #[test]
    fn test_plasma_beta_non_negative() {
        let beta = plasma_beta(1e5, [0.01, 0.02, 0.0]);
        assert!(beta >= 0.0);
    }

    #[test]
    fn test_plasma_beta_zero_b_infinite() {
        let beta = plasma_beta(1e5, [0.0; 3]);
        assert!(beta.is_infinite());
    }

    // ---- SPH operators (simple two-particle sanity checks) ------------------

    fn simple_kernel(r: f64, h: f64) -> f64 {
        if h < 1e-300 {
            return 0.0;
        }
        let q = r / h;
        if q < 1.0 { 1.0 - q } else { 0.0 }
    }

    #[test]
    fn test_div_b_zero_field_gives_zero() {
        let particles = vec![
            MhdParticle::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 0.0, 0.001, 0.5),
            MhdParticle::new([0.1, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 0.0, 0.001, 0.5),
        ];
        let db = div_b(&particles, 0, simple_kernel);
        assert!(db.abs() < 1e-10);
    }

    #[test]
    fn test_curl_b_zero_field_gives_zero() {
        let particles = vec![
            MhdParticle::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 0.0, 0.001, 0.5),
            MhdParticle::new([0.1, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 0.0, 0.001, 0.5),
        ];
        let curl = curl_b(&particles, 0, simple_kernel);
        for &c in &curl {
            assert!(c.abs() < 1e-10);
        }
    }

    #[test]
    fn test_induction_zero_vel_gives_zero() {
        let particles = vec![
            MhdParticle::new(
                [0.0, 0.0, 0.0],
                [0.0; 3],
                [0.01, 0.0, 0.0],
                1.0,
                0.0,
                0.001,
                0.5,
            ),
            MhdParticle::new(
                [0.1, 0.0, 0.0],
                [0.0; 3],
                [0.01, 0.0, 0.0],
                1.0,
                0.0,
                0.001,
                0.5,
            ),
        ];
        let rhs = induction_equation_rhs(&particles, 0, simple_kernel);
        for &r in &rhs {
            assert!(r.abs() < 1e-10);
        }
    }

    #[test]
    fn test_pi_constant_sanity() {
        // Ensure PI is available
        assert!((PI - std::f64::consts::PI).abs() < 1e-15);
    }
}
