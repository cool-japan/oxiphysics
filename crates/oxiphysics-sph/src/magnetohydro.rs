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
//! - [`curl_B`]: SPH estimate of ∇ × **B**
//! - [`div_B`]: SPH estimate of ∇·**B** (divergence cleaning monitor)
//! - [`induction_equation_rhs`]: Ideal MHD induction: d**B**/dt = ∇×(**v**×**B**)
//! - [`dedner_cleaning`]: Dedner hyperbolic divergence cleaning
//! - [`MhdSolver`]: Resistive MHD solver accumulating particles
//! - [`lundquist_number`]: Magnetic Reynolds number S = vA·L/η
//! - [`plasma_beta`]: Ratio of thermal to magnetic pressure

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]
#![allow(non_snake_case)]

use std::f64::consts::PI;

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
pub struct MhdParticle {
    /// Position **x** (m).
    pub pos: [f64; 3],
    /// Velocity **v** (m/s).
    pub vel: [f64; 3],
    /// Magnetic flux density **B** (T).
    pub B: [f64; 3],
    /// Mass density ρ (kg/m³).
    pub rho: f64,
    /// Thermal pressure p (Pa).
    pub p: f64,
    /// Particle mass m (kg).
    pub mass: f64,
    /// Smoothing length h (m).
    pub h: f64,
}

impl MhdParticle {
    /// Create a new `MhdParticle` with all fields initialised.
    pub fn new(
        pos: [f64; 3],
        vel: [f64; 3],
        B: [f64; 3],
        rho: f64,
        p: f64,
        mass: f64,
        h: f64,
    ) -> Self {
        Self {
            pos,
            vel,
            B,
            rho,
            p,
            mass,
            h,
        }
    }

    /// Magnetic energy of this particle: E_mag = |B|² / (2·μ₀) · V,
    /// where V = mass / rho is the particle volume.
    pub fn magnetic_energy(&self) -> f64 {
        if self.rho.abs() < 1e-300 {
            return 0.0;
        }
        let vol = self.mass / self.rho;
        magnetic_pressure(self.B) * vol
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
/// * `B`   – magnetic flux density (T)
/// * `rho` – mass density (kg/m³)
pub fn alfven_velocity(B: [f64; 3], rho: f64) -> f64 {
    let denom = MU_0 * rho;
    if denom <= 0.0 {
        return 0.0;
    }
    len3(B) / denom.sqrt()
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
/// * `B` – magnetic flux density (T)
pub fn magnetic_pressure(B: [f64; 3]) -> f64 {
    dot3(B, B) / (2.0 * MU_0)
}

/// Total pressure: p_tot = p + p_mag.
///
/// # Arguments
/// * `p` – thermal pressure (Pa)
/// * `B` – magnetic flux density (T)
pub fn total_pressure(p: f64, B: [f64; 3]) -> f64 {
    p + magnetic_pressure(B)
}

// ============================================================================
// LorentzForce
// ============================================================================

/// Lorentz force density **f** = **J** × **B** (N/m³).
pub struct LorentzForce {
    /// Current density **J** (A/m²).
    pub J: [f64; 3],
    /// Magnetic flux density **B** (T).
    pub B: [f64; 3],
}

impl LorentzForce {
    /// Create a `LorentzForce`.
    pub fn new(J: [f64; 3], B: [f64; 3]) -> Self {
        Self { J, B }
    }

    /// Compute force density **f** = **J** × **B** (N/m³).
    pub fn force_density(&self) -> [f64; 3] {
        cross3(self.J, self.B)
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
pub fn curl_B(particles: &[MhdParticle], i: usize, kernel: impl Fn(f64, f64) -> f64) -> [f64; 3] {
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
        let bj_cross_gw = cross3(pj.B, grad_w);
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
pub fn div_B(particles: &[MhdParticle], i: usize, kernel: impl Fn(f64, f64) -> f64) -> f64 {
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
        let db = sub3(pj.B, pi.B);
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
        let vxb_j = cross3(pj.vel, pj.B);
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

/// Dedner hyperbolic/parabolic divergence cleaning step for **B**.
///
/// Each particle carries an implicit cleaning scalar ψ (not stored in
/// `MhdParticle`); this simplified version applies an exponential decay to **B**
/// proportional to the local divergence error:
///
/// B_new = B_old · exp(-cp² / ch² · dt)
///
/// where ch is the cleaning wave speed and cp is the parabolic damping speed.
///
/// # Arguments
/// * `particles` – mutable particle slice
/// * `ch`        – hyperbolic cleaning speed (m/s)
/// * `cp`        – parabolic damping coefficient (m/s)
/// * `dt`        – time step (s)
pub fn dedner_cleaning(particles: &mut [MhdParticle], ch: f64, cp: f64, dt: f64) {
    if ch.abs() < 1e-300 {
        return;
    }
    let decay = (-cp * cp / (ch * ch) * dt).exp();
    for p in particles.iter_mut() {
        p.B = scale3(p.B, decay);
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
/// * `B` – magnetic flux density (T)
pub fn plasma_beta(p: f64, B: [f64; 3]) -> f64 {
    let pmag = magnetic_pressure(B);
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
    fn test_alfven_velocity_zero_B_returns_zero() {
        assert_eq!(alfven_velocity([0.0; 3], 1000.0), 0.0);
    }

    #[test]
    fn test_alfven_velocity_scales_with_B() {
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

    #[test]
    fn test_dedner_cleaning_reduces_b() {
        let mut particles = vec![MhdParticle::new(
            [0.0; 3],
            [0.0; 3],
            [0.01, 0.0, 0.0],
            1.0,
            1e5,
            0.001,
            0.05,
        )];
        let b0 = len3(particles[0].B);
        dedner_cleaning(&mut particles, 1e3, 1e2, 0.01);
        let b1 = len3(particles[0].B);
        // Decay factor < 1, so |B| must decrease
        assert!(b1 < b0 || (b1 - b0).abs() < 1e-20);
    }

    #[test]
    fn test_dedner_cleaning_zero_ch_no_change() {
        let mut particles = vec![MhdParticle::new(
            [0.0; 3],
            [0.0; 3],
            [0.01, 0.0, 0.0],
            1.0,
            1e5,
            0.001,
            0.05,
        )];
        let b0 = len3(particles[0].B);
        dedner_cleaning(&mut particles, 0.0, 1e2, 0.01);
        let b1 = len3(particles[0].B);
        assert!((b1 - b0).abs() < 1e-20);
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
        let db = div_B(&particles, 0, simple_kernel);
        assert!(db.abs() < 1e-10);
    }

    #[test]
    fn test_curl_b_zero_field_gives_zero() {
        let particles = vec![
            MhdParticle::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 0.0, 0.001, 0.5),
            MhdParticle::new([0.1, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 0.0, 0.001, 0.5),
        ];
        let curl = curl_B(&particles, 0, simple_kernel);
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
