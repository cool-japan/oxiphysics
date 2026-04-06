#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electromagnetic FEM — eddy currents and Maxwell's equations.
//!
//! Provides tools for eddy-current analysis, transmission lines and inductors:
//!
//! - [`ElectromagneticFEM`]: EM material properties container
//! - [`skin_depth`]: Classical skin depth formula
//! - [`impedance_surface`]: Surface impedance (R_s, X_s)
//! - [`EddyCurrent`]: Eddy-current model with power loss and EMF
//! - [`helmholtz_1d_solution`]: 1-D vector potential in a conductor
//! - [`CoaxialLine`]: Coaxial transmission line parameters
//! - [`SolenoidInductor`]: Solenoid inductance and field
//! - [`mutual_inductance_coaxial`]: Neumann formula approximation
//! - [`maxwell_stress_tensor`]: Maxwell stress tensor T_{ij}
//! - [`poynting_vector`]: Poynting vector S = (1/μ₀) E × B
//! - [`em_energy_density`]: Electromagnetic energy density
//! - [`waveguide_te10_cutoff`]: TE10 cutoff frequency

#![allow(dead_code)]
#![allow(non_snake_case)]

use std::f64::consts::PI;

// ============================================================================
// Physical constants
// ============================================================================

/// Permeability of free space μ₀ (H/m).
const MU_0: f64 = 4.0 * PI * 1e-7;
/// Permittivity of free space ε₀ (F/m).
const EPS_0: f64 = 8.854_187_817e-12;
/// Speed of light in vacuum c (m/s).
const C_LIGHT: f64 = 2.997_924_58e8;

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
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

// ============================================================================
// ElectromagneticFEM
// ============================================================================

/// Material and frequency parameters for an electromagnetic FEM problem.
pub struct ElectromagneticFEM {
    /// Number of nodes in the FEM mesh.
    pub n_nodes: usize,
    /// Electrical conductivity σ (S/m).
    pub sigma: f64,
    /// Relative permeability μ_r (dimensionless).
    pub mu_r: f64,
    /// Relative permittivity ε_r (dimensionless).
    pub epsilon_r: f64,
    /// Angular frequency ω = 2πf (rad/s).
    pub omega: f64,
}

impl ElectromagneticFEM {
    /// Create a new `ElectromagneticFEM`.
    pub fn new(n_nodes: usize, sigma: f64, mu_r: f64, epsilon_r: f64, omega: f64) -> Self {
        Self {
            n_nodes,
            sigma,
            mu_r,
            epsilon_r,
            omega,
        }
    }

    /// Skin depth δ for this material and frequency.
    pub fn skin_depth(&self) -> f64 {
        skin_depth(self.omega, self.sigma, self.mu_r)
    }

    /// Surface impedance (R_s, X_s) for this material and frequency.
    pub fn surface_impedance(&self) -> (f64, f64) {
        impedance_surface(self.omega, self.sigma, self.mu_r)
    }
}

// ============================================================================
// Free functions
// ============================================================================

/// Skin depth δ = sqrt(2 / (ω·σ·μ_r·μ₀)).
///
/// # Arguments
/// * `omega` – angular frequency ω (rad/s)
/// * `sigma` – electrical conductivity σ (S/m)
/// * `mu_r`  – relative permeability μ_r
pub fn skin_depth(omega: f64, sigma: f64, mu_r: f64) -> f64 {
    let denom = omega * sigma * mu_r * MU_0;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    (2.0 / denom).sqrt()
}

/// Surface impedance Z_s = (1 + i)·R_s where R_s = 1/(σ·δ).
///
/// Returns `(R_s, X_s)` where X_s = R_s (for good conductors).
///
/// # Arguments
/// * `omega` – angular frequency ω (rad/s)
/// * `sigma` – electrical conductivity σ (S/m)
/// * `mu_r`  – relative permeability μ_r
pub fn impedance_surface(omega: f64, sigma: f64, mu_r: f64) -> (f64, f64) {
    let delta = skin_depth(omega, sigma, mu_r);
    if delta.is_infinite() || sigma.abs() < 1e-300 {
        return (0.0, 0.0);
    }
    let rs = 1.0 / (sigma * delta);
    (rs, rs)
}

// ============================================================================
// EddyCurrent
// ============================================================================

/// Eddy current model for a conductive loop or bulk conductor.
pub struct EddyCurrent {
    /// Driving frequency f (Hz).
    pub frequency: f64,
    /// Electrical conductivity σ (S/m).
    pub conductivity: f64,
    /// Magnetic permeability μ = μ_r·μ₀ (H/m).
    pub permeability: f64,
    /// Applied magnetic field amplitude B₀ (T).
    pub field_amplitude: f64,
}

impl EddyCurrent {
    /// Create a new `EddyCurrent`.
    pub fn new(frequency: f64, conductivity: f64, permeability: f64, field_amplitude: f64) -> Self {
        Self {
            frequency,
            conductivity,
            permeability,
            field_amplitude,
        }
    }

    /// Time-averaged power loss per unit volume P = σ·|E|²/2 ≈ σ·(ω·B₀·δ)²/4.
    pub fn power_loss_per_volume(&self) -> f64 {
        let omega = 2.0 * PI * self.frequency;
        let delta = skin_depth(omega, self.conductivity, self.permeability / MU_0);
        // Induced electric field amplitude inside conductor E ≈ ω·B₀·δ
        let e_amp = omega * self.field_amplitude * delta;
        0.5 * self.conductivity * e_amp * e_amp
    }

    /// Induced EMF in a loop of area A: |ε| = ω·B₀·A.
    pub fn induced_emf(&self, area: f64) -> f64 {
        let omega = 2.0 * PI * self.frequency;
        omega * self.field_amplitude * area
    }

    /// Quality factor Q = ω·μ / (2·σ·Rs²) (simplified figure of merit).
    pub fn quality_factor(&self) -> f64 {
        let omega = 2.0 * PI * self.frequency;
        let delta = skin_depth(omega, self.conductivity, self.permeability / MU_0);
        if delta.is_infinite() || self.conductivity.abs() < 1e-300 {
            return 0.0;
        }
        let rs = 1.0 / (self.conductivity * delta);
        if rs.abs() < 1e-300 {
            return f64::INFINITY;
        }
        omega * self.permeability / (2.0 * rs * rs * self.conductivity)
    }
}

// ============================================================================
// Helmholtz 1-D solution
// ============================================================================

/// 1-D Helmholtz solution for vector potential in a semi-infinite conductor:
///
/// A(x) = A₀·exp(-(1+i)·x/δ)
///
/// Returns `(Re[A], Im[A])`.
///
/// # Arguments
/// * `x`     – depth into conductor (m, x ≥ 0)
/// * `delta` – skin depth δ (m)
/// * `a0`    – surface amplitude A₀ (T·m)
pub fn helmholtz_1d_solution(x: f64, delta: f64, a0: f64) -> (f64, f64) {
    if delta.abs() < 1e-300 {
        return (0.0, 0.0);
    }
    let exp_decay = (-x / delta).exp();
    let phase = -x / delta;
    (a0 * exp_decay * phase.cos(), a0 * exp_decay * phase.sin())
}

// ============================================================================
// CoaxialLine
// ============================================================================

/// Coaxial transmission line with inner and outer conductors.
pub struct CoaxialLine {
    /// Inner conductor radius a (m).
    pub inner_r: f64,
    /// Outer conductor inner radius b (m).
    pub outer_r: f64,
    /// Conductor conductivity σ (S/m).
    pub sigma: f64,
    /// Relative permittivity ε_r of the dielectric.
    pub epsilon_r: f64,
    /// Relative permeability μ_r of the dielectric.
    pub mu_r: f64,
}

impl CoaxialLine {
    /// Create a `CoaxialLine`.
    pub fn new(inner_r: f64, outer_r: f64, sigma: f64, epsilon_r: f64, mu_r: f64) -> Self {
        Self {
            inner_r,
            outer_r,
            sigma,
            epsilon_r,
            mu_r,
        }
    }

    /// Capacitance per unit length C = 2π·ε₀·ε_r / ln(b/a)  (F/m).
    pub fn capacitance_per_length(&self) -> f64 {
        if self.inner_r <= 0.0 || self.outer_r <= self.inner_r {
            return 0.0;
        }
        2.0 * PI * EPS_0 * self.epsilon_r / (self.outer_r / self.inner_r).ln()
    }

    /// Inductance per unit length L = μ₀·μ_r / (2π) · ln(b/a)  (H/m).
    pub fn inductance_per_length(&self) -> f64 {
        if self.inner_r <= 0.0 || self.outer_r <= self.inner_r {
            return 0.0;
        }
        MU_0 * self.mu_r / (2.0 * PI) * (self.outer_r / self.inner_r).ln()
    }

    /// Characteristic impedance Z₀ = sqrt(L/C) = (η / 2π) · ln(b/a)  (Ω).
    pub fn impedance(&self) -> f64 {
        let l = self.inductance_per_length();
        let c = self.capacitance_per_length();
        if c.abs() < 1e-300 {
            return f64::INFINITY;
        }
        (l / c).sqrt()
    }

    /// Velocity factor v/c = 1 / sqrt(ε_r · μ_r).
    pub fn velocity_factor(&self) -> f64 {
        (self.epsilon_r * self.mu_r).sqrt().recip()
    }
}

// ============================================================================
// SolenoidInductor
// ============================================================================

/// Air-core solenoid inductor.
pub struct SolenoidInductor {
    /// Number of turns N.
    pub n_turns: usize,
    /// Solenoid length l (m).
    pub length: f64,
    /// Solenoid radius r (m).
    pub radius: f64,
}

impl SolenoidInductor {
    /// Create a `SolenoidInductor`.
    pub fn new(n_turns: usize, length: f64, radius: f64) -> Self {
        Self {
            n_turns,
            length,
            radius,
        }
    }

    /// Inductance L = μ₀·N²·A / l  (H).
    pub fn inductance(&self) -> f64 {
        if self.length.abs() < 1e-300 {
            return 0.0;
        }
        let area = PI * self.radius * self.radius;
        let n = self.n_turns as f64;
        MU_0 * n * n * area / self.length
    }

    /// Axial magnetic field at the solenoid centre: B = μ₀·N·I / l  (T).
    pub fn field_at_center(&self, current: f64) -> f64 {
        if self.length.abs() < 1e-300 {
            return 0.0;
        }
        let n = self.n_turns as f64;
        MU_0 * n * current / self.length
    }
}

// ============================================================================
// mutual_inductance_coaxial
// ============================================================================

/// Approximate mutual inductance between two coaxial circular loops.
///
/// Uses the Neumann approximation valid when the loops are well separated (dist ≫ r):
///
/// M ≈ μ₀·π·r₁²·r₂² / (2·(dist² + ((r₁+r₂)/2)²)^{3/2})
///
/// # Arguments
/// * `r1`   – radius of loop 1 (m)
/// * `r2`   – radius of loop 2 (m)
/// * `dist` – separation distance (m)
pub fn mutual_inductance_coaxial(r1: f64, r2: f64, dist: f64) -> f64 {
    let avg_r = (r1 + r2) / 2.0;
    let denom = (dist * dist + avg_r * avg_r).powf(1.5);
    if denom < 1e-300 {
        return 0.0;
    }
    MU_0 * PI * r1 * r1 * r2 * r2 / (2.0 * denom)
}

// ============================================================================
// maxwell_stress_tensor
// ============================================================================

/// Maxwell stress tensor T_{ij} = ε₀(E_i·E_j - δ_{ij}|E|²/2) + (1/μ₀)(B_i·B_j - δ_{ij}|B|²/2).
///
/// # Arguments
/// * `e` – electric field vector (V/m)
/// * `b` – magnetic flux density vector (T)
pub fn maxwell_stress_tensor(e: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 3] {
    let e2 = dot3(e, e);
    let b2 = dot3(b, b);
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            t[i][j] =
                EPS_0 * (e[i] * e[j] - 0.5 * delta * e2) + (b[i] * b[j] - 0.5 * delta * b2) / MU_0;
        }
    }
    t
}

// ============================================================================
// poynting_vector
// ============================================================================

/// Poynting vector **S** = (1/μ₀) **E** × **B**  (W/m²).
///
/// # Arguments
/// * `e` – electric field (V/m)
/// * `b` – magnetic flux density (T)
pub fn poynting_vector(e: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let exb = cross3(e, b);
    [exb[0] / MU_0, exb[1] / MU_0, exb[2] / MU_0]
}

// ============================================================================
// em_energy_density
// ============================================================================

/// Electromagnetic energy density u = ε₀|E|²/2 + |B|²/(2μ₀)  (J/m³).
///
/// # Arguments
/// * `e` – electric field (V/m)
/// * `b` – magnetic flux density (T)
pub fn em_energy_density(e: [f64; 3], b: [f64; 3]) -> f64 {
    0.5 * EPS_0 * dot3(e, e) + dot3(b, b) / (2.0 * MU_0)
}

// ============================================================================
// waveguide_te10_cutoff
// ============================================================================

/// TE₁₀ cutoff frequency for a rectangular waveguide of width a:
///
/// f_c = c / (2a)
///
/// # Arguments
/// * `a` – broad-wall dimension of the waveguide (m)
pub fn waveguide_te10_cutoff(a: f64) -> f64 {
    if a.abs() < 1e-300 {
        return f64::INFINITY;
    }
    C_LIGHT / (2.0 * a)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ---- skin_depth ---------------------------------------------------------

    #[test]
    fn test_skin_depth_positive() {
        let delta = skin_depth(2.0 * PI * 50.0, 5.8e7, 1.0);
        assert!(delta > 0.0);
    }

    #[test]
    fn test_skin_depth_decreases_with_frequency() {
        let delta1 = skin_depth(2.0 * PI * 50.0, 5.8e7, 1.0);
        let delta2 = skin_depth(2.0 * PI * 1e6, 5.8e7, 1.0);
        assert!(delta2 < delta1, "Skin depth must decrease with frequency");
    }

    #[test]
    fn test_skin_depth_zero_omega_infinite() {
        let delta = skin_depth(0.0, 5.8e7, 1.0);
        assert!(delta.is_infinite());
    }

    #[test]
    fn test_skin_depth_known_value_copper_50hz() {
        // δ ≈ 9.35 mm for copper at 50 Hz
        let delta = skin_depth(2.0 * PI * 50.0, 5.8e7, 1.0);
        assert!(
            (delta - 9.35e-3).abs() < 0.5e-3,
            "Expected ~9.35 mm, got {delta:.4e}"
        );
    }

    #[test]
    fn test_skin_depth_decreases_with_conductivity() {
        let delta1 = skin_depth(2.0 * PI * 1e6, 1e6, 1.0);
        let delta2 = skin_depth(2.0 * PI * 1e6, 1e7, 1.0);
        assert!(delta2 < delta1);
    }

    // ---- impedance_surface --------------------------------------------------

    #[test]
    fn test_surface_impedance_equal_components() {
        let (rs, xs) = impedance_surface(2.0 * PI * 1e6, 5.8e7, 1.0);
        assert!(
            (rs - xs).abs() < 1e-15,
            "R_s must equal X_s for good conductors"
        );
    }

    #[test]
    fn test_surface_impedance_positive() {
        let (rs, xs) = impedance_surface(2.0 * PI * 1e6, 5.8e7, 1.0);
        assert!(rs > 0.0 && xs > 0.0);
    }

    // ---- EddyCurrent --------------------------------------------------------

    #[test]
    fn test_eddy_current_power_loss_positive() {
        let ec = EddyCurrent::new(1e3, 1e6, MU_0, 1e-3);
        assert!(ec.power_loss_per_volume() >= 0.0);
    }

    #[test]
    fn test_eddy_current_induced_emf_positive() {
        let ec = EddyCurrent::new(50.0, 5.8e7, MU_0, 0.01);
        assert!(ec.induced_emf(0.01) > 0.0);
    }

    #[test]
    fn test_eddy_current_emf_scales_with_area() {
        let ec = EddyCurrent::new(50.0, 5.8e7, MU_0, 0.01);
        let emf1 = ec.induced_emf(0.01);
        let emf2 = ec.induced_emf(0.02);
        assert!((emf2 - 2.0 * emf1).abs() < 1e-15);
    }

    // ---- helmholtz_1d_solution ----------------------------------------------

    #[test]
    fn test_helmholtz_surface_value() {
        let (re, im) = helmholtz_1d_solution(0.0, 0.01, 1.0);
        assert!((re - 1.0).abs() < 1e-10, "At x=0 Re = A0");
        assert!(im.abs() < 1e-10, "At x=0 Im = 0");
    }

    #[test]
    fn test_helmholtz_decays_with_depth() {
        let a0 = 1.0;
        let delta = 0.01;
        let (re1, im1) = helmholtz_1d_solution(0.001, delta, a0);
        let (re2, im2) = helmholtz_1d_solution(0.05, delta, a0);
        let amp1 = (re1 * re1 + im1 * im1).sqrt();
        let amp2 = (re2 * re2 + im2 * im2).sqrt();
        assert!(amp2 < amp1, "Amplitude must decay with depth");
    }

    #[test]
    fn test_helmholtz_zero_delta_returns_zero() {
        let (re, im) = helmholtz_1d_solution(0.01, 0.0, 1.0);
        assert_eq!(re, 0.0);
        assert_eq!(im, 0.0);
    }

    // ---- CoaxialLine --------------------------------------------------------

    #[test]
    fn test_coaxial_impedance_50_ohm() {
        // Standard 50Ω coax: ln(b/a) = 50*2π / (η ≈ 377/√ε_r)
        // For ε_r = 2.25 (polyethylene), b/a ≈ 3.48
        let inner = 0.445e-3;
        let outer = inner * 3.48;
        let coax = CoaxialLine::new(inner, outer, 5.8e7, 2.25, 1.0);
        let z = coax.impedance();
        assert!((z - 50.0).abs() < 2.0, "Expected ~50Ω, got {z:.2}");
    }

    #[test]
    fn test_coaxial_capacitance_positive() {
        let coax = CoaxialLine::new(1e-3, 5e-3, 5.8e7, 2.0, 1.0);
        assert!(coax.capacitance_per_length() > 0.0);
    }

    #[test]
    fn test_coaxial_inductance_positive() {
        let coax = CoaxialLine::new(1e-3, 5e-3, 5.8e7, 2.0, 1.0);
        assert!(coax.inductance_per_length() > 0.0);
    }

    #[test]
    fn test_coaxial_velocity_factor_less_than_one() {
        let coax = CoaxialLine::new(1e-3, 5e-3, 5.8e7, 2.25, 1.0);
        let vf = coax.velocity_factor();
        assert!(
            vf > 0.0 && vf < 1.0,
            "Velocity factor must be in (0, 1) for ε_r > 1"
        );
    }

    #[test]
    fn test_coaxial_inner_zero_returns_zero() {
        let coax = CoaxialLine::new(0.0, 5e-3, 5.8e7, 2.0, 1.0);
        assert_eq!(coax.capacitance_per_length(), 0.0);
        assert_eq!(coax.inductance_per_length(), 0.0);
    }

    // ---- SolenoidInductor ---------------------------------------------------

    #[test]
    fn test_solenoid_inductance_positive() {
        let sol = SolenoidInductor::new(100, 0.1, 0.02);
        assert!(sol.inductance() > 0.0);
    }

    #[test]
    fn test_solenoid_inductance_scales_n_squared() {
        let sol1 = SolenoidInductor::new(100, 0.1, 0.02);
        let sol2 = SolenoidInductor::new(200, 0.1, 0.02);
        let ratio = sol2.inductance() / sol1.inductance();
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "L ∝ N², ratio should be 4, got {ratio}"
        );
    }

    #[test]
    fn test_solenoid_field_positive() {
        let sol = SolenoidInductor::new(100, 0.1, 0.02);
        assert!(sol.field_at_center(1.0) > 0.0);
    }

    #[test]
    fn test_solenoid_field_scales_with_current() {
        let sol = SolenoidInductor::new(100, 0.1, 0.02);
        let b1 = sol.field_at_center(1.0);
        let b2 = sol.field_at_center(2.0);
        assert!((b2 - 2.0 * b1).abs() < 1e-20);
    }

    #[test]
    fn test_solenoid_zero_length_returns_zero() {
        let sol = SolenoidInductor::new(100, 0.0, 0.02);
        assert_eq!(sol.inductance(), 0.0);
        assert_eq!(sol.field_at_center(1.0), 0.0);
    }

    // ---- mutual_inductance_coaxial ------------------------------------------

    #[test]
    fn test_mutual_inductance_positive() {
        let m = mutual_inductance_coaxial(0.05, 0.05, 0.1);
        assert!(m > 0.0);
    }

    #[test]
    fn test_mutual_inductance_decreases_with_distance() {
        let m1 = mutual_inductance_coaxial(0.05, 0.05, 0.1);
        let m2 = mutual_inductance_coaxial(0.05, 0.05, 0.5);
        assert!(m2 < m1);
    }

    // ---- maxwell_stress_tensor ----------------------------------------------

    #[test]
    fn test_maxwell_stress_tensor_zero_fields() {
        let t = maxwell_stress_tensor([0.0; 3], [0.0; 3]);
        for row in &t {
            for &val in row.iter() {
                assert_eq!(val, 0.0);
            }
        }
    }

    #[test]
    fn test_maxwell_stress_tensor_diagonal_positive_e_field() {
        // E along x, B = 0: T_xx = ε₀·Ex²/2
        let e = [1000.0, 0.0, 0.0];
        let t = maxwell_stress_tensor(e, [0.0; 3]);
        assert!(t[0][0] > 0.0);
    }

    #[test]
    fn test_maxwell_stress_tensor_symmetric() {
        let e = [1e3, 2e3, 3e3];
        let b = [1e-3, 2e-3, 3e-3];
        let t = maxwell_stress_tensor(e, b);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (t[i][j] - t[j][i]).abs() < 1e-20,
                    "Tensor must be symmetric"
                );
            }
        }
    }

    // ---- poynting_vector ----------------------------------------------------

    #[test]
    fn test_poynting_vector_direction() {
        // E along x, B along y → S along z
        let e = [1.0, 0.0, 0.0];
        let b = [0.0, 1e-6, 0.0];
        let s = poynting_vector(e, b);
        assert!(s[2] > 0.0, "S_z must be positive for E_x × B_y");
        assert!(s[0].abs() < 1e-30);
        assert!(s[1].abs() < 1e-30);
    }

    #[test]
    fn test_poynting_vector_zero_fields() {
        let s = poynting_vector([0.0; 3], [0.0; 3]);
        assert_eq!(s, [0.0; 3]);
    }

    // ---- em_energy_density --------------------------------------------------

    #[test]
    fn test_energy_density_non_negative() {
        let u = em_energy_density([1e3, 0.0, 0.0], [1e-3, 0.0, 0.0]);
        assert!(u >= 0.0);
    }

    #[test]
    fn test_energy_density_zero_fields() {
        let u = em_energy_density([0.0; 3], [0.0; 3]);
        assert_eq!(u, 0.0);
    }

    #[test]
    fn test_energy_density_increases_with_field() {
        let u1 = em_energy_density([1e3, 0.0, 0.0], [0.0; 3]);
        let u2 = em_energy_density([2e3, 0.0, 0.0], [0.0; 3]);
        assert!(u2 > u1);
    }

    // ---- waveguide_te10_cutoff ----------------------------------------------

    #[test]
    fn test_waveguide_te10_cutoff_standard_x_band() {
        // X-band WR-90: a = 22.86 mm → fc ≈ 6.557 GHz
        let fc = waveguide_te10_cutoff(22.86e-3);
        assert!(
            (fc - 6.557e9).abs() < 0.01e9,
            "Expected ~6.56 GHz, got {fc:.3e}"
        );
    }

    #[test]
    fn test_waveguide_te10_cutoff_positive() {
        let fc = waveguide_te10_cutoff(0.05);
        assert!(fc > 0.0);
    }

    #[test]
    fn test_waveguide_te10_cutoff_decreases_with_width() {
        let fc1 = waveguide_te10_cutoff(0.02);
        let fc2 = waveguide_te10_cutoff(0.04);
        assert!(fc2 < fc1);
    }

    #[test]
    fn test_waveguide_te10_zero_width_infinite() {
        let fc = waveguide_te10_cutoff(0.0);
        assert!(fc.is_infinite());
    }

    // ---- ElectromagneticFEM struct ------------------------------------------

    #[test]
    fn test_electromagnetic_fem_skin_depth_matches_free_fn() {
        let fem = ElectromagneticFEM::new(100, 5.8e7, 1.0, 1.0, 2.0 * PI * 1e6);
        let delta_struct = fem.skin_depth();
        let delta_free = skin_depth(2.0 * PI * 1e6, 5.8e7, 1.0);
        assert!((delta_struct - delta_free).abs() < 1e-20);
    }

    #[test]
    fn test_electromagnetic_fem_surface_impedance_matches_free_fn() {
        let fem = ElectromagneticFEM::new(100, 5.8e7, 1.0, 1.0, 2.0 * PI * 1e6);
        let (rs_struct, xs_struct) = fem.surface_impedance();
        let (rs_free, xs_free) = impedance_surface(2.0 * PI * 1e6, 5.8e7, 1.0);
        assert!((rs_struct - rs_free).abs() < 1e-20);
        assert!((xs_struct - xs_free).abs() < 1e-20);
    }
}
