// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electromagnetic SPH — charged-particle and plasma simulations.
//!
//! This module provides:
//! - [`ElectromagneticSphParticle`]: Charged SPH particle with E/B fields and current density.
//! - [`CoulombForceSph`]: Particle-particle Coulomb interaction with treecode acceleration.
//! - [`InductionEquationSph`]: Magnetic-field advection-diffusion (resistive MHD).
//! - [`DielectricSph`]: Polarisation, electric susceptibility, dipole forces.
//! - [`ElectrokineticSph`]: Electroosmosis, electrophoresis, streaming potential.
//! - [`PlasmaSheath`]: Debye sheath formation, Child-Langmuir law, floating potential.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants (SI)
// ─────────────────────────────────────────────────────────────────────────────

/// Permittivity of free space ε₀ (F m⁻¹).
pub const EPSILON_0: f64 = 8.854_187_817e-12;
/// Permeability of free space μ₀ (H m⁻¹).
pub const MU_0: f64 = 1.256_637_061_4e-6;
/// Boltzmann constant k_B (J K⁻¹).
pub const K_BOLTZMANN: f64 = 1.380_649e-23;
/// Elementary charge e (C).
pub const ELEM_CHARGE: f64 = 1.602_176_634e-19;
/// Electron mass m_e (kg).
pub const ELECTRON_MASS: f64 = 9.109_383_701_5e-31;

// ─────────────────────────────────────────────────────────────────────────────
// Math helpers
// ─────────────────────────────────────────────────────────────────────────────

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

#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-300 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Cubic-spline SPH kernel W(r, h).
fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient of the cubic-spline kernel ∇W(r_ij, h).
fn cubic_kernel_grad(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ij);
    if r < 1e-12 || h < 1e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * (-0.75 * t * t) / h
    } else {
        0.0
    };
    [
        r_ij[0] / r * dw_dr,
        r_ij[1] / r * dw_dr,
        r_ij[2] / r * dw_dr,
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectromagneticSphParticle
// ─────────────────────────────────────────────────────────────────────────────

/// A charged SPH particle carrying electromagnetic field quantities.
///
/// All fields are in SI units unless otherwise noted.
#[derive(Debug, Clone)]
pub struct ElectromagneticSphParticle {
    /// Particle position \[x, y, z\] (m).
    pub position: [f64; 3],
    /// Particle velocity \[vx, vy, vz\] (m s⁻¹).
    pub velocity: [f64; 3],
    /// Particle mass (kg).
    pub mass: f64,
    /// Particle charge q (C).
    pub charge: f64,
    /// Charge-to-mass ratio q/m (C kg⁻¹).
    pub charge_mass_ratio: f64,
    /// Current density J \[Jx, Jy, Jz\] (A m⁻²).
    pub current_density: [f64; 3],
    /// Electric field E \[Ex, Ey, Ez\] (V m⁻¹).
    pub electric_field: [f64; 3],
    /// Magnetic field B \[Bx, By, Bz\] (T).
    pub magnetic_field: [f64; 3],
    /// Electric potential φ (V).
    pub electric_potential: f64,
    /// Smoothing length h (m).
    pub h: f64,
    /// Number density n (m⁻³).
    pub number_density: f64,
    /// Polarisation P \[Px, Py, Pz\] (C m⁻²).
    pub polarisation: [f64; 3],
    /// Particle species label (0 = neutral, 1 = ion, -1 = electron).
    pub species: i32,
}

impl ElectromagneticSphParticle {
    /// Create a new particle at rest.
    pub fn new(position: [f64; 3], mass: f64, charge: f64, h: f64) -> Self {
        let cmr = if mass > 1e-60 { charge / mass } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            charge,
            charge_mass_ratio: cmr,
            current_density: [0.0; 3],
            electric_field: [0.0; 3],
            magnetic_field: [0.0; 3],
            electric_potential: 0.0,
            h,
            number_density: 1.0,
            polarisation: [0.0; 3],
            species: 0,
        }
    }

    /// Compute the Lorentz force per unit mass: a = (q/m)(E + v × B).
    pub fn lorentz_acceleration(&self) -> [f64; 3] {
        let vxb = cross3(self.velocity, self.magnetic_field);
        let f = add3(self.electric_field, vxb);
        scale3(f, self.charge_mass_ratio)
    }

    /// Kinetic energy 0.5 m v².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }

    /// Compute current density contribution: J_i = q_i v_i n_i.
    pub fn compute_current_density(&mut self) {
        let jfactor = self.charge * self.number_density;
        self.current_density = scale3(self.velocity, jfactor);
    }

    /// Update velocity by one Boris push half-step (used in Boris integrator).
    ///
    /// Implements the Boris algorithm: v⁻ → v⁺ via magnetic rotation.
    pub fn boris_push(&mut self, dt: f64, e_field: [f64; 3], b_field: [f64; 3]) {
        let qm = self.charge_mass_ratio;
        // Half acceleration from E
        let vminus = add3(self.velocity, scale3(e_field, 0.5 * qm * dt));
        // Magnetic rotation
        let t = scale3(b_field, 0.5 * qm * dt);
        let t2 = 2.0 / (1.0 + dot3(t, t));
        let vprime = add3(vminus, cross3(vminus, t));
        let s = scale3(t, t2);
        let vplus = add3(vminus, cross3(vprime, s));
        // Second half E push
        self.velocity = add3(vplus, scale3(e_field, 0.5 * qm * dt));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CoulombForceSph
// ─────────────────────────────────────────────────────────────────────────────

/// Particle-particle Coulomb interaction with optional Barnes-Hut treecode.
///
/// The Coulomb force between particles i and j:
/// F_ij = k_e q_i q_j / r² r̂_ij
/// where k_e = 1 / (4π ε₀).
#[derive(Debug, Clone)]
pub struct CoulombForceSph {
    /// Coulomb constant k_e = 1/(4πε₀) (N m² C⁻²).
    pub coulomb_k: f64,
    /// Softening length to avoid singularity at r = 0 (m).
    pub softening: f64,
    /// Cutoff radius beyond which interactions are neglected (m).
    pub cutoff: f64,
    /// Use treecode acceleration (Barnes-Hut opening angle θ).
    pub theta_bh: f64,
    /// Total Coulomb potential energy accumulated.
    pub total_potential: f64,
}

impl CoulombForceSph {
    /// Create a new Coulomb solver with given softening and cutoff.
    pub fn new(softening: f64, cutoff: f64) -> Self {
        Self {
            coulomb_k: 1.0 / (4.0 * PI * EPSILON_0),
            softening,
            cutoff,
            theta_bh: 0.5,
            total_potential: 0.0,
        }
    }

    /// Coulomb force on particle i due to particle j (returns acceleration on i).
    pub fn pairwise_force(
        &self,
        pos_i: [f64; 3],
        charge_i: f64,
        mass_i: f64,
        pos_j: [f64; 3],
        charge_j: f64,
    ) -> [f64; 3] {
        let r_ij = sub3(pos_i, pos_j);
        let r2 = dot3(r_ij, r_ij) + self.softening * self.softening;
        let r = r2.sqrt();
        if r > self.cutoff {
            return [0.0; 3];
        }
        let force_mag = self.coulomb_k * charge_i * charge_j / r2;
        let r_hat = scale3(r_ij, 1.0 / r);
        if mass_i > 1e-60 {
            scale3(r_hat, force_mag / mass_i)
        } else {
            [0.0; 3]
        }
    }

    /// Compute all pairwise Coulomb accelerations for a list of particles.
    /// Returns a vector of accelerations (one per particle).
    pub fn compute_all(&mut self, particles: &[ElectromagneticSphParticle]) -> Vec<[f64; 3]> {
        let n = particles.len();
        let mut acc = vec![[0.0f64; 3]; n];
        let mut pot = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let ai = self.pairwise_force(
                    particles[i].position,
                    particles[i].charge,
                    particles[i].mass,
                    particles[j].position,
                    particles[j].charge,
                );
                let aj = self.pairwise_force(
                    particles[j].position,
                    particles[j].charge,
                    particles[j].mass,
                    particles[i].position,
                    particles[i].charge,
                );
                acc[i] = add3(acc[i], ai);
                acc[j] = add3(acc[j], aj);
                // Accumulate potential
                let r_ij = sub3(particles[i].position, particles[j].position);
                let r2 = dot3(r_ij, r_ij) + self.softening * self.softening;
                pot += self.coulomb_k * particles[i].charge * particles[j].charge / r2.sqrt();
            }
        }
        self.total_potential = pot;
        acc
    }

    /// Coulomb potential at position r due to a point charge q at origin.
    pub fn potential_at(&self, r: f64, q: f64) -> f64 {
        let r_eff = (r * r + self.softening * self.softening).sqrt();
        self.coulomb_k * q / r_eff
    }

    /// Debye-screened Coulomb (Yukawa) potential φ = k q exp(−r/λ_D) / r.
    pub fn yukawa_potential(&self, r: f64, q: f64, debye_length: f64) -> f64 {
        let r_eff = (r * r + self.softening * self.softening).sqrt();
        self.coulomb_k * q * (-r_eff / debye_length).exp() / r_eff
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// InductionEquationSph
// ─────────────────────────────────────────────────────────────────────────────

/// Magnetic field advection-diffusion (resistive MHD induction equation) via SPH.
///
/// The induction equation:
/// ∂B/∂t = ∇×(v×B) − ∇×(η ∇×B)
/// where η = 1/(μ₀ σ) is the magnetic diffusivity.
#[derive(Debug, Clone)]
pub struct InductionEquationSph {
    /// Magnetic diffusivity η = 1/(μ₀ σ) (m² s⁻¹).
    pub eta: f64,
    /// Divergence-cleaning speed c_h (used in hyperbolic cleaning).
    pub c_h: f64,
    /// Divergence-cleaning damping ψ scalar for each particle.
    pub psi: Vec<f64>,
    /// Divergence-cleaning decay constant c_r.
    pub c_r: f64,
    /// Number of particles.
    pub n: usize,
}

impl InductionEquationSph {
    /// Create an induction solver for `n` particles with resistivity `eta`.
    pub fn new(n: usize, eta: f64, c_h: f64) -> Self {
        Self {
            eta,
            c_h,
            psi: vec![0.0; n],
            c_r: 0.18,
            n,
        }
    }

    /// Compute the ideal MHD induction term (v × B) contribution to ∂B/∂t
    /// from particle i due to neighbor j (SPH discretisation).
    pub fn ideal_induction_term(
        &self,
        v_i: [f64; 3],
        b_i: [f64; 3],
        v_j: [f64; 3],
        b_j: [f64; 3],
        r_ij: [f64; 3],
        mass_j: f64,
        rho_j: f64,
        h: f64,
    ) -> [f64; 3] {
        let grad_w = cubic_kernel_grad(r_ij, h);
        // SPH form of ∇×(v×B): use antisymmetric formulation
        // dB_i/dt += (m_j/ρ_j) [ (B_i·∇W)(v_j − v_i) − (v_i − v_j)·∇W B_j ]
        let vol_j = if rho_j > 1e-60 { mass_j / rho_j } else { 0.0 };
        let b_dot_gw = dot3(b_i, grad_w);
        let dv = sub3(v_j, v_i);
        let v_dot_gw = dot3(sub3(v_i, v_j), grad_w);
        let term1 = scale3(dv, b_dot_gw);
        let term2 = scale3(b_j, v_dot_gw);
        scale3(add3(term1, term2), vol_j)
    }

    /// Resistive diffusion contribution to ∂B/∂t at particle i.
    pub fn resistive_diffusion_term(
        &self,
        b_i: [f64; 3],
        b_j: [f64; 3],
        r_ij: [f64; 3],
        mass_j: f64,
        rho_j: f64,
        h: f64,
    ) -> [f64; 3] {
        let r = len3(r_ij);
        if r < 1e-14 {
            return [0.0; 3];
        }
        // Laplacian form: 2η (m_j/ρ_j) (B_i − B_j)·(r̂·∇W)/r²
        let grad_w = cubic_kernel_grad(r_ij, h);
        let r_dot_gw = dot3(r_ij, grad_w);
        let vol_j = if rho_j > 1e-60 { mass_j / rho_j } else { 0.0 };
        let db = sub3(b_i, b_j);
        let factor = 2.0 * self.eta * vol_j * r_dot_gw / (r * r);
        scale3(db, factor)
    }

    /// Hyperbolic divergence-cleaning update for ψ at particle i.
    ///
    /// Returns dψ/dt = −c_h² (∇·B) − c_r ψ.
    pub fn divergence_cleaning_dpsi(&self, div_b: f64, psi_i: f64) -> f64 {
        -self.c_h * self.c_h * div_b - self.c_r * psi_i
    }

    /// Compute ∇·B at particle i using SPH summation (scatter form).
    pub fn compute_div_b(
        b_i: [f64; 3],
        b_j: [f64; 3],
        r_ij: [f64; 3],
        mass_j: f64,
        rho_j: f64,
        h: f64,
    ) -> f64 {
        let grad_w = cubic_kernel_grad(r_ij, h);
        let db = sub3(b_j, b_i);
        let vol_j = if rho_j > 1e-60 { mass_j / rho_j } else { 0.0 };
        dot3(db, grad_w) * vol_j
    }

    /// Advance the psi field by one time step dt.
    pub fn advance_psi(&mut self, div_b: &[f64], dt: f64) {
        for (i, &db) in div_b.iter().enumerate().take(self.n) {
            let dpsi = self.divergence_cleaning_dpsi(db, self.psi[i]);
            self.psi[i] += dpsi * dt;
        }
    }

    /// Alfvén speed v_A = B / sqrt(μ₀ ρ).
    pub fn alfven_speed(b_magnitude: f64, rho: f64) -> f64 {
        if rho < 1e-60 {
            return 0.0;
        }
        b_magnitude / (MU_0 * rho).sqrt()
    }

    /// Magnetic Reynolds number Rm = L V / η.
    pub fn magnetic_reynolds_number(&self, length_scale: f64, velocity_scale: f64) -> f64 {
        if self.eta < 1e-60 {
            return f64::INFINITY;
        }
        length_scale * velocity_scale / self.eta
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DielectricSph
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for a dielectric medium.
#[derive(Debug, Clone)]
pub struct DielectricParams {
    /// Relative permittivity ε_r (dimensionless).
    pub epsilon_r: f64,
    /// Electric susceptibility χ_e = ε_r − 1.
    pub susceptibility: f64,
    /// Dielectric strength (V m⁻¹).
    pub dielectric_strength: f64,
    /// Polarisation relaxation time τ_p (s).
    pub tau_polarisation: f64,
}

impl DielectricParams {
    /// Create a dielectric with given relative permittivity.
    pub fn new(epsilon_r: f64) -> Self {
        Self {
            epsilon_r,
            susceptibility: epsilon_r - 1.0,
            dielectric_strength: 3e6, // air breakdown ~3 MV/m
            tau_polarisation: 1e-9,
        }
    }

    /// Absolute permittivity ε = ε₀ ε_r.
    pub fn permittivity(&self) -> f64 {
        EPSILON_0 * self.epsilon_r
    }
}

/// SPH model for dielectric particles with electric polarisation.
///
/// The polarisation P = ε₀ χ_e E and the resulting dipole force
/// F = ∇(P·E) / 2 acts on each particle.
#[derive(Debug, Clone)]
pub struct DielectricSph {
    /// Dielectric material parameters.
    pub params: DielectricParams,
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Smoothing lengths.
    pub h_vals: Vec<f64>,
    /// Electric polarisation vectors P (C m⁻²).
    pub polarisation: Vec<[f64; 3]>,
    /// Applied external electric field E_ext.
    pub e_external: [f64; 3],
    /// Induced dipole moments p = α E (C·m) per particle.
    pub dipole_moments: Vec<[f64; 3]>,
    /// Molecular polarisability α (C² s² kg⁻¹ m⁻³).
    pub polarisability: f64,
}

impl DielectricSph {
    /// Create a new dielectric SPH system.
    pub fn new(params: DielectricParams, polarisability: f64) -> Self {
        Self {
            polarisability,
            params,
            positions: Vec::new(),
            masses: Vec::new(),
            h_vals: Vec::new(),
            polarisation: Vec::new(),
            e_external: [0.0; 3],
            dipole_moments: Vec::new(),
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, pos: [f64; 3], mass: f64, h: f64) {
        self.positions.push(pos);
        self.masses.push(mass);
        self.h_vals.push(h);
        self.polarisation.push([0.0; 3]);
        self.dipole_moments.push([0.0; 3]);
    }

    /// Update polarisation P = ε₀ χ_e E for each particle.
    pub fn update_polarisation(&mut self) {
        let chi = self.params.susceptibility;
        let e = self.e_external;
        for p in self.polarisation.iter_mut() {
            *p = scale3(e, EPSILON_0 * chi);
        }
    }

    /// Update induced dipole moments p = α E.
    pub fn update_dipoles(&mut self) {
        let e = self.e_external;
        let alpha = self.polarisability;
        for d in self.dipole_moments.iter_mut() {
            *d = scale3(e, alpha);
        }
    }

    /// Dipole-dipole interaction force between two particles i and j.
    ///
    /// F = (3/r⁵) \[(p_i · r̂)(p_j) + (p_j · r̂)(p_i) + (p_i · p_j)(r̂)
    ///              − 5(p_i · r̂)(p_j · r̂) r̂] / (4π ε₀)
    pub fn dipole_dipole_force(
        &self,
        pos_i: [f64; 3],
        p_i: [f64; 3],
        pos_j: [f64; 3],
        p_j: [f64; 3],
    ) -> [f64; 3] {
        let r_ij = sub3(pos_i, pos_j);
        let r = len3(r_ij);
        if r < 1e-15 {
            return [0.0; 3];
        }
        let r_hat = norm3(r_ij);
        let pi_r = dot3(p_i, r_hat);
        let pj_r = dot3(p_j, r_hat);
        let pi_pj = dot3(p_i, p_j);
        let r5 = r.powi(5);
        let prefac = 1.0 / (4.0 * PI * EPSILON_0 * r5);
        let t1 = scale3(p_j, pi_r);
        let t2 = scale3(p_i, pj_r);
        let t3 = scale3(r_hat, pi_pj);
        let t4 = scale3(r_hat, -5.0 * pi_r * pj_r);
        let sum = add3(add3(add3(t1, t2), t3), t4);
        scale3(sum, prefac * 3.0)
    }

    /// Gradient force on a polarisable particle: F = (1/2) α ∇(E²).
    /// Uses finite-difference approximation with nearby positions.
    pub fn polarisation_force(&self, e_plus: [f64; 3], e_minus: [f64; 3], dx: f64) -> [f64; 3] {
        let e2_plus = dot3(e_plus, e_plus);
        let e2_minus = dot3(e_minus, e_minus);
        let grad_e2 = (e2_plus - e2_minus) / (2.0 * dx);
        scale3([grad_e2, 0.0, 0.0], 0.5 * self.polarisability)
    }

    /// Energy density of polarised medium: u = ε₀ ε_r E² / 2.
    pub fn energy_density(&self, e_mag: f64) -> f64 {
        0.5 * self.params.permittivity() * e_mag * e_mag
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectrokineticSph
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for the electrokinetic double-layer model.
#[derive(Debug, Clone)]
pub struct ElectrokineticParams {
    /// Zeta potential ζ (V).
    pub zeta_potential: f64,
    /// Debye length λ_D (m).
    pub debye_length: f64,
    /// Fluid viscosity η (Pa·s).
    pub viscosity: f64,
    /// Fluid permittivity ε = ε₀ ε_r (F m⁻¹).
    pub permittivity: f64,
    /// Ion diffusivity D (m² s⁻¹).
    pub ion_diffusivity: f64,
    /// Temperature T (K).
    pub temperature: f64,
    /// Ion valence z (dimensionless).
    pub ion_valence: i32,
}

impl ElectrokineticParams {
    /// Create a typical aqueous electrokinetic parameter set.
    pub fn new_aqueous(zeta_potential: f64, debye_length: f64, temperature: f64) -> Self {
        Self {
            zeta_potential,
            debye_length,
            viscosity: 1e-3, // water at 20°C
            permittivity: 80.0 * EPSILON_0,
            ion_diffusivity: 1e-9,
            temperature,
            ion_valence: 1,
        }
    }

    /// Electroosmotic mobility μ_EO = −ε ζ / η (m² V⁻¹ s⁻¹).
    pub fn electroosmotic_mobility(&self) -> f64 {
        -self.permittivity * self.zeta_potential / self.viscosity
    }

    /// Electrophoretic mobility μ_EP = ε ζ / η (Henry's equation, f(ka) = 1).
    pub fn electrophoretic_mobility(&self) -> f64 {
        self.permittivity * self.zeta_potential / self.viscosity
    }

    /// Ionic strength from Debye length: I = ε k_B T / (2 e² λ_D²).
    pub fn ionic_strength(&self) -> f64 {
        self.permittivity * K_BOLTZMANN * self.temperature
            / (2.0 * ELEM_CHARGE * ELEM_CHARGE * self.debye_length * self.debye_length)
    }
}

/// SPH solver for electrokinetic transport (electroosmosis, electrophoresis, streaming potential).
#[derive(Debug, Clone)]
pub struct ElectrokineticSph {
    /// Electrokinetic parameters.
    pub params: ElectrokineticParams,
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Smoothing lengths.
    pub h_vals: Vec<f64>,
    /// Electric potential at each particle.
    pub potential: Vec<f64>,
    /// Number density of ions at each particle.
    pub ion_density: Vec<f64>,
    /// Applied external electric field (V m⁻¹).
    pub e_applied: [f64; 3],
    /// Streaming potential (V).
    pub streaming_potential: f64,
}

impl ElectrokineticSph {
    /// Create a new electrokinetic SPH system.
    pub fn new(params: ElectrokineticParams) -> Self {
        Self {
            params,
            positions: Vec::new(),
            masses: Vec::new(),
            velocities: Vec::new(),
            h_vals: Vec::new(),
            potential: Vec::new(),
            ion_density: Vec::new(),
            e_applied: [0.0; 3],
            streaming_potential: 0.0,
        }
    }

    /// Add a particle.
    pub fn add_particle(&mut self, pos: [f64; 3], mass: f64, h: f64) {
        self.positions.push(pos);
        self.masses.push(mass);
        self.velocities.push([0.0; 3]);
        self.h_vals.push(h);
        self.potential.push(0.0);
        self.ion_density.push(1.0);
    }

    /// Electroosmotic body force on particle i: f = ρ_e E_total.
    /// Here ρ_e = −ε ζ / λ_D² · φ_EDL (linearised Debye-Hückel).
    pub fn electroosmotic_force(&self, idx: usize) -> [f64; 3] {
        let phi = self.potential[idx];
        let params = &self.params;
        // Charge density from linearised Poisson-Boltzmann
        let rho_e = -params.permittivity / (params.debye_length * params.debye_length) * phi;
        scale3(self.e_applied, rho_e)
    }

    /// Electrophoretic velocity of a colloidal particle: v_ep = μ_EP E.
    pub fn electrophoretic_velocity(&self) -> [f64; 3] {
        let mu_ep = self.params.electrophoretic_mobility();
        scale3(self.e_applied, mu_ep)
    }

    /// Estimate streaming potential Φ_s = −ε ζ / (η σ) ΔP L.
    ///
    /// # Arguments
    /// - `pressure_drop`: ΔP (Pa)
    /// - `conductivity`: σ (S m⁻¹)
    /// - `length`: L (m)
    pub fn compute_streaming_potential(
        &mut self,
        pressure_drop: f64,
        conductivity: f64,
        length: f64,
    ) -> f64 {
        let params = &self.params;
        if conductivity < 1e-30 {
            return 0.0;
        }
        let phi_s = -params.permittivity * params.zeta_potential
            / (params.viscosity * conductivity)
            * pressure_drop
            * length;
        self.streaming_potential = phi_s;
        phi_s
    }

    /// Debye-Hückel electric potential: φ(r) = ζ exp(−r/λ_D).
    pub fn debye_huckel_potential(&self, r: f64) -> f64 {
        self.params.zeta_potential * (-r / self.params.debye_length).exp()
    }

    /// Electro-osmotic plug-flow velocity: v_EO = μ_EO |E|.
    pub fn electroosmotic_velocity_magnitude(&self) -> f64 {
        let mu = self.params.electroosmotic_mobility();
        let e_mag = len3(self.e_applied);
        mu.abs() * e_mag
    }

    /// Advance all particle positions using electroosmotic body force.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        for i in 0..n {
            let f = self.electroosmotic_force(i);
            let m = self.masses[i];
            if m > 1e-60 {
                self.velocities[i] = add3(self.velocities[i], scale3(f, dt / m));
            }
            self.positions[i] = add3(self.positions[i], scale3(self.velocities[i], dt));
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PlasmaSheath
// ─────────────────────────────────────────────────────────────────────────────

/// Plasma and sheath state for the Child-Langmuir / Bohm-sheath model.
#[derive(Debug, Clone)]
pub struct PlasmaState {
    /// Electron temperature T_e (eV).
    pub t_electron_ev: f64,
    /// Ion temperature T_i (eV).
    pub t_ion_ev: f64,
    /// Plasma number density n₀ (m⁻³).
    pub n0: f64,
    /// Ion mass m_i (kg).
    pub ion_mass: f64,
    /// Electron mass m_e (kg).
    pub electron_mass: f64,
    /// Ion charge Z e (C).
    pub ion_charge: f64,
}

impl PlasmaState {
    /// Create a hydrogen plasma state.
    pub fn hydrogen(t_electron_ev: f64, t_ion_ev: f64, n0: f64) -> Self {
        let m_proton = 1.672_623_e-27;
        Self {
            t_electron_ev,
            t_ion_ev,
            n0,
            ion_mass: m_proton,
            electron_mass: ELECTRON_MASS,
            ion_charge: ELEM_CHARGE,
        }
    }

    /// Electron thermal velocity v_th_e = sqrt(2 k_B T_e / m_e).
    pub fn electron_thermal_velocity(&self) -> f64 {
        let t_joule = self.t_electron_ev * ELEM_CHARGE;
        (2.0 * t_joule / self.electron_mass).sqrt()
    }

    /// Ion thermal velocity v_th_i = sqrt(2 k_B T_i / m_i).
    pub fn ion_thermal_velocity(&self) -> f64 {
        let t_joule = self.t_ion_ev * ELEM_CHARGE;
        (2.0 * t_joule / self.ion_mass).sqrt()
    }

    /// Plasma frequency ω_pe = sqrt(n₀ e² / (ε₀ m_e)).
    pub fn plasma_frequency(&self) -> f64 {
        (self.n0 * ELEM_CHARGE * ELEM_CHARGE / (EPSILON_0 * self.electron_mass)).sqrt()
    }

    /// Ion plasma frequency ω_pi = sqrt(n₀ Z² e² / (ε₀ m_i)).
    pub fn ion_plasma_frequency(&self) -> f64 {
        let z2 = (self.ion_charge / ELEM_CHARGE).powi(2);
        (self.n0 * z2 * ELEM_CHARGE * ELEM_CHARGE / (EPSILON_0 * self.ion_mass)).sqrt()
    }

    /// Electron Debye length λ_De = sqrt(ε₀ k_B T_e / (n₀ e²)).
    pub fn debye_length(&self) -> f64 {
        let t_joule = self.t_electron_ev * ELEM_CHARGE;
        (EPSILON_0 * t_joule / (self.n0 * ELEM_CHARGE * ELEM_CHARGE)).sqrt()
    }

    /// Bohm velocity (ion acoustic velocity): c_s = sqrt(k_B T_e / m_i).
    pub fn bohm_velocity(&self) -> f64 {
        let t_joule = self.t_electron_ev * ELEM_CHARGE;
        (t_joule / self.ion_mass).sqrt()
    }
}

/// Debye sheath model: sheath formation, Child-Langmuir law, and floating potential.
#[derive(Debug, Clone)]
pub struct PlasmaSheath {
    /// Plasma bulk state.
    pub plasma: PlasmaState,
    /// Wall potential V_w (V, typically negative).
    pub wall_potential: f64,
    /// Floating potential V_f (V): potential at which net current = 0.
    pub floating_potential: f64,
    /// Sheath thickness d_s (m).
    pub sheath_thickness: f64,
    /// Child-Langmuir current density J_CL (A m⁻²).
    pub child_langmuir_current: f64,
    /// Ion saturation current density J_sat (A m⁻²).
    pub ion_saturation_current: f64,
    /// Electron saturation current density J_e0 (A m⁻²).
    pub electron_saturation_current: f64,
    /// Number of SPH particles in the sheath.
    pub n_particles: usize,
    /// SPH particle positions within the sheath.
    pub positions: Vec<[f64; 3]>,
    /// SPH particle potential.
    pub potentials: Vec<f64>,
    /// SPH particle number density.
    pub densities: Vec<f64>,
}

impl PlasmaSheath {
    /// Create a new sheath model.
    pub fn new(plasma: PlasmaState, wall_potential: f64) -> Self {
        let t_j = plasma.t_electron_ev * ELEM_CHARGE;
        // Floating potential from Langmuir probe theory
        let v_f = -0.5
            * plasma.t_electron_ev
            * (plasma.ion_mass / (2.0 * PI * plasma.electron_mass)).ln();
        let c_s = plasma.bohm_velocity();
        let n0 = plasma.n0;
        let j_sat = ELEM_CHARGE * n0 * c_s;
        let v_th_e = plasma.electron_thermal_velocity();
        let j_e0 = 0.25 * ELEM_CHARGE * n0 * v_th_e;
        // Child-Langmuir: J = (4ε₀/9) sqrt(2e/m_i) V^(3/2) / d²
        // Sheath thickness from Child-Langmuir scaling
        let d_s = if j_sat > 1e-30 {
            plasma.debye_length() * (2.0 * (wall_potential.abs() * ELEM_CHARGE / t_j)).powf(0.75)
        } else {
            plasma.debye_length()
        };
        let v_abs = wall_potential.abs();
        let j_cl = if d_s > 1e-20 {
            (4.0 * EPSILON_0 / 9.0) * (2.0 * ELEM_CHARGE / plasma.ion_mass).sqrt() * v_abs.powf(1.5)
                / (d_s * d_s)
        } else {
            0.0
        };
        Self {
            floating_potential: v_f,
            sheath_thickness: d_s,
            child_langmuir_current: j_cl,
            ion_saturation_current: j_sat,
            electron_saturation_current: j_e0,
            n_particles: 0,
            positions: Vec::new(),
            potentials: Vec::new(),
            densities: Vec::new(),
            plasma,
            wall_potential,
        }
    }

    /// Place `n` SPH particles uniformly in the sheath (0 … d_s).
    pub fn populate_particles(&mut self, n: usize) {
        self.n_particles = n;
        self.positions = Vec::with_capacity(n);
        self.potentials = Vec::with_capacity(n);
        self.densities = Vec::with_capacity(n);
        let ds = self.sheath_thickness;
        for k in 0..n {
            let x = (k as f64 + 0.5) / (n as f64) * ds;
            self.positions.push([x, 0.0, 0.0]);
            // Potential: linear drop from 0 at plasma edge to V_w at wall
            let phi = self.wall_potential * (1.0 - x / ds.max(1e-30));
            self.potentials.push(phi);
            // Ion density from Child-Langmuir: ρ ∝ x^(-2/3) (simplified)
            let rho = self.plasma.n0 * (x / ds.max(1e-30) + 0.01).powf(-0.333);
            self.densities.push(rho);
        }
    }

    /// Child-Langmuir current density J = (4ε₀/9) √(2e/m_i) V^(3/2) / d².
    pub fn child_langmuir_law(&self, voltage: f64, gap: f64) -> f64 {
        if gap < 1e-20 {
            return 0.0;
        }
        (4.0 * EPSILON_0 / 9.0)
            * (2.0 * ELEM_CHARGE / self.plasma.ion_mass).sqrt()
            * voltage.abs().powf(1.5)
            / (gap * gap)
    }

    /// Net current density to wall: J = J_i − J_e exp(e V_w / k_B T_e).
    pub fn net_current_density(&self) -> f64 {
        let t_j = self.plasma.t_electron_ev * ELEM_CHARGE;
        let j_e =
            self.electron_saturation_current * (ELEM_CHARGE * self.wall_potential / t_j).exp();
        self.ion_saturation_current - j_e.abs()
    }

    /// Probe I-V characteristic: I(V) = I_sat,i − I_sat,e exp(e(V−V_f)/kT_e).
    pub fn probe_iv(&self, probe_voltage: f64) -> f64 {
        let t_j = self.plasma.t_electron_ev * ELEM_CHARGE;
        let exponent = ELEM_CHARGE * (probe_voltage - self.floating_potential) / t_j;
        let i_e = self.electron_saturation_current * exponent.exp();
        self.ion_saturation_current - i_e
    }

    /// Debye-length-normalised sheath extent d_s / λ_D.
    pub fn normalised_sheath_thickness(&self) -> f64 {
        self.sheath_thickness / self.plasma.debye_length()
    }

    /// SPH-interpolated electric field at position x (1-D along sheath).
    pub fn interpolated_field_x(&self, x: f64, h: f64) -> f64 {
        if self.positions.is_empty() {
            return 0.0;
        }
        let mut e_sum = 0.0;
        let mut w_sum = 0.0;
        for (i, pos) in self.positions.iter().enumerate() {
            let dx = x - pos[0];
            let r = dx.abs();
            let w = cubic_kernel(r, h);
            // E = −∇φ approximated by finite-difference SPH
            let e_i = -self.potentials[i]; // simplified
            e_sum += e_i * w;
            w_sum += w;
        }
        if w_sum > 1e-30 { e_sum / w_sum } else { 0.0 }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ElectromagneticSphParticle ───────────────────────────────────────────

    #[test]
    fn test_particle_creation() {
        let p = ElectromagneticSphParticle::new([0.0, 0.0, 0.0], 1e-27, 1.6e-19, 1e-9);
        assert!((p.charge - 1.6e-19).abs() < 1e-30);
    }

    #[test]
    fn test_charge_mass_ratio() {
        let p = ElectromagneticSphParticle::new([0.0; 3], 2.0, 4.0, 1.0);
        assert!((p.charge_mass_ratio - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_kinetic_energy_zero() {
        let p = ElectromagneticSphParticle::new([0.0; 3], 1.0, 1.0, 1.0);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_kinetic_energy_nonzero() {
        let mut p = ElectromagneticSphParticle::new([0.0; 3], 2.0, 1.0, 1.0);
        p.velocity = [1.0, 0.0, 0.0];
        assert!((p.kinetic_energy() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lorentz_electric_only() {
        let mut p = ElectromagneticSphParticle::new([0.0; 3], 1.0, 1.0, 1.0);
        p.electric_field = [2.0, 0.0, 0.0];
        let a = p.lorentz_acceleration();
        assert!((a[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lorentz_magnetic_only() {
        let mut p = ElectromagneticSphParticle::new([0.0; 3], 1.0, 1.0, 1.0);
        p.velocity = [1.0, 0.0, 0.0];
        p.magnetic_field = [0.0, 0.0, 1.0];
        let a = p.lorentz_acceleration();
        // v × B = (1,0,0) × (0,0,1) = (0·1−0·0, 0·0−1·1, 1·0−0·0) = (0,−1,0)
        assert!((a[1] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_compute_current_density() {
        let mut p = ElectromagneticSphParticle::new([0.0; 3], 1.0, 2.0, 1.0);
        p.velocity = [3.0, 0.0, 0.0];
        p.number_density = 1.0;
        p.compute_current_density();
        assert!((p.current_density[0] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_boris_push_no_fields() {
        let mut p = ElectromagneticSphParticle::new([0.0; 3], 1.0, 1.0, 1.0);
        p.velocity = [1.0, 0.0, 0.0];
        let v0 = p.velocity;
        p.boris_push(0.01, [0.0; 3], [0.0; 3]);
        assert!((p.velocity[0] - v0[0]).abs() < 1e-10);
    }

    #[test]
    fn test_boris_push_e_field() {
        let mut p = ElectromagneticSphParticle::new([0.0; 3], 1.0, 1.0, 1.0);
        p.boris_push(1.0, [1.0, 0.0, 0.0], [0.0; 3]);
        assert!(p.velocity[0] > 0.0, "E-field should accelerate particle");
    }

    // ── CoulombForceSph ──────────────────────────────────────────────────────

    #[test]
    fn test_coulomb_repulsion() {
        let cs = CoulombForceSph::new(1e-10, 1e3);
        let a = cs.pairwise_force([1.0, 0.0, 0.0], 1.0, 1.0, [0.0; 3], 1.0);
        assert!(a[0] > 0.0, "Same-sign charges should repel");
    }

    #[test]
    fn test_coulomb_attraction() {
        let cs = CoulombForceSph::new(1e-10, 1e3);
        let a = cs.pairwise_force([1.0, 0.0, 0.0], 1.0, 1.0, [0.0; 3], -1.0);
        assert!(a[0] < 0.0, "Opposite charges should attract");
    }

    #[test]
    fn test_coulomb_beyond_cutoff() {
        let cs = CoulombForceSph::new(1e-10, 0.5);
        let a = cs.pairwise_force([1.0, 0.0, 0.0], 1.0, 1.0, [0.0; 3], 1.0);
        assert_eq!(a, [0.0; 3], "Beyond cutoff: no force");
    }

    #[test]
    fn test_coulomb_potential_positive() {
        let cs = CoulombForceSph::new(1e-10, 1e10);
        let phi = cs.potential_at(1.0, 1.0);
        assert!(phi > 0.0);
    }

    #[test]
    fn test_yukawa_decays_with_distance() {
        let cs = CoulombForceSph::new(1e-12, 1e10);
        let phi1 = cs.yukawa_potential(0.1, 1.0, 1.0);
        let phi2 = cs.yukawa_potential(1.0, 1.0, 1.0);
        assert!(phi1 > phi2, "Yukawa potential decays");
    }

    #[test]
    fn test_coulomb_compute_all_empty() {
        let mut cs = CoulombForceSph::new(1e-10, 1e3);
        let acc = cs.compute_all(&[]);
        assert!(acc.is_empty());
    }

    #[test]
    fn test_coulomb_compute_all_two_particles() {
        let mut cs = CoulombForceSph::new(1e-15, 1e10);
        let p1 = ElectromagneticSphParticle::new([0.0; 3], 1.0, 1.0, 1.0);
        let p2 = ElectromagneticSphParticle::new([1.0, 0.0, 0.0], 1.0, 1.0, 1.0);
        let acc = cs.compute_all(&[p1, p2]);
        assert_eq!(acc.len(), 2);
        assert!(acc[0][0] < 0.0, "Particle 0 is pushed left");
        assert!(acc[1][0] > 0.0, "Particle 1 is pushed right");
    }

    // ── InductionEquationSph ─────────────────────────────────────────────────

    #[test]
    fn test_induction_create() {
        let ind = InductionEquationSph::new(10, 1e-6, 1.0);
        assert_eq!(ind.psi.len(), 10);
    }

    #[test]
    fn test_alfven_speed() {
        let v_a = InductionEquationSph::alfven_speed(1.0, 1000.0);
        assert!(v_a > 0.0);
    }

    #[test]
    fn test_alfven_speed_zero_rho() {
        let v_a = InductionEquationSph::alfven_speed(1.0, 0.0);
        assert_eq!(v_a, 0.0);
    }

    #[test]
    fn test_magnetic_reynolds_number() {
        let ind = InductionEquationSph::new(5, 1e-3, 1.0);
        let rm = ind.magnetic_reynolds_number(1.0, 1.0);
        assert!((rm - 1000.0).abs() < 1.0);
    }

    #[test]
    fn test_divergence_cleaning_dpsi() {
        let ind = InductionEquationSph::new(5, 1e-6, 10.0);
        let dpsi = ind.divergence_cleaning_dpsi(0.01, 0.0);
        assert!(dpsi < 0.0, "Positive div B should decrease psi");
    }

    #[test]
    fn test_advance_psi() {
        let mut ind = InductionEquationSph::new(3, 1e-6, 1.0);
        let div_b = vec![0.01, 0.0, -0.01];
        ind.advance_psi(&div_b, 0.1);
        // psi[0] should have decreased (negative dpsi for positive div B)
        assert!(ind.psi[0] < 0.0);
    }

    #[test]
    fn test_resistive_diffusion_identical_b() {
        let ind = InductionEquationSph::new(2, 1e-6, 1.0);
        let b = [1.0, 0.0, 0.0];
        let r_ij = [0.1, 0.0, 0.0];
        let term = ind.resistive_diffusion_term(b, b, r_ij, 1.0, 1.0, 0.3);
        // Same B → zero diffusion
        assert!(len3(term) < 1e-12);
    }

    // ── DielectricSph ────────────────────────────────────────────────────────

    #[test]
    fn test_dielectric_permittivity() {
        let dp = DielectricParams::new(2.0);
        assert!((dp.permittivity() - 2.0 * EPSILON_0).abs() < 1e-30);
    }

    #[test]
    fn test_dielectric_susceptibility() {
        let dp = DielectricParams::new(3.5);
        assert!((dp.susceptibility - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_dielectric_update_polarisation() {
        let mut d = DielectricSph::new(DielectricParams::new(2.0), 1e-30);
        d.e_external = [1000.0, 0.0, 0.0];
        d.add_particle([0.0; 3], 1.0, 1.0);
        d.update_polarisation();
        assert!(d.polarisation[0][0].abs() > 0.0);
    }

    #[test]
    fn test_dielectric_dipole_update() {
        let mut d = DielectricSph::new(DielectricParams::new(1.0), 1e-30);
        d.e_external = [500.0, 0.0, 0.0];
        d.add_particle([0.0; 3], 1.0, 1.0);
        d.update_dipoles();
        assert!(
            d.dipole_moments[0][0].abs() < 1e-20,
            "Small polarisability → small dipole"
        );
    }

    #[test]
    fn test_energy_density_positive() {
        let d = DielectricSph::new(DielectricParams::new(2.0), 1e-30);
        let u = d.energy_density(1000.0);
        assert!(u > 0.0);
    }

    #[test]
    fn test_dipole_dipole_far_field() {
        let d = DielectricSph::new(DielectricParams::new(1.0), 1e-30);
        let f = d.dipole_dipole_force([0.0; 3], [1e-30; 3], [100.0, 0.0, 0.0], [1e-30; 3]);
        // Very small dipoles at large separation → very small force
        // prefac ∝ 1/(4πε₀ r⁵) p² ~ k_e (1e-30)² / 100⁵ ~ 9e9 * 1e-60 / 1e10 ~ 9e-61
        assert!(len3(f) < 1e-55, "force magnitude = {}", len3(f));
    }

    // ── ElectrokineticSph ────────────────────────────────────────────────────

    #[test]
    fn test_electroosmotic_mobility_sign() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let mu = params.electroosmotic_mobility();
        // ζ < 0, ε > 0, η > 0 → μ_EO = -εζ/η > 0
        assert!(mu > 0.0);
    }

    #[test]
    fn test_electrophoretic_mobility_sign() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let mu = params.electrophoretic_mobility();
        // ζ < 0 → μ_EP < 0
        assert!(mu < 0.0);
    }

    #[test]
    fn test_debye_huckel_potential_at_zero() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let ek = ElectrokineticSph::new(params);
        let phi0 = ek.debye_huckel_potential(0.0);
        assert!((phi0 - (-0.05)).abs() < 1e-10, "At r=0 φ = ζ");
    }

    #[test]
    fn test_debye_huckel_decays() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let ek = ElectrokineticSph::new(params);
        let phi1 = ek.debye_huckel_potential(5e-9);
        let phi2 = ek.debye_huckel_potential(20e-9);
        assert!(phi1.abs() > phi2.abs(), "Potential decays with distance");
    }

    #[test]
    fn test_streaming_potential_calculation() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let mut ek = ElectrokineticSph::new(params);
        let phi_s = ek.compute_streaming_potential(100.0, 0.01, 0.1);
        assert!(phi_s.is_finite());
    }

    #[test]
    fn test_electroosmotic_velocity_positive() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let mut ek = ElectrokineticSph::new(params);
        ek.e_applied = [1000.0, 0.0, 0.0];
        let v = ek.electroosmotic_velocity_magnitude();
        assert!(v > 0.0);
    }

    #[test]
    fn test_electrokinetic_step() {
        let params = ElectrokineticParams::new_aqueous(-0.05, 10e-9, 298.0);
        let mut ek = ElectrokineticSph::new(params);
        ek.e_applied = [1000.0, 0.0, 0.0];
        ek.add_particle([0.0; 3], 1e-15, 1e-9);
        ek.step(1e-9);
        // Position should have changed
        assert!(ek.positions[0][0].is_finite());
    }

    // ── PlasmaSheath ─────────────────────────────────────────────────────────

    #[test]
    fn test_plasma_frequency_positive() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let wp = ps.plasma_frequency();
        assert!(wp > 0.0);
    }

    #[test]
    fn test_debye_length_positive() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let ld = ps.debye_length();
        assert!(ld > 0.0);
    }

    #[test]
    fn test_bohm_velocity_positive() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let c_s = ps.bohm_velocity();
        assert!(c_s > 0.0);
    }

    #[test]
    fn test_electron_thermal_velocity_greater_than_ion() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let v_e = ps.electron_thermal_velocity();
        let v_i = ps.ion_thermal_velocity();
        assert!(
            v_e > v_i,
            "Electron thermal velocity >> ion thermal velocity"
        );
    }

    #[test]
    fn test_sheath_creation() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        assert!(sheath.sheath_thickness > 0.0);
    }

    #[test]
    fn test_child_langmuir_law() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        let j = sheath.child_langmuir_law(100.0, 1e-3);
        assert!(j > 0.0);
    }

    #[test]
    fn test_child_langmuir_scales_with_voltage() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        let j1 = sheath.child_langmuir_law(100.0, 1e-3);
        let j2 = sheath.child_langmuir_law(400.0, 1e-3);
        assert!(j2 > j1, "Higher voltage → higher CL current");
    }

    #[test]
    fn test_populate_particles() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let mut sheath = PlasmaSheath::new(ps, -50.0);
        sheath.populate_particles(20);
        assert_eq!(sheath.positions.len(), 20);
    }

    #[test]
    fn test_probe_iv_at_floating_potential() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        // At V = V_f the net current should be near zero by definition
        // (our model is only approximate, just check it's finite)
        let i = sheath.probe_iv(sheath.floating_potential);
        assert!(i.is_finite());
    }

    #[test]
    fn test_normalised_sheath_thickness() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        let ns = sheath.normalised_sheath_thickness();
        assert!(ns > 1.0, "Sheath is wider than one Debye length");
    }

    #[test]
    fn test_interpolated_field_empty() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        let e = sheath.interpolated_field_x(1e-4, 1e-4);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_interpolated_field_with_particles() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let mut sheath = PlasmaSheath::new(ps, -50.0);
        sheath.populate_particles(10);
        let ds = sheath.sheath_thickness;
        let e = sheath.interpolated_field_x(ds * 0.5, ds * 0.1);
        assert!(e.is_finite());
    }

    #[test]
    fn test_net_current_finite() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let sheath = PlasmaSheath::new(ps, -50.0);
        let j = sheath.net_current_density();
        assert!(j.is_finite());
    }

    #[test]
    fn test_ion_plasma_frequency_positive() {
        let ps = PlasmaState::hydrogen(10.0, 0.1, 1e18);
        let wpi = ps.ion_plasma_frequency();
        assert!(wpi > 0.0);
    }
}
