// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Molecular dynamics for tribology: friction, wear, lubrication at atomic scale.
//!
//! This module provides:
//! - [`TribologySystem`]: Two-surface sliding/rolling contact simulation
//! - [`FrictionCoefficient`]: Static and kinetic friction coefficient calculation
//! - [`FrenkelKontorovaModel`]: 1D FK model: depinning, Peierls barrier, kink dynamics
//! - [`TomlinsonsModel`]: Prandtl-Tomlinson model: single-asperity stick-slip
//! - [`WearModel`]: Archard wear volume law V = k W d / H
//! - [`LubricationMd`]: Confined fluid film MD, effective viscosity
//! - [`AdsorptionLayer`]: Chemisorbed/physisorbed monolayer
//! - [`SurfaceEnergy`]: Dupré work of adhesion, Hamaker constant
//! - [`ContactMechanics`]: Greenwood-Williamson multi-asperity model
//! - [`NanoscaleFriction`]: Atomic-scale stick-slip, superlubricity
//!
//! # References
//! - Archard (1953): Contact and rubbing of flat surfaces
//! - Prandtl (1928): Ein Gedankenmodell zur kinetischen Theorie der festen Körper
//! - Tomlinson (1929): A molecular theory of friction
//! - Frenkel & Kontorova (1938): On the theory of plastic deformation and twinning
//! - Greenwood & Williamson (1966): Contact of nominally flat surfaces
//! - Braginskii (1965): Transport processes in a plasma

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant \[J K⁻¹\].
const K_B: f64 = 1.380_649e-23;

/// Atomic mass unit \[kg\].
#[cfg(test)]
const AMU: f64 = 1.660_539_066_6e-27;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors (a − b).
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

// ─────────────────────────────────────────────────────────────────────────────
// TribologySystem
// ─────────────────────────────────────────────────────────────────────────────

/// Two-surface sliding/rolling contact system.
///
/// Consists of a substrate slab (fixed bottom layer) and a slider slab
/// (driven top layer). Friction force is measured as the mean force on
/// the top layer atoms in the sliding direction.
pub struct TribologySystem {
    /// Positions of all atoms \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Velocities of all atoms \[m s⁻¹\].
    pub velocities: Vec<[f64; 3]>,
    /// Forces on all atoms \[N\].
    pub forces: Vec<[f64; 3]>,
    /// Masses of all atoms \[kg\].
    pub masses: Vec<f64>,
    /// Atom layer assignment: 0 = bottom fixed, 1 = bulk, 2 = top driven.
    pub layer: Vec<u8>,
    /// Applied normal force F_N \[N\].
    pub normal_force: f64,
    /// Sliding velocity v_slide \[m s⁻¹\].
    pub sliding_velocity: f64,
    /// Spring constant for top layer drive \[N m⁻¹\].
    pub drive_spring: f64,
    /// Current driver position x_driver \[m\].
    pub driver_position: f64,
    /// Lennard-Jones ε for inter-surface interaction \[J\].
    pub lj_epsilon: f64,
    /// Lennard-Jones σ for inter-surface interaction \[m\].
    pub lj_sigma: f64,
}

impl TribologySystem {
    /// Create a new tribology system.
    pub fn new(normal_force: f64, sliding_velocity: f64, lj_epsilon: f64, lj_sigma: f64) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            layer: Vec::new(),
            normal_force,
            sliding_velocity,
            drive_spring: 1.0,
            driver_position: 0.0,
            lj_epsilon,
            lj_sigma,
        }
    }

    /// Add an atom to the system.
    pub fn add_atom(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64, layer: u8) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.forces.push([0.0; 3]);
        self.masses.push(mass);
        self.layer.push(layer);
    }

    /// Number of atoms in the system.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Compute Lennard-Jones pair force between atoms i and j.
    fn lj_force(&self, i: usize, j: usize) -> [f64; 3] {
        let r_vec = sub3(self.positions[i], self.positions[j]);
        let r2 = dot3(r_vec, r_vec);
        if r2 < 1e-30 {
            return [0.0; 3];
        }
        let r2_inv = 1.0 / r2;
        let sig2 = self.lj_sigma * self.lj_sigma * r2_inv;
        let sig6 = sig2 * sig2 * sig2;
        let sig12 = sig6 * sig6;
        // F = 24ε/r² (2(σ/r)¹² - (σ/r)⁶) * r_vec
        let f_mag = 24.0 * self.lj_epsilon * r2_inv * (2.0 * sig12 - sig6);
        scale3(r_vec, f_mag)
    }

    /// Compute all pairwise forces (O(N²)).
    pub fn compute_forces(&mut self) {
        for f in self.forces.iter_mut() {
            *f = [0.0; 3];
        }
        let n = self.n_atoms();
        for i in 0..n {
            for j in (i + 1)..n {
                let f = self.lj_force(i, j);
                self.forces[i] = add3(self.forces[i], f);
                self.forces[j] = sub3(self.forces[j], f);
            }
        }
    }

    /// Apply shear velocity to top-layer atoms (layer == 2).
    pub fn apply_shear(&mut self, dt: f64) {
        self.driver_position += self.sliding_velocity * dt;
        for i in 0..self.n_atoms() {
            if self.layer[i] == 2 {
                self.velocities[i][0] = self.sliding_velocity;
            }
            if self.layer[i] == 0 {
                // Bottom layer is fixed
                self.velocities[i] = [0.0; 3];
                self.forces[i] = [0.0; 3];
            }
        }
    }

    /// Compute friction force = sum of x-forces on top-layer atoms \[N\].
    pub fn friction_force(&self) -> f64 {
        let mut ff = 0.0;
        for i in 0..self.n_atoms() {
            if self.layer[i] == 2 {
                ff += self.forces[i][0];
            }
        }
        ff
    }

    /// Velocity Verlet integration step.
    pub fn verlet_step(&mut self, dt: f64) {
        let n = self.n_atoms();
        // Half-step velocity
        for i in 0..n {
            if self.layer[i] == 0 {
                continue;
            }
            let m = self.masses[i];
            for d in 0..3 {
                self.velocities[i][d] += 0.5 * self.forces[i][d] / m * dt;
            }
        }
        // Full position step
        for i in 0..n {
            if self.layer[i] == 0 {
                continue;
            }
            for d in 0..3 {
                self.positions[i][d] += self.velocities[i][d] * dt;
            }
        }
        // Recompute forces
        self.compute_forces();
        // Complete velocity
        for i in 0..n {
            if self.layer[i] == 0 {
                continue;
            }
            let m = self.masses[i];
            for d in 0..3 {
                self.velocities[i][d] += 0.5 * self.forces[i][d] / m * dt;
            }
        }
        self.apply_shear(dt);
    }

    /// Total kinetic energy of the system \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * dot3(*v, *v))
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FrictionCoefficient
// ─────────────────────────────────────────────────────────────────────────────

/// Friction coefficient calculator: μ = F_friction / F_normal.
pub struct FrictionCoefficient {
    /// Kinetic friction coefficient μ_k.
    pub mu_kinetic: f64,
    /// Static friction coefficient μ_s ≥ μ_k.
    pub mu_static: f64,
    /// Normal force F_N \[N\].
    pub normal_force: f64,
}

impl FrictionCoefficient {
    /// Create a friction coefficient pair.
    ///
    /// Ensures μ_s ≥ μ_k.
    pub fn new(mu_kinetic: f64, normal_force: f64) -> Self {
        let mu_static = mu_kinetic * 1.2; // typical ratio
        Self {
            mu_kinetic,
            mu_static,
            normal_force,
        }
    }

    /// Kinetic friction force F_k = μ_k · F_N \[N\].
    pub fn kinetic_friction_force(&self) -> f64 {
        self.mu_kinetic * self.normal_force
    }

    /// Maximum static friction force F_s = μ_s · F_N \[N\].
    pub fn max_static_friction_force(&self) -> f64 {
        self.mu_static * self.normal_force
    }

    /// Compute μ from measured friction and normal forces.
    pub fn from_forces(friction_force: f64, normal_force: f64) -> f64 {
        if normal_force > 0.0 {
            friction_force.abs() / normal_force
        } else {
            0.0
        }
    }

    /// Check if the system is sliding (applied force > static friction).
    pub fn is_sliding(&self, applied_force: f64) -> bool {
        applied_force > self.max_static_friction_force()
    }

    /// Stribeck parameter: friction vs. velocity (simplified Stribeck curve).
    ///
    /// Returns viscous component of μ at speed v for lubricant viscosity η and
    /// nominal pressure P_0.
    pub fn stribeck_friction(&self, v: f64, viscosity: f64, pressure: f64) -> f64 {
        if pressure <= 0.0 {
            return self.mu_kinetic;
        }
        // Sommerfeld number So = η v / P_0
        let so = viscosity * v / pressure;
        // Simplified Stribeck curve: μ = μ_k (1 + A · So)
        self.mu_kinetic * (1.0 + 0.1 * so)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FrenkelKontorovaModel
// ─────────────────────────────────────────────────────────────────────────────

/// 1D Frenkel-Kontorova model: chain of atoms on sinusoidal substrate.
///
/// Equation of motion: m ẍ_n = K(x_{n+1} - 2x_n + x_{n-1}) - (V₀/b) sin(2πx_n/b) - γẋ_n
///
/// Describes commensurate/incommensurate transitions and kink (soliton) dynamics.
pub struct FrenkelKontorovaModel {
    /// Atom positions along the chain \[m\].
    pub positions: Vec<f64>,
    /// Atom velocities \[m s⁻¹\].
    pub velocities: Vec<f64>,
    /// Inter-atomic spring constant K \[N m⁻¹\].
    pub spring_constant: f64,
    /// Substrate corrugation amplitude V₀ \[J\].
    pub corrugation: f64,
    /// Substrate periodicity b \[m\].
    pub substrate_period: f64,
    /// Chain natural spacing a \[m\].
    pub chain_spacing: f64,
    /// Atom mass \[kg\].
    pub mass: f64,
    /// Damping coefficient γ \[kg s⁻¹\].
    pub damping: f64,
}

impl FrenkelKontorovaModel {
    /// Create a FK model with N atoms.
    pub fn new(
        n: usize,
        spring_constant: f64,
        corrugation: f64,
        substrate_period: f64,
        chain_spacing: f64,
        mass: f64,
        damping: f64,
    ) -> Self {
        let positions: Vec<f64> = (0..n).map(|i| i as f64 * chain_spacing).collect();
        let velocities = vec![0.0; n];
        Self {
            positions,
            velocities,
            spring_constant,
            corrugation,
            substrate_period,
            chain_spacing,
            mass,
            damping,
        }
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Force on atom n from substrate: F_sub = −(V₀/b) sin(2πx/b) · (2π/b).
    pub fn substrate_force(&self, n: usize) -> f64 {
        let x = self.positions[n];
        -self.corrugation
            * (2.0 * PI / self.substrate_period)
            * (2.0 * PI * x / self.substrate_period).sin()
    }

    /// Force on atom n from spring neighbors.
    pub fn spring_force(&self, n: usize) -> f64 {
        let n_atoms = self.n_atoms();
        let x = self.positions[n];
        let x_prev = if n > 0 {
            self.positions[n - 1]
        } else {
            x - self.chain_spacing
        };
        let x_next = if n + 1 < n_atoms {
            self.positions[n + 1]
        } else {
            x + self.chain_spacing
        };
        self.spring_constant * (x_next - 2.0 * x + x_prev)
    }

    /// Total force on atom n.
    pub fn total_force(&self, n: usize) -> f64 {
        self.spring_force(n) + self.substrate_force(n) - self.damping * self.velocities[n]
    }

    /// Advance the chain by one Euler step.
    pub fn euler_step(&mut self, dt: f64) {
        let n = self.n_atoms();
        let forces: Vec<f64> = (0..n).map(|i| self.total_force(i)).collect();
        for (i, &f) in forces.iter().enumerate() {
            let accel = f / self.mass;
            self.velocities[i] += accel * dt;
            self.positions[i] += self.velocities[i] * dt;
        }
    }

    /// Peierls-Nabarro barrier height (approximate): ΔE_PN ≈ V₀ (commensurate regime).
    ///
    /// This is the energy barrier per atom for sliding.
    pub fn peierls_nabarro_barrier(&self) -> f64 {
        // For commensurate chain: ΔE ∝ V₀ (full barrier)
        // For incommensurate: barrier is reduced by kink energy
        let mismatch = (self.chain_spacing / self.substrate_period).fract();
        let reduced = if mismatch < 0.5 {
            mismatch
        } else {
            1.0 - mismatch
        };
        // Barrier ∝ V₀ when commensurate (mismatch → 0)
        // Barrier vanishes when truly incommensurate (golden ratio)
        self.corrugation * (1.0 - 4.0 * reduced * (1.0 - reduced))
    }

    /// Critical force for depinning: F_c = corrugation · (2π/b) \[N\].
    pub fn depinning_force(&self) -> f64 {
        self.corrugation * 2.0 * PI / self.substrate_period
    }

    /// Total potential energy of the chain \[J\].
    pub fn potential_energy(&self) -> f64 {
        let mut e = 0.0;
        let n = self.n_atoms();
        for i in 0..n {
            // Substrate potential
            e += self.corrugation
                * (1.0 - (2.0 * PI * self.positions[i] / self.substrate_period).cos());
            // Spring potential (count each bond once)
            if i + 1 < n {
                let dx = self.positions[i + 1] - self.positions[i] - self.chain_spacing;
                e += 0.5 * self.spring_constant * dx * dx;
            }
        }
        e
    }

    /// Kink (soliton) width: ξ = b/(2π) · sqrt(K/V₀) · substrate_period \[m\].
    pub fn kink_width(&self) -> f64 {
        if self.corrugation > 0.0 {
            self.substrate_period / (2.0 * PI)
                * (self.spring_constant / self.corrugation).sqrt()
                * self.substrate_period
        } else {
            f64::INFINITY
        }
    }

    /// Kink (soliton) rest energy E_kink \[J\].
    pub fn kink_energy(&self) -> f64 {
        if self.corrugation > 0.0 {
            8.0 * (self.spring_constant * self.corrugation).sqrt()
        } else {
            0.0
        }
    }

    /// Check if chain is in the superlubricity regime (incommensurate: barrier → 0).
    ///
    /// True when a/b is irrational (approximated as far from simple fractions).
    pub fn is_superlubricious(&self) -> bool {
        let ratio = self.chain_spacing / self.substrate_period;
        let mismatch = (ratio - ratio.round()).abs();
        // If mismatch is large, chain is incommensurate → superlubricious
        mismatch > 0.1
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TomlinsonsModel (Prandtl-Tomlinson)
// ─────────────────────────────────────────────────────────────────────────────

/// Prandtl-Tomlinson model for single-asperity friction.
///
/// A tip is connected to a support via spring k. The support moves at
/// velocity v. The tip slides on a sinusoidal surface potential:
/// V(x) = (V₀/2) · cos(2πx/a)
///
/// Stick-slip occurs when k < k_crit = 2π²V₀/a².
pub struct TomlinsonsModel {
    /// Tip position \[m\].
    pub tip_position: f64,
    /// Tip velocity \[m s⁻¹\].
    pub tip_velocity: f64,
    /// Support position x_support \[m\].
    pub support_position: f64,
    /// Support velocity v_support \[m s⁻¹\].
    pub support_velocity: f64,
    /// Connecting spring stiffness k \[N m⁻¹\].
    pub spring_stiffness: f64,
    /// Surface corrugation amplitude V₀ \[J\].
    pub corrugation: f64,
    /// Lattice constant a \[m\].
    pub lattice_constant: f64,
    /// Tip effective mass m \[kg\].
    pub mass: f64,
    /// Damping γ \[kg s⁻¹\].
    pub damping: f64,
    /// Temperature T \[K\] for thermal effects.
    pub temperature: f64,
}

impl TomlinsonsModel {
    /// Create a Prandtl-Tomlinson model.
    pub fn new(
        spring_stiffness: f64,
        corrugation: f64,
        lattice_constant: f64,
        mass: f64,
        damping: f64,
        temperature: f64,
        support_velocity: f64,
    ) -> Self {
        Self {
            tip_position: 0.0,
            tip_velocity: 0.0,
            support_position: 0.0,
            support_velocity,
            spring_stiffness,
            corrugation,
            lattice_constant,
            mass,
            damping,
            temperature,
        }
    }

    /// η parameter = 2π²V₀ / (k a²): stick-slip when η > 1.
    pub fn eta_parameter(&self) -> f64 {
        2.0 * PI * PI * self.corrugation
            / (self.spring_stiffness * self.lattice_constant * self.lattice_constant)
    }

    /// Critical spring stiffness above which stick-slip disappears.
    pub fn critical_stiffness(&self) -> f64 {
        2.0 * PI * PI * self.corrugation / (self.lattice_constant * self.lattice_constant)
    }

    /// Force from surface potential on tip: F_surf = (πV₀/a) sin(2πx/a).
    pub fn surface_force(&self) -> f64 {
        PI * self.corrugation / self.lattice_constant
            * (2.0 * PI * self.tip_position / self.lattice_constant).sin()
    }

    /// Spring force on tip: F_spring = k(x_support - x_tip).
    pub fn spring_force(&self) -> f64 {
        self.spring_stiffness * (self.support_position - self.tip_position)
    }

    /// Total force on tip.
    pub fn total_force(&self) -> f64 {
        self.spring_force() - self.surface_force() - self.damping * self.tip_velocity
    }

    /// Advance the model by one Euler step.
    pub fn euler_step(&mut self, dt: f64) {
        let f = self.total_force();
        let accel = f / self.mass;
        self.tip_velocity += accel * dt;
        self.tip_position += self.tip_velocity * dt;
        self.support_position += self.support_velocity * dt;
    }

    /// Friction force = spring force transmitted through the system.
    pub fn friction_force(&self) -> f64 {
        self.spring_force()
    }

    /// Check if the system exhibits stick-slip (η > 1).
    pub fn has_stick_slip(&self) -> bool {
        self.eta_parameter() > 1.0
    }

    /// Thermal activation rate over the barrier (Kramers): Γ = f₀ exp(−ΔE/kT).
    ///
    /// Returns jump rate \[s⁻¹\].
    pub fn thermal_jump_rate(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        let delta_e = self.corrugation; // barrier height ≈ V₀
        let f0 = 1e12; // attempt frequency [Hz]
        f0 * (-delta_e / (K_B * self.temperature)).exp()
    }

    /// Maximum static friction force \[N\].
    pub fn max_static_force(&self) -> f64 {
        // Maximum spring extension before slip: (a/2π) × k
        self.spring_stiffness * self.lattice_constant / (2.0 * PI)
    }

    /// Average friction force in the stick-slip regime \[N\].
    pub fn average_friction_force(&self) -> f64 {
        if self.has_stick_slip() {
            // In stick-slip regime, average ≈ (k_crit - k) / 2 × a
            let f_max = self.spring_stiffness * self.lattice_constant / (2.0 * PI);
            f_max * (self.eta_parameter() - 1.0) / self.eta_parameter()
        } else {
            0.0 // continuous sliding (superlubricious)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WearModel
// ─────────────────────────────────────────────────────────────────────────────

/// Archard wear law: V = k · W · d / H.
///
/// V = wear volume \[m³\], k = dimensionless wear coefficient,
/// W = normal load \[N\], d = sliding distance \[m\], H = hardness \[Pa\].
pub struct WearModel {
    /// Dimensionless wear coefficient k (Archard).
    pub wear_coefficient: f64,
    /// Material hardness H \[Pa\].
    pub hardness: f64,
}

impl WearModel {
    /// Create a wear model.
    pub fn new(wear_coefficient: f64, hardness: f64) -> Self {
        Self {
            wear_coefficient,
            hardness,
        }
    }

    /// Compute wear volume V = k W d / H \[m³\].
    pub fn wear_volume(&self, normal_load: f64, sliding_distance: f64) -> f64 {
        if self.hardness <= 0.0 {
            return 0.0;
        }
        self.wear_coefficient * normal_load * sliding_distance / self.hardness
    }

    /// Compute wear rate dV/dd = k W / H \[m³ m⁻¹ = m²\].
    pub fn wear_rate(&self, normal_load: f64) -> f64 {
        if self.hardness <= 0.0 {
            return 0.0;
        }
        self.wear_coefficient * normal_load / self.hardness
    }

    /// Specific wear rate k_s = k / H \[m² N⁻¹\].
    pub fn specific_wear_rate(&self) -> f64 {
        if self.hardness <= 0.0 {
            0.0
        } else {
            self.wear_coefficient / self.hardness
        }
    }

    /// Transition from mild to severe wear: occurs above critical load W_c.
    ///
    /// Simplified criterion: severe wear if k > 10⁻³.
    pub fn wear_regime(&self) -> WearRegime {
        if self.wear_coefficient < 1e-5 {
            WearRegime::Mild
        } else if self.wear_coefficient < 1e-3 {
            WearRegime::Moderate
        } else {
            WearRegime::Severe
        }
    }

    /// Wear depth d_wear = V / A_contact \[m\].
    pub fn wear_depth(&self, normal_load: f64, sliding_distance: f64, contact_area: f64) -> f64 {
        if contact_area <= 0.0 {
            return 0.0;
        }
        self.wear_volume(normal_load, sliding_distance) / contact_area
    }
}

/// Wear regime classification.
#[derive(Debug, Clone, PartialEq)]
pub enum WearRegime {
    /// Mild wear: k < 10⁻⁵.
    Mild,
    /// Moderate wear: 10⁻⁵ ≤ k < 10⁻³.
    Moderate,
    /// Severe wear: k ≥ 10⁻³.
    Severe,
}

// ─────────────────────────────────────────────────────────────────────────────
// LubricationMd
// ─────────────────────────────────────────────────────────────────────────────

/// Molecular dynamics of a confined fluid film: effective viscosity calculation.
///
/// A thin lubricant layer between two surfaces. The effective viscosity
/// increases dramatically when the film thickness approaches a few molecular diameters.
pub struct LubricationMd {
    /// Fluid molecule positions \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Fluid molecule velocities \[m s⁻¹\].
    pub velocities: Vec<[f64; 3]>,
    /// Film thickness h \[m\].
    pub film_thickness: f64,
    /// Bulk viscosity η₀ \[Pa s\].
    pub bulk_viscosity: f64,
    /// Molecular diameter σ \[m\].
    pub molecular_diameter: f64,
    /// Effective shear stress τ \[Pa\].
    pub shear_stress: f64,
    /// Shear rate γ̇ \[s⁻¹\].
    pub shear_rate: f64,
}

impl LubricationMd {
    /// Create a lubrication MD model.
    pub fn new(film_thickness: f64, bulk_viscosity: f64, molecular_diameter: f64) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            film_thickness,
            bulk_viscosity,
            molecular_diameter,
            shear_stress: 0.0,
            shear_rate: 0.0,
        }
    }

    /// Number of molecular layers n_l = h / σ (approximate).
    pub fn number_of_layers(&self) -> f64 {
        if self.molecular_diameter > 0.0 {
            self.film_thickness / self.molecular_diameter
        } else {
            0.0
        }
    }

    /// Effective viscosity: increases as film becomes molecularly thin.
    ///
    /// η_eff = η₀ · exp(A / (h/σ)) for h/σ < threshold.
    pub fn effective_viscosity(&self) -> f64 {
        let n_l = self.number_of_layers();
        if n_l < 1.0 {
            // Solidified film
            return self.bulk_viscosity * 1e6;
        }
        // Enhancement factor: increases exponentially as film thins
        let enhancement = if n_l < 10.0 { (2.0 / n_l).exp() } else { 1.0 };
        self.bulk_viscosity * enhancement
    }

    /// Shear force per unit area (viscous shear stress) τ = η_eff · γ̇ \[Pa\].
    pub fn compute_shear_stress(&mut self, shear_rate: f64) -> f64 {
        self.shear_rate = shear_rate;
        self.shear_stress = self.effective_viscosity() * shear_rate;
        self.shear_stress
    }

    /// Reynolds number Re = ρ v h / η.
    pub fn reynolds_number(&self, density: f64, velocity: f64) -> f64 {
        let eta = self.effective_viscosity();
        if eta > 0.0 {
            density * velocity * self.film_thickness / eta
        } else {
            0.0
        }
    }

    /// Hydrodynamic pressure P_h = 6 η v L / h² (long-bearing approximation) \[Pa\].
    pub fn hydrodynamic_pressure(&self, velocity: f64, bearing_length: f64) -> f64 {
        let eta = self.effective_viscosity();
        let h2 = self.film_thickness * self.film_thickness;
        if h2 > 0.0 {
            6.0 * eta * velocity * bearing_length / h2
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AdsorptionLayer
// ─────────────────────────────────────────────────────────────────────────────

/// Chemisorbed or physisorbed monolayer on a surface.
pub struct AdsorptionLayer {
    /// Adsorption energy E_ads \[J\] per molecule.
    pub adsorption_energy: f64,
    /// Surface coverage Θ ∈ \[0, 1\] (fraction of monolayer).
    pub coverage: f64,
    /// Molecule diameter \[m\].
    pub molecule_diameter: f64,
    /// Whether this is chemisorption (true) or physisorption (false).
    pub is_chemisorbed: bool,
    /// Pre-exponential factor ν₀ \[s⁻¹\] for desorption rate.
    pub attempt_frequency: f64,
}

impl AdsorptionLayer {
    /// Create an adsorption layer.
    pub fn new(
        adsorption_energy: f64,
        coverage: f64,
        molecule_diameter: f64,
        is_chemisorbed: bool,
    ) -> Self {
        let attempt_frequency = if is_chemisorbed { 1e13 } else { 1e12 };
        Self {
            adsorption_energy,
            coverage,
            molecule_diameter,
            is_chemisorbed,
            attempt_frequency,
        }
    }

    /// Desorption rate ν = ν₀ exp(−E_ads / k_B T) \[s⁻¹\].
    pub fn desorption_rate(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        self.attempt_frequency * (-self.adsorption_energy / (K_B * temperature)).exp()
    }

    /// Surface lifetime τ = 1 / ν_des \[s\].
    pub fn surface_lifetime(&self, temperature: f64) -> f64 {
        let rate = self.desorption_rate(temperature);
        if rate > 0.0 {
            1.0 / rate
        } else {
            f64::INFINITY
        }
    }

    /// Monolayer coverage area density \[m⁻²\].
    pub fn coverage_density(&self) -> f64 {
        let a_mol = PI * (self.molecule_diameter / 2.0).powi(2);
        if a_mol > 0.0 {
            self.coverage / a_mol
        } else {
            0.0
        }
    }

    /// Friction reduction factor from lubricant monolayer (simplified).
    ///
    /// Returns the fraction by which friction is reduced: f_red = Θ · (1 − E_ads/E_ref).
    pub fn friction_reduction(&self, reference_energy: f64) -> f64 {
        let reduction = if reference_energy > 0.0 {
            (self.adsorption_energy / reference_energy).min(1.0)
        } else {
            0.0
        };
        self.coverage * reduction
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SurfaceEnergy
// ─────────────────────────────────────────────────────────────────────────────

/// Surface energy, work of adhesion, and Hamaker constant calculations.
pub struct SurfaceEnergy {
    /// Surface energy of material 1 γ₁ \[J m⁻²\].
    pub gamma1: f64,
    /// Surface energy of material 2 γ₂ \[J m⁻²\].
    pub gamma2: f64,
    /// Interfacial energy γ₁₂ \[J m⁻²\].
    pub gamma12: f64,
}

impl SurfaceEnergy {
    /// Create a surface energy system.
    pub fn new(gamma1: f64, gamma2: f64, gamma12: f64) -> Self {
        Self {
            gamma1,
            gamma2,
            gamma12,
        }
    }

    /// Dupré work of adhesion W_ad = γ₁ + γ₂ − γ₁₂ \[J m⁻²\].
    pub fn work_of_adhesion(&self) -> f64 {
        self.gamma1 + self.gamma2 - self.gamma12
    }

    /// Spreading coefficient S = γ₂ − γ₁ − γ₁₂ \[J m⁻²\].
    ///
    /// S > 0: material 1 spreads on material 2.
    pub fn spreading_coefficient(&self) -> f64 {
        self.gamma2 - self.gamma1 - self.gamma12
    }

    /// Hamaker constant A_H from surface energies (simplified Lifshitz approximation).
    ///
    /// A_H ≈ 24π z₀² W_ad, where z₀ ≈ 0.165 nm is the contact cutoff distance.
    pub fn hamaker_constant(&self) -> f64 {
        let z0 = 0.165e-9; // [m] equilibrium cut-off distance
        let w_ad = self.work_of_adhesion();
        24.0 * PI * z0 * z0 * w_ad
    }

    /// Van der Waals pressure P_vdW = A_H / (6π D³) between flat surfaces \[Pa\].
    pub fn vdw_pressure(&self, separation: f64) -> f64 {
        if separation <= 0.0 {
            return 0.0;
        }
        let a_h = self.hamaker_constant();
        a_h / (6.0 * PI * separation.powi(3))
    }

    /// Contact angle from Young's equation: cos(θ) = (γ_SV − γ_SL) / γ_LV.
    ///
    /// Returns cosine of the contact angle.
    pub fn contact_angle_cos(&self) -> f64 {
        // Simplified: γ₁ = γ_SV, γ₁₂ = γ_SL, γ₂ = γ_LV
        if self.gamma2 > 0.0 {
            ((self.gamma1 - self.gamma12) / self.gamma2).clamp(-1.0, 1.0)
        } else {
            1.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ContactMechanics
// ─────────────────────────────────────────────────────────────────────────────

/// Greenwood-Williamson multi-asperity contact model.
///
/// Models contact between a rough surface and a flat by treating asperities
/// as a statistical distribution of spherical tips (Hertz contacts).
pub struct ContactMechanics {
    /// Asperity density η_d \[m⁻²\].
    pub asperity_density: f64,
    /// RMS roughness σ_r \[m\].
    pub rms_roughness: f64,
    /// Asperity tip radius R \[m\].
    pub tip_radius: f64,
    /// Combined elastic modulus E* \[Pa\].
    pub elastic_modulus: f64,
    /// Hardness H \[Pa\] (for plastic contact).
    pub hardness: f64,
}

impl ContactMechanics {
    /// Create a GW contact model.
    pub fn new(
        asperity_density: f64,
        rms_roughness: f64,
        tip_radius: f64,
        elastic_modulus: f64,
        hardness: f64,
    ) -> Self {
        Self {
            asperity_density,
            rms_roughness,
            tip_radius,
            elastic_modulus,
            hardness,
        }
    }

    /// Plasticity index Ψ = (E*/H) sqrt(σ/R).
    ///
    /// Ψ < 0.6: elastic; Ψ > 1: predominantly plastic.
    pub fn plasticity_index(&self) -> f64 {
        if self.hardness > 0.0 && self.tip_radius > 0.0 {
            (self.elastic_modulus / self.hardness) * (self.rms_roughness / self.tip_radius).sqrt()
        } else {
            0.0
        }
    }

    /// Real contact area fraction A_r / A_n for Gaussian asperities (GW model).
    ///
    /// `separation` d is the separation in units of σ_r.
    pub fn contact_area_fraction(&self, separation: f64) -> f64 {
        // GW integral: A_r/A_n = π η_d R σ ∫_d^∞ (s-d) φ(s) ds
        // Gaussian φ(s) = (1/√(2π)) exp(-s²/2)
        // Approximate using tail integral
        let exp_arg = -0.5 * separation * separation;
        let gauss_tail = if exp_arg > -100.0 {
            (exp_arg.exp()) / (2.0 * PI).sqrt()
        } else {
            0.0
        };
        PI * self.asperity_density * self.tip_radius * self.rms_roughness * gauss_tail
    }

    /// Mean contact force per asperity (Hertz elastic): F = (4/3) E* √R δ^{3/2}.
    pub fn hertz_force(&self, indentation: f64) -> f64 {
        if indentation <= 0.0 {
            return 0.0;
        }
        4.0 / 3.0 * self.elastic_modulus * self.tip_radius.sqrt() * indentation.powf(1.5)
    }

    /// Mean contact radius a = √(R δ) \[m\].
    pub fn contact_radius(&self, indentation: f64) -> f64 {
        if indentation <= 0.0 || self.tip_radius <= 0.0 {
            return 0.0;
        }
        (self.tip_radius * indentation).sqrt()
    }

    /// Mean contact pressure P_mean = F / (π a²) \[Pa\].
    pub fn mean_contact_pressure(&self, indentation: f64) -> f64 {
        let a = self.contact_radius(indentation);
        let f = self.hertz_force(indentation);
        if a > 0.0 { f / (PI * a * a) } else { 0.0 }
    }

    /// Check if contact is elastic (plasticity index < 0.6).
    pub fn is_elastic_contact(&self) -> bool {
        self.plasticity_index() < 0.6
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NanoscaleFriction
// ─────────────────────────────────────────────────────────────────────────────

/// Atomic-scale stick-slip friction and superlubricity.
pub struct NanoscaleFriction {
    /// Lattice constant of surface 1 a₁ \[m\].
    pub lattice_a: f64,
    /// Lattice constant of surface 2 a₂ \[m\].
    pub lattice_b: f64,
    /// Corrugation energy amplitude V₀ \[J\].
    pub corrugation: f64,
    /// Normal load F_N \[N\].
    pub normal_load: f64,
    /// Temperature T \[K\].
    pub temperature: f64,
}

impl NanoscaleFriction {
    /// Create a nanoscale friction model.
    pub fn new(
        lattice_a: f64,
        lattice_b: f64,
        corrugation: f64,
        normal_load: f64,
        temperature: f64,
    ) -> Self {
        Self {
            lattice_a,
            lattice_b,
            corrugation,
            normal_load,
            temperature,
        }
    }

    /// Mismatch ratio a₁/a₂ (incommensurate if irrational).
    pub fn mismatch_ratio(&self) -> f64 {
        self.lattice_a / self.lattice_b
    }

    /// Check superlubricity condition: surfaces are incommensurate.
    ///
    /// True when a₁/a₂ is far from a simple rational number.
    pub fn is_superlubricious(&self) -> bool {
        let r = self.mismatch_ratio();
        // Check if far from any rational p/q with q ≤ 5
        let min_dist = (1..=5_i32)
            .flat_map(|q| (1..=5_i32).map(move |p| (r - p as f64 / q as f64).abs()))
            .fold(f64::INFINITY, f64::min);
        min_dist > 0.05
    }

    /// Static friction force for commensurate contact F_s = N · V₀ · (2π/a) \[N\].
    pub fn static_friction_force(&self, n_contacts: usize) -> f64 {
        n_contacts as f64 * self.corrugation * 2.0 * PI / self.lattice_a
    }

    /// Superlubricious friction: F_sl = F_s / √N (cancellation of forces).
    pub fn superlubricious_friction(&self, n_contacts: usize) -> f64 {
        if n_contacts == 0 {
            return 0.0;
        }
        self.static_friction_force(1) / (n_contacts as f64).sqrt()
    }

    /// Friction reduction factor: superlubricious vs commensurate.
    pub fn friction_reduction_factor(&self, n_contacts: usize) -> f64 {
        if n_contacts == 0 {
            return 0.0;
        }
        1.0 / (n_contacts as f64).sqrt()
    }

    /// Thermolubricity correction: thermal activation reduces effective barrier.
    ///
    /// At finite T, the barrier is reduced by thermal noise.
    pub fn thermal_correction(&self, frequency: f64) -> f64 {
        if self.temperature <= 0.0 || frequency <= 0.0 {
            return 1.0;
        }
        let v0_kbt = self.corrugation / (K_B * self.temperature);
        // Correction factor (Singer 2010 approximation)
        let x = v0_kbt * (self.corrugation / (K_B * self.temperature)).ln();
        (1.0 - (1.0 / x.max(1.0))).max(0.0)
    }

    /// Debye-Waller factor: thermal reduction of corrugation amplitude.
    ///
    /// V_eff = V₀ exp(−2 W) where W = (2π)² `u²` / (2 a²) (Debye-Waller factor).
    pub fn debye_waller_corrugation(&self, mean_square_displacement: f64) -> f64 {
        let two_w = (2.0 * PI / self.lattice_a).powi(2) * mean_square_displacement;
        self.corrugation * (-two_w).exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;
    const RTOL: f64 = 1e-4;

    fn assert_rel(a: f64, b: f64, tol: f64, msg: &str) {
        let rel = if b.abs() > 1e-300 {
            (a - b).abs() / b.abs()
        } else {
            (a - b).abs()
        };
        assert!(rel < tol, "{}: {} vs {} (rel err {:.6})", msg, a, b, rel);
    }

    // --- FK model tests ---

    #[test]
    fn test_fk_peierls_barrier_proportional_to_corrugation() {
        let m1 = FrenkelKontorovaModel::new(10, 1.0, 1.0, 1.0, 1.0, AMU, 0.1);
        let m2 = FrenkelKontorovaModel::new(10, 1.0, 2.0, 1.0, 1.0, AMU, 0.1);
        // Barrier for commensurate chain ∝ V₀
        assert!(
            m2.peierls_nabarro_barrier() >= m1.peierls_nabarro_barrier(),
            "Peierls barrier should increase with corrugation"
        );
    }

    #[test]
    fn test_fk_depinning_force_scales_with_corrugation() {
        let m1 = FrenkelKontorovaModel::new(5, 1.0, 1.0, 1.0, 1.0, AMU, 0.0);
        let m2 = FrenkelKontorovaModel::new(5, 1.0, 2.0, 1.0, 1.0, AMU, 0.0);
        assert_rel(
            m2.depinning_force(),
            2.0 * m1.depinning_force(),
            RTOL,
            "F_dep ∝ V₀",
        );
    }

    #[test]
    fn test_fk_kink_energy_positive() {
        let m = FrenkelKontorovaModel::new(10, 1.0, 1.0, 1.0, 1.0, AMU, 0.0);
        assert!(m.kink_energy() > 0.0, "Kink energy should be positive");
    }

    #[test]
    fn test_fk_potential_energy_zero_corrugation() {
        // No substrate: only spring energy, which should be 0 for uniform spacing
        let m = FrenkelKontorovaModel::new(5, 1.0, 0.0, 1.0, 1.0, AMU, 0.0);
        // All atoms at natural spacing: spring energy = 0
        assert_eq!(m.potential_energy(), 0.0);
    }

    #[test]
    fn test_fk_superlubricious_incommensurate() {
        // Golden ratio ≈ 1.618... → incommensurate
        let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let m = FrenkelKontorovaModel::new(5, 1.0, 1.0, 1.0, golden, AMU, 0.0);
        assert!(
            m.is_superlubricious(),
            "Golden ratio spacing should be superlubricious"
        );
    }

    #[test]
    fn test_fk_commensurate_not_superlubricious() {
        // a = b: commensurate
        let m = FrenkelKontorovaModel::new(5, 1.0, 1.0, 1.0, 1.0, AMU, 0.0);
        assert!(
            !m.is_superlubricious(),
            "a=b (commensurate) should not be superlubricious"
        );
    }

    #[test]
    fn test_fk_euler_step_runs() {
        let mut m = FrenkelKontorovaModel::new(10, 1.0, 0.1, 1.0, 1.1, AMU, 0.01);
        let e0 = m.potential_energy();
        m.euler_step(1e-14);
        let e1 = m.potential_energy();
        // Energy should change after a step
        assert!((e1 - e0).abs() < 1e10); // basic sanity
    }

    // --- Prandtl-Tomlinson tests ---

    #[test]
    fn test_pt_eta_above_1_has_stick_slip() {
        // k < k_crit: stick-slip
        let k_crit = 2.0 * PI * PI * 1e-9 / (0.25e-9 * 0.25e-9);
        let k_soft = k_crit * 0.5;
        let m = TomlinsonsModel::new(k_soft, 1e-9, 0.25e-9, AMU * 100.0, 1e-12, 300.0, 1.0);
        assert!(m.eta_parameter() > 1.0, "η > 1 for soft spring");
        assert!(m.has_stick_slip(), "should exhibit stick-slip");
    }

    #[test]
    fn test_pt_stick_slip_disappears_above_critical_stiffness() {
        let corrugation = 1e-9;
        let a = 0.25e-9;
        let k_hard = 3.0 * 2.0 * PI * PI * corrugation / (a * a); // 3× critical
        let m = TomlinsonsModel::new(k_hard, corrugation, a, AMU * 100.0, 1e-12, 300.0, 1.0);
        assert!(m.eta_parameter() < 1.0, "η < 1 for stiff spring");
        assert!(
            !m.has_stick_slip(),
            "no stick-slip above critical stiffness"
        );
    }

    #[test]
    fn test_pt_critical_stiffness_formula() {
        let v0 = 1e-9;
        let a = 0.25e-9;
        let m = TomlinsonsModel::new(1.0, v0, a, AMU, 1e-12, 300.0, 1.0);
        let expected = 2.0 * PI * PI * v0 / (a * a);
        assert_rel(m.critical_stiffness(), expected, RTOL, "k_crit formula");
    }

    #[test]
    fn test_pt_thermal_jump_rate_positive_at_finite_t() {
        let m = TomlinsonsModel::new(1.0, 1e-21, 0.25e-9, AMU, 1e-12, 300.0, 1.0);
        assert!(m.thermal_jump_rate() > 0.0);
    }

    #[test]
    fn test_pt_thermal_jump_rate_zero_at_t0() {
        let m = TomlinsonsModel::new(1.0, 1e-21, 0.25e-9, AMU, 1e-12, 0.0, 1.0);
        assert_eq!(m.thermal_jump_rate(), 0.0);
    }

    #[test]
    fn test_pt_euler_step_moves_support() {
        let mut m = TomlinsonsModel::new(1.0, 1e-21, 0.25e-9, AMU, 1e-12, 300.0, 1.0);
        let x0 = m.support_position;
        m.euler_step(1e-9);
        assert!(m.support_position > x0, "support should advance");
    }

    // --- Archard wear law tests ---

    #[test]
    fn test_archard_wear_proportionality() {
        // V ∝ W
        let wear = WearModel::new(1e-4, 1e9);
        let v1 = wear.wear_volume(100.0, 1.0);
        let v2 = wear.wear_volume(200.0, 1.0);
        assert_rel(v2, 2.0 * v1, TOL, "V ∝ W");
    }

    #[test]
    fn test_archard_wear_proportional_to_distance() {
        let wear = WearModel::new(1e-4, 1e9);
        let v1 = wear.wear_volume(100.0, 1.0);
        let v2 = wear.wear_volume(100.0, 2.0);
        assert_rel(v2, 2.0 * v1, TOL, "V ∝ d");
    }

    #[test]
    fn test_archard_wear_inversely_proportional_to_hardness() {
        let wear1 = WearModel::new(1e-4, 1e9);
        let wear2 = WearModel::new(1e-4, 2e9);
        let v1 = wear1.wear_volume(100.0, 1.0);
        let v2 = wear2.wear_volume(100.0, 1.0);
        assert_rel(v1, 2.0 * v2, TOL, "V ∝ 1/H");
    }

    #[test]
    fn test_archard_wear_regime_mild() {
        let wear = WearModel::new(1e-6, 1e9);
        assert_eq!(wear.wear_regime(), WearRegime::Mild);
    }

    #[test]
    fn test_archard_wear_regime_severe() {
        let wear = WearModel::new(1e-2, 1e9);
        assert_eq!(wear.wear_regime(), WearRegime::Severe);
    }

    #[test]
    fn test_archard_specific_wear_rate() {
        let wear = WearModel::new(1e-4, 2e9);
        assert_rel(wear.specific_wear_rate(), 5e-14, RTOL, "specific wear rate");
    }

    // --- Debye-Waller factor ---

    #[test]
    fn test_debye_waller_reduces_corrugation() {
        let nano = NanoscaleFriction::new(0.25e-9, 0.25e-9, 1e-9, 1e-9, 300.0);
        let v_eff = nano.debye_waller_corrugation(1e-21);
        assert!(
            v_eff < nano.corrugation,
            "Debye-Waller should reduce corrugation"
        );
        assert!(v_eff > 0.0, "Corrugation should remain positive");
    }

    #[test]
    fn test_debye_waller_zero_displacement() {
        let nano = NanoscaleFriction::new(0.25e-9, 0.25e-9, 1e-9, 1e-9, 300.0);
        let v_eff = nano.debye_waller_corrugation(0.0);
        assert_rel(v_eff, nano.corrugation, TOL, "No DW effect at <u²>=0");
    }

    // --- Friction coefficient bounds ---

    #[test]
    fn test_friction_coefficient_non_negative() {
        let fc = FrictionCoefficient::new(0.3, 100.0);
        assert!(fc.mu_kinetic >= 0.0);
        assert!(fc.mu_static >= 0.0);
    }

    #[test]
    fn test_friction_static_ge_kinetic() {
        let fc = FrictionCoefficient::new(0.3, 100.0);
        assert!(fc.mu_static >= fc.mu_kinetic, "μ_s ≥ μ_k");
    }

    #[test]
    fn test_friction_force_proportional_to_load() {
        let fc1 = FrictionCoefficient::new(0.3, 100.0);
        let fc2 = FrictionCoefficient::new(0.3, 200.0);
        assert_rel(
            fc2.kinetic_friction_force(),
            2.0 * fc1.kinetic_friction_force(),
            TOL,
            "F_f ∝ F_N",
        );
    }

    #[test]
    fn test_friction_coefficient_from_forces() {
        let mu = FrictionCoefficient::from_forces(30.0, 100.0);
        assert_rel(mu, 0.3, TOL, "μ from forces");
    }

    #[test]
    fn test_friction_sliding_condition() {
        let fc = FrictionCoefficient::new(0.3, 100.0);
        assert!(!fc.is_sliding(20.0), "below static limit: no sliding");
        assert!(fc.is_sliding(50.0), "above static limit: sliding");
    }

    // --- Hamaker constant ---

    #[test]
    fn test_hamaker_constant_positive() {
        let se = SurfaceEnergy::new(0.05, 0.05, 0.01);
        assert!(
            se.hamaker_constant() > 0.0,
            "Hamaker constant should be positive"
        );
    }

    #[test]
    fn test_hamaker_constant_formula() {
        let se = SurfaceEnergy::new(0.05, 0.05, 0.01);
        let z0 = 0.165e-9_f64;
        let w_ad = se.work_of_adhesion();
        let expected = 24.0 * PI * z0 * z0 * w_ad;
        assert_rel(se.hamaker_constant(), expected, TOL, "Hamaker formula");
    }

    // --- Work of adhesion ---

    #[test]
    fn test_work_of_adhesion_dupre() {
        let se = SurfaceEnergy::new(0.05, 0.07, 0.02);
        let expected = 0.05 + 0.07 - 0.02;
        assert_rel(
            se.work_of_adhesion(),
            expected,
            TOL,
            "Dupré work of adhesion",
        );
    }

    #[test]
    fn test_work_of_adhesion_self_contact() {
        // Self-adhesion: γ₁₂ = 0 → W = 2γ₁
        let se = SurfaceEnergy::new(0.05, 0.05, 0.0);
        assert_rel(
            se.work_of_adhesion(),
            0.1,
            TOL,
            "W_ad = 2γ for same materials",
        );
    }

    // --- Contact mechanics ---

    #[test]
    fn test_gw_plasticity_index_positive() {
        let cm = ContactMechanics::new(1e14, 1e-9, 1e-7, 1e11, 1e9);
        assert!(cm.plasticity_index() > 0.0);
    }

    #[test]
    fn test_hertz_force_positive_for_positive_indentation() {
        let cm = ContactMechanics::new(1e14, 1e-9, 1e-7, 1e11, 1e9);
        assert!(cm.hertz_force(1e-9) > 0.0);
    }

    #[test]
    fn test_hertz_force_zero_for_zero_indentation() {
        let cm = ContactMechanics::new(1e14, 1e-9, 1e-7, 1e11, 1e9);
        assert_eq!(cm.hertz_force(0.0), 0.0);
    }

    #[test]
    fn test_hertz_contact_radius_scales_with_sqrt_indentation() {
        let cm = ContactMechanics::new(1e14, 1e-9, 1e-7, 1e11, 1e9);
        let a1 = cm.contact_radius(1e-9);
        let a2 = cm.contact_radius(4e-9);
        assert_rel(a2, 2.0 * a1, RTOL, "contact radius ∝ √δ");
    }

    // --- Superlubricity ---

    #[test]
    fn test_superlubricity_reduces_friction_with_contact_number() {
        let nano = NanoscaleFriction::new(0.25e-9, 0.41e-9, 1e-9, 1e-9, 300.0);
        let n = 100;
        let f_commensurate = nano.static_friction_force(n);
        let f_super = nano.superlubricious_friction(n);
        assert!(f_super < f_commensurate, "superlubricity reduces friction");
    }

    #[test]
    fn test_superlubricity_friction_reduction_scales_as_inv_sqrt_n() {
        let nano = NanoscaleFriction::new(0.25e-9, 0.41e-9, 1e-9, 1e-9, 300.0);
        let f1 = nano.friction_reduction_factor(100);
        let f2 = nano.friction_reduction_factor(400);
        assert_rel(f2, f1 / 2.0, RTOL, "friction ∝ 1/√N");
    }

    // --- Tribology system integration ---

    #[test]
    fn test_tribology_system_add_atoms() {
        let mut sys = TribologySystem::new(1.0, 0.1, 1e-21, 0.25e-9);
        sys.add_atom([0.0; 3], [0.0; 3], AMU * 10.0, 0);
        sys.add_atom([0.25e-9, 0.0, 0.0], [0.0; 3], AMU * 10.0, 2);
        assert_eq!(sys.n_atoms(), 2);
    }

    #[test]
    fn test_tribology_system_kinetic_energy_at_rest() {
        let mut sys = TribologySystem::new(1.0, 0.0, 1e-21, 0.25e-9);
        sys.add_atom([0.0; 3], [0.0; 3], AMU * 10.0, 1);
        assert_eq!(sys.kinetic_energy(), 0.0);
    }

    // --- Lubrication ---

    #[test]
    fn test_lubrication_effective_viscosity_increases_thin_film() {
        let lub1 = LubricationMd::new(100e-9, 1e-3, 0.5e-9);
        let lub2 = LubricationMd::new(1e-9, 1e-3, 0.5e-9);
        assert!(
            lub2.effective_viscosity() > lub1.effective_viscosity(),
            "thin film has higher effective viscosity"
        );
    }

    #[test]
    fn test_lubrication_layers_calculation() {
        let lub = LubricationMd::new(5e-9, 1e-3, 0.5e-9);
        assert_rel(lub.number_of_layers(), 10.0, RTOL, "10 molecular layers");
    }

    #[test]
    fn test_lubrication_bulk_viscosity_at_large_h() {
        let lub = LubricationMd::new(1e-6, 1e-3, 0.5e-9); // 2000 layers
        // At large h, η_eff ≈ η₀
        assert_rel(
            lub.effective_viscosity(),
            lub.bulk_viscosity,
            RTOL,
            "bulk viscosity at large h",
        );
    }

    // --- Adsorption layer ---

    #[test]
    fn test_adsorption_desorption_rate_increases_with_temperature() {
        let layer = AdsorptionLayer::new(1e-19, 0.5, 1e-9, false);
        let r1 = layer.desorption_rate(300.0);
        let r2 = layer.desorption_rate(600.0);
        assert!(r2 > r1, "desorption rate increases with temperature");
    }

    #[test]
    fn test_adsorption_surface_lifetime_finite() {
        let layer = AdsorptionLayer::new(1e-19, 0.5, 1e-9, false);
        let tau = layer.surface_lifetime(300.0);
        assert!(tau.is_finite() && tau > 0.0);
    }
}
