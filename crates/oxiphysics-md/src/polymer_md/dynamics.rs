#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dynamics models: RouseModel, ZimmModel, RouseChain, TubeModel, RingPolymer,
//! LangevinIntegrator, RouseModeTracker, ZimmHydrodynamics, ChainDiffusion,
//! ReptationDynamics.

use super::rouse_eigenvalue;
use std::f64::consts::PI;

/// Rouse model for unentangled polymer.
///
/// Computes normal modes, eigenvalues, mode amplitudes, and relaxation spectrum.
#[derive(Debug, Clone)]
pub struct RouseModel {
    /// Number of beads N.
    pub n: usize,
    /// Bond spring constant k.
    pub k_spring: f64,
    /// Friction coefficient zeta.
    pub zeta: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Mode amplitudes X_p.
    pub modes: Vec<f64>,
}

impl RouseModel {
    /// Create a new Rouse model.
    pub fn new(n: usize, k_spring: f64, zeta: f64, kt: f64) -> Self {
        let modes = vec![0.0f64; n];
        Self {
            n,
            k_spring,
            zeta,
            kt,
            modes,
        }
    }

    /// Rouse eigenvalue lambda_p = 2*k*(1-cos(p*pi/N)).
    pub fn eigenvalue(&self, p: usize) -> f64 {
        rouse_eigenvalue(p, self.n, self.k_spring)
    }

    /// Relaxation time tau_p = zeta / lambda_p.
    pub fn relaxation_time(&self, p: usize) -> f64 {
        let lam = self.eigenvalue(p);
        if lam < 1e-30 {
            return f64::INFINITY;
        }
        self.zeta / lam
    }

    /// Rouse relaxation time tau_1 (longest mode).
    pub fn tau_rouse(&self) -> f64 {
        self.relaxation_time(1)
    }

    /// Mean-square displacement at time t (diffusion regime).
    pub fn msd_center_of_mass(&self, t: f64) -> f64 {
        let d = self.kt / (self.n as f64 * self.zeta);
        6.0 * d * t
    }

    /// Mode amplitude <X_p^2> from equipartition.
    pub fn mode_amplitude_sq(&self, p: usize) -> f64 {
        let lam = self.eigenvalue(p);
        if lam < 1e-30 {
            return 0.0;
        }
        self.kt / lam
    }

    /// Viscosity contribution from Rouse model.
    pub fn rouse_viscosity(&self) -> f64 {
        let tau_r = self.tau_rouse();
        self.n as f64 * self.kt * tau_r
    }
}

/// Zimm model with hydrodynamic interactions.
///
/// Uses pre-averaging approximation with Oseen tensor.
#[derive(Debug, Clone)]
pub struct ZimmModel {
    /// Number of beads.
    pub n: usize,
    /// Bond length b.
    pub b: f64,
    /// Solvent viscosity eta_s.
    pub eta_s: f64,
    /// Thermal energy kT.
    pub kt: f64,
}

impl ZimmModel {
    /// Create a new Zimm model.
    pub fn new(n: usize, b: f64, eta_s: f64, kt: f64) -> Self {
        Self { n, b, eta_s, kt }
    }

    /// Zimm relaxation time tau_z = eta_s * R^3 / kT where R = b*N^(3/5).
    pub fn zimm_time(&self) -> f64 {
        let r = self.b * (self.n as f64).powf(0.5887); // Flory exponent 3/5 ~ 0.5887
        self.eta_s * r * r * r / self.kt
    }

    /// Intrinsic viscosity \[eta\] = eta_s * R^3 / N.
    pub fn intrinsic_viscosity(&self) -> f64 {
        let r = self.b * (self.n as f64).powf(0.5887);
        r * r * r / self.n as f64
    }

    /// Diffusion coefficient D = kT / (6*pi*eta_s*R).
    pub fn diffusion_coefficient(&self) -> f64 {
        let r = self.b * (self.n as f64).powf(0.5887);
        self.kt / (6.0 * PI * self.eta_s * r)
    }

    /// Translational friction coefficient.
    pub fn friction_coeff(&self) -> f64 {
        self.kt / self.diffusion_coefficient()
    }
}

/// Rouse chain dynamics simulator.
///
/// Propagates bead positions under harmonic connectivity and Langevin thermostat.
#[derive(Debug, Clone)]
pub struct RouseChain {
    /// Number of beads.
    pub n: usize,
    /// Spring constant k.
    pub k_spring: f64,
    /// Friction coefficient ζ per bead.
    pub zeta: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Time step dt.
    pub dt: f64,
    /// Bead positions.
    pub positions: Vec<f64>,
    /// Bead velocities.
    pub velocities: Vec<f64>,
}

impl RouseChain {
    /// Create a new Rouse chain.
    pub fn new(n: usize, k_spring: f64, zeta: f64, kt: f64, dt: f64) -> Self {
        let positions: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let velocities = vec![0.0f64; n];
        RouseChain {
            n,
            k_spring,
            zeta,
            kt,
            dt,
            positions,
            velocities,
        }
    }

    /// Harmonic force on bead i from its neighbors.
    pub fn spring_force(&self, i: usize) -> f64 {
        let mut f = 0.0f64;
        if i > 0 {
            f += self.k_spring * (self.positions[i - 1] - self.positions[i]);
        }
        if i + 1 < self.n {
            f += self.k_spring * (self.positions[i + 1] - self.positions[i]);
        }
        f
    }

    /// Rouse normal mode amplitude X_p.
    pub fn normal_mode(&self, p: usize) -> f64 {
        let n = self.n as f64;
        let norm = if p == 0 {
            1.0 / n.sqrt()
        } else {
            (2.0 / n).sqrt()
        };
        norm * self
            .positions
            .iter()
            .enumerate()
            .map(|(i, &x)| x * (PI * p as f64 * (i as f64 + 0.5) / n).cos())
            .sum::<f64>()
    }

    /// Relaxation time for mode p.
    pub fn tau_p(&self, p: usize) -> f64 {
        if p == 0 {
            return f64::INFINITY;
        }
        let lam = rouse_eigenvalue(p, self.n, self.k_spring);
        if lam < 1e-30 {
            f64::INFINITY
        } else {
            self.zeta / lam
        }
    }

    /// Mean-square displacement of center of mass at time t.
    pub fn com_msd(&self, t: f64) -> f64 {
        let d_cm = self.kt / (self.n as f64 * self.zeta);
        2.0 * d_cm * t
    }

    /// Internal MSD (Rouse regime): <|r_n(t) - r_n(0)|²> ~ t^(1/2) for t < τ_R.
    pub fn internal_msd_rouse(&self, t: f64) -> f64 {
        let b2 = 1.0; // assumed unit bond length
        let tau_1 = self.tau_p(1);
        let w = self.k_spring / self.zeta;
        let prefactor = b2 * (2.0 * w * t / PI).sqrt();
        prefactor.min(self.n as f64 * b2 * tau_1 / t)
    }

    /// Perform one Brownian dynamics (Euler-Maruyama) step without noise.
    pub fn deterministic_step(&mut self) {
        let forces: Vec<f64> = (0..self.n).map(|i| self.spring_force(i)).collect();
        for i in 0..self.n {
            self.positions[i] += forces[i] / self.zeta * self.dt;
        }
    }
}

/// Tube model for reptation dynamics (de Gennes / Doi-Edwards).
///
/// Describes the dynamics of a long, entangled polymer chain in a network
/// of obstacles (the "tube").
#[derive(Debug, Clone)]
pub struct TubeModel {
    /// Degree of polymerization N.
    pub n: usize,
    /// Entanglement strand length N_e.
    pub n_e: usize,
    /// Primitive path length L_tube.
    pub l_tube: f64,
    /// Tube diameter a.
    pub a: f64,
    /// Monomeric friction coefficient ζ.
    pub zeta: f64,
    /// Bond length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
}

impl TubeModel {
    /// Create a new tube model.
    pub fn new(n: usize, n_e: usize, a: f64, zeta: f64, b: f64, kt: f64) -> Self {
        let z = n as f64 / n_e as f64;
        let l_tube = z * a;
        TubeModel {
            n,
            n_e,
            l_tube,
            a,
            zeta,
            b,
            kt,
        }
    }

    /// Number of entanglements Z = N / N_e.
    pub fn z_entanglements(&self) -> f64 {
        self.n as f64 / self.n_e as f64
    }

    /// Rouse time of an entanglement strand: τ_e = ζ N_e² b² / (3π² kT).
    pub fn tau_e(&self) -> f64 {
        let ne = self.n_e as f64;
        self.zeta * ne * ne * self.b * self.b / (3.0 * PI * PI * self.kt)
    }

    /// Reptation (disengagement) time: τ_d = 3 Z³ τ_e.
    pub fn tau_d(&self) -> f64 {
        let z = self.z_entanglements();
        3.0 * z * z * z * self.tau_e()
    }

    /// Rouse time of the entire chain: τ_R = Z² τ_e.
    pub fn tau_rouse(&self) -> f64 {
        let z = self.z_entanglements();
        z * z * self.tau_e()
    }

    /// Reptation diffusion coefficient: D = kT / (3 ζ N Z).
    pub fn diffusion_coefficient(&self) -> f64 {
        let z = self.z_entanglements();
        let n = self.n as f64;
        self.kt / (3.0 * self.zeta * n * z)
    }

    /// Plateau modulus: G_N = ρ kT / N_e (in units where volume = 1).
    pub fn plateau_modulus(&self) -> f64 {
        self.kt / self.n_e as f64
    }

    /// Viscosity in reptation regime: η₀ ~ Z³ τ_e.
    pub fn viscosity(&self) -> f64 {
        let z = self.z_entanglements();
        self.plateau_modulus() * z * z * z * self.tau_e()
    }

    /// Mean-square displacement of chain end in 3 regimes.
    ///
    /// Returns (g1, g2, g3) for t < τ_e, τ_e < t < τ_R, t > τ_d.
    pub fn msd_regimes(&self, t: f64) -> f64 {
        let te = self.tau_e();
        let tr = self.tau_rouse();
        let td = self.tau_d();
        let b2 = self.b * self.b;
        if t < te {
            // Rouse within entanglement strand: ~ t^(1/2)
            b2 * (t / te).sqrt()
        } else if t < tr {
            // Constrained Rouse: ~ t^(1/4)
            self.a * self.b * (t / te).powf(0.25)
        } else if t < td {
            // Reptation along tube: ~ t^(1/2)
            self.a * self.a * (t / td).sqrt()
        } else {
            // Free diffusion: ~ t
            6.0 * self.diffusion_coefficient() * t
        }
    }
}

/// Ring polymer molecular dynamics (RPMD) for quantum effects.
///
/// Represents quantum particle as a ring of P classical beads connected
/// by harmonic springs with frequency ω_P = P kT / ℏ.
#[derive(Debug, Clone)]
pub struct RingPolymer {
    /// Number of beads P (Trotter number).
    pub n_beads: usize,
    /// Physical mass m.
    pub mass: f64,
    /// Temperature (in units kT).
    pub kt: f64,
    /// Planck constant ℏ (reduced).
    pub hbar: f64,
    /// Bead positions.
    pub positions: Vec<f64>,
    /// Bead momenta.
    pub momenta: Vec<f64>,
}

impl RingPolymer {
    /// Create a new ring polymer.
    pub fn new(n_beads: usize, mass: f64, kt: f64, hbar: f64) -> Self {
        let positions = vec![0.0f64; n_beads];
        let momenta = vec![0.0f64; n_beads];
        RingPolymer {
            n_beads,
            mass,
            kt,
            hbar,
            positions,
            momenta,
        }
    }

    /// Frequency of inter-bead springs: ω_P = P kT / ℏ.
    pub fn omega_p(&self) -> f64 {
        self.n_beads as f64 * self.kt / self.hbar
    }

    /// Spring constant for inter-bead harmonic potential: k_P = m ω_P².
    pub fn spring_constant(&self) -> f64 {
        self.mass * self.omega_p().powi(2)
    }

    /// Internal spring force on bead i (ring topology: i connects to i±1 mod P).
    pub fn spring_force(&self, i: usize) -> f64 {
        let p = self.n_beads;
        let next = (i + 1) % p;
        let prev = (i + p - 1) % p;
        let k = self.spring_constant();
        k * (self.positions[prev] + self.positions[next] - 2.0 * self.positions[i])
    }

    /// Centroid position (quantum mechanical average position).
    pub fn centroid(&self) -> f64 {
        self.positions.iter().sum::<f64>() / self.n_beads as f64
    }

    /// Centroid momentum.
    pub fn centroid_momentum(&self) -> f64 {
        self.momenta.iter().sum::<f64>() / self.n_beads as f64
    }

    /// Ring polymer internal energy (spring term only).
    pub fn internal_energy(&self) -> f64 {
        let k = self.spring_constant();
        let p = self.n_beads;
        (0..p)
            .map(|i| {
                let next = (i + 1) % p;
                let dr = self.positions[next] - self.positions[i];
                0.5 * k * dr * dr
            })
            .sum()
    }

    /// Estimator for kinetic energy using virial theorem.
    pub fn kinetic_energy_estimator(&self, external_force: &[f64]) -> f64 {
        let centroid = self.centroid();
        let p = self.n_beads as f64;
        let virial: f64 = self
            .positions
            .iter()
            .zip(external_force.iter())
            .map(|(&x, &f)| (x - centroid) * f)
            .sum();
        self.kt / 2.0 - virial / (2.0 * p)
    }

    /// Normal mode transformation (DFT-based for ring polymer).
    pub fn normal_mode_positions(&self) -> Vec<f64> {
        let p = self.n_beads;
        (0..p)
            .map(|k| {
                let re: f64 = self
                    .positions
                    .iter()
                    .enumerate()
                    .map(|(j, &x)| x * (2.0 * PI * j as f64 * k as f64 / p as f64).cos())
                    .sum();
                re / p as f64
            })
            .collect()
    }
}

/// Langevin dynamics integrator (BAOAB splitting scheme).
///
/// Implements the BAOAB Langevin integrator for canonical ensemble MD.
#[derive(Debug, Clone)]
pub struct LangevinIntegrator {
    /// Friction coefficient γ.
    pub gamma: f64,
    /// Temperature kT.
    pub kt: f64,
    /// Time step dt.
    pub dt: f64,
    /// Bead masses.
    pub masses: Vec<f64>,
}

impl LangevinIntegrator {
    /// Create a new Langevin integrator.
    pub fn new(gamma: f64, kt: f64, dt: f64, masses: Vec<f64>) -> Self {
        LangevinIntegrator {
            gamma,
            kt,
            dt,
            masses,
        }
    }

    /// BAOAB O-step: stochastic velocity update with friction and noise.
    ///
    /// v_new = v * exp(-γ dt) + sqrt((1 - exp(-2γdt)) kT/m) * ξ.
    pub fn o_step(&self, velocities: &mut Vec<[f64; 3]>, noise: &[[f64; 3]]) {
        let c1 = (-self.gamma * self.dt).exp();
        for (i, v) in velocities.iter_mut().enumerate() {
            let m = self.masses[i];
            let c2 = ((1.0 - c1 * c1) * self.kt / m).max(0.0).sqrt();
            let xi = if i < noise.len() { noise[i] } else { [0.0; 3] };
            v[0] = c1 * v[0] + c2 * xi[0];
            v[1] = c1 * v[1] + c2 * xi[1];
            v[2] = c1 * v[2] + c2 * xi[2];
        }
    }

    /// A-step: position update x += v * dt/2.
    pub fn a_step(&self, positions: &mut Vec<[f64; 3]>, velocities: &[[f64; 3]]) {
        let half_dt = 0.5 * self.dt;
        for (x, v) in positions.iter_mut().zip(velocities.iter()) {
            x[0] += v[0] * half_dt;
            x[1] += v[1] * half_dt;
            x[2] += v[2] * half_dt;
        }
    }

    /// B-step: velocity update v += F/m * dt/2.
    pub fn b_step(&self, velocities: &mut Vec<[f64; 3]>, forces: &[[f64; 3]]) {
        let half_dt = 0.5 * self.dt;
        for (i, (v, f)) in velocities.iter_mut().zip(forces.iter()).enumerate() {
            let m = self.masses[i];
            v[0] += f[0] / m * half_dt;
            v[1] += f[1] / m * half_dt;
            v[2] += f[2] / m * half_dt;
        }
    }

    /// Instantaneous kinetic energy.
    pub fn kinetic_energy(&self, velocities: &[[f64; 3]]) -> f64 {
        velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// Instantaneous temperature from kinetic energy.
    pub fn temperature(&self, velocities: &[[f64; 3]]) -> f64 {
        let n_dof = 3 * velocities.len();
        if n_dof == 0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy(velocities) / n_dof as f64
    }
}

/// Chain diffusion coefficient calculator for Rouse and Zimm models.
///
/// Provides predictions for the center-of-mass diffusion coefficient
/// in both unentangled (Rouse) and good-solvent (Zimm) limits.
#[derive(Debug, Clone)]
pub struct ChainDiffusion {
    /// Degree of polymerization N.
    pub n: usize,
    /// Bond length b.
    pub b: f64,
    /// Friction coefficient per bead ζ.
    pub zeta_bead: f64,
    /// Solvent viscosity η_s.
    pub eta_s: f64,
    /// Thermal energy kT.
    pub kt: f64,
}

impl ChainDiffusion {
    /// Create a chain diffusion calculator.
    pub fn new(n: usize, b: f64, zeta_bead: f64, eta_s: f64, kt: f64) -> Self {
        Self {
            n,
            b,
            zeta_bead,
            eta_s,
            kt,
        }
    }

    /// Rouse diffusion: D_Rouse = kT / (N ζ).
    pub fn rouse_diffusion(&self) -> f64 {
        self.kt / (self.n as f64 * self.zeta_bead)
    }

    /// Zimm diffusion: D_Zimm = kT / (6π η_s R_H) where R_H ~ b N^(3/5).
    pub fn zimm_diffusion(&self) -> f64 {
        let r_h = self.b * (self.n as f64).powf(0.588);
        self.kt / (6.0 * PI * self.eta_s * r_h)
    }

    /// Rouse relaxation time τ_R = ζ N² b² / (3π² kT).
    pub fn rouse_relaxation_time(&self) -> f64 {
        let n = self.n as f64;
        self.zeta_bead * n * n * self.b * self.b / (3.0 * PI * PI * self.kt)
    }

    /// Zimm relaxation time τ_Z = η_s b³ N^(3ν) / kT, ν = 3/5.
    pub fn zimm_relaxation_time(&self) -> f64 {
        let n = self.n as f64;
        self.eta_s * self.b.powi(3) * n.powf(1.764) / self.kt
    }

    /// Mean-square displacement of center of mass at time t.
    ///
    /// In Rouse regime: <Δr²_cm(t)> = 6 D_Rouse t.
    pub fn cm_msd_rouse(&self, t: f64) -> f64 {
        6.0 * self.rouse_diffusion() * t
    }

    /// MSD in Zimm regime.
    pub fn cm_msd_zimm(&self, t: f64) -> f64 {
        6.0 * self.zimm_diffusion() * t
    }

    /// Crossover time between Rouse and diffusion regimes.
    pub fn rouse_crossover_time(&self) -> f64 {
        self.rouse_relaxation_time()
    }

    /// Hydrodynamic radius R_H.
    pub fn hydrodynamic_radius(&self) -> f64 {
        self.b * (self.n as f64).powf(0.588)
    }
}

/// Extended reptation / tube model with contour length fluctuations (CLF)
/// and constraint release (CR) corrections.
///
/// Based on Doi-Edwards-Milner-McLeish (DEMM) theory.
///
/// References:
/// - Doi, M. & Edwards, S.F. (1986). *The Theory of Polymer Dynamics*.
/// - Milner, S.T. & McLeish, T.C.B. (1998). *Phys. Rev. Lett.* 81, 725.
#[derive(Debug, Clone)]
pub struct ReptationDynamics {
    /// Degree of polymerization N.
    pub n: usize,
    /// Entanglement strand length N_e.
    pub n_e: usize,
    /// Monomeric friction coefficient ζ.
    pub zeta: f64,
    /// Bond length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Constraint release (CR) factor β_CR (0 = no CR, 1 = full CR).
    pub cr_factor: f64,
}

impl ReptationDynamics {
    /// Create a new reptation dynamics model.
    pub fn new(n: usize, n_e: usize, zeta: f64, b: f64, kt: f64, cr_factor: f64) -> Self {
        Self {
            n,
            n_e,
            zeta,
            b,
            kt,
            cr_factor,
        }
    }

    /// Number of entanglements Z = N / N_e.
    pub fn z(&self) -> f64 {
        self.n as f64 / self.n_e as f64
    }

    /// Rouse time of an entanglement segment: τ_e.
    pub fn tau_e(&self) -> f64 {
        let ne = self.n_e as f64;
        self.zeta * ne * ne * self.b * self.b / (3.0 * PI * PI * self.kt)
    }

    /// Reptation (bare) time τ_d,0 = 3 Z³ τ_e.
    pub fn tau_d0(&self) -> f64 {
        3.0 * self.z().powi(3) * self.tau_e()
    }

    /// CLF correction factor: reduces τ_d below bare value.
    ///
    /// f_CLF ≈ (1 - C/Z^(1/2))² for large Z, C ≈ 1.47.
    pub fn clf_correction(&self) -> f64 {
        let z = self.z();
        if z < 1.0 {
            return 1.0;
        }
        let c = 1.47;
        let f = 1.0 - c / z.sqrt();
        if f < 0.0 { 0.01 } else { f * f }
    }

    /// Corrected reptation time with CLF: τ_d = τ_d,0 * f_CLF.
    pub fn tau_d_clf(&self) -> f64 {
        self.tau_d0() * self.clf_correction()
    }

    /// Constraint release correction: τ_d,CR = τ_d * (1 + β_CR / Z).
    pub fn tau_d_cr(&self) -> f64 {
        let z = self.z();
        self.tau_d_clf() * (1.0 + self.cr_factor / z)
    }

    /// Diffusion coefficient D = kT / (3 ζ N Z).
    pub fn diffusion_coefficient(&self) -> f64 {
        let z = self.z();
        let n = self.n as f64;
        self.kt / (3.0 * self.zeta * n * z)
    }

    /// Zero-shear viscosity η_0 ~ Z³.
    ///
    /// η₀ = G_N^0 * τ_d / (1 + β_CR/Z)
    pub fn viscosity_zero_shear(&self) -> f64 {
        let gn0 = self.kt / self.n_e as f64;
        gn0 * self.tau_d_cr()
    }

    /// Steady-state compliance J_e^0 = 6/(5 G_N^0 Z).
    pub fn steady_state_compliance(&self) -> f64 {
        let z = self.z();
        let gn0 = self.kt / self.n_e as f64;
        6.0 / (5.0 * gn0 * z)
    }

    /// Tube survival probability P(t) in reptation regime (pure reptation).
    ///
    /// P(t) = (8/π²) Σ_{p odd} (1/p²) exp(-p² t / τ_d)
    pub fn tube_survival_prob(&self, t: f64, n_terms: usize) -> f64 {
        let tau = self.tau_d_cr().max(1e-30);
        let mut sum = 0.0f64;
        for k in 0..n_terms {
            let p = (2 * k + 1) as f64;
            sum += (-p * p * t / tau).exp() / (p * p);
        }
        8.0 / (PI * PI) * sum
    }

    /// Relaxation modulus G(t) = G_N^0 * P(t).
    pub fn relaxation_modulus(&self, t: f64, n_terms: usize) -> f64 {
        let gn0 = self.kt / self.n_e as f64;
        gn0 * self.tube_survival_prob(t, n_terms)
    }
}

/// Rouse model with explicit normal mode dynamics.
///
/// Tracks all N normal modes X_p and their time evolution under
/// the Langevin thermostat. Provides autocorrelation functions.
#[derive(Debug, Clone)]
pub struct RouseModeTracker {
    /// Number of beads N.
    pub n: usize,
    /// Spring constant k.
    pub k: f64,
    /// Friction coefficient ζ.
    pub zeta: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Time step dt.
    pub dt: f64,
    /// Mode amplitudes X_p (n values).
    pub amplitudes: Vec<f64>,
    /// Mode amplitude history for correlation functions.
    pub history: Vec<Vec<f64>>,
    /// Max history length.
    pub max_history: usize,
}

impl RouseModeTracker {
    /// Create a Rouse mode tracker.
    pub fn new(n: usize, k: f64, zeta: f64, kt: f64, dt: f64, max_history: usize) -> Self {
        Self {
            n,
            k,
            zeta,
            kt,
            dt,
            amplitudes: vec![0.0f64; n],
            history: Vec::new(),
            max_history,
        }
    }

    /// Rouse eigenvalue λ_p = 2k(1 - cos(pπ/N)).
    pub fn lambda_p(&self, p: usize) -> f64 {
        2.0 * self.k * (1.0 - (p as f64 * PI / self.n as f64).cos())
    }

    /// Relaxation time τ_p = ζ / λ_p.
    pub fn tau_p(&self, p: usize) -> f64 {
        let lam = self.lambda_p(p);
        if lam < 1e-30 {
            f64::INFINITY
        } else {
            self.zeta / lam
        }
    }

    /// Equipartition amplitude <X_p²> = kT / λ_p.
    pub fn equilibrium_amplitude_sq(&self, p: usize) -> f64 {
        let lam = self.lambda_p(p);
        if lam < 1e-30 { 0.0 } else { self.kt / lam }
    }

    /// Deterministic relaxation of mode p by one step.
    ///
    /// dX_p/dt = -λ_p X_p / ζ
    pub fn relax_mode(&mut self, p: usize) {
        if p >= self.n {
            return;
        }
        let tau = self.tau_p(p);
        self.amplitudes[p] *= (-self.dt / tau.max(1e-30)).exp();
    }

    /// Relax all modes by one step.
    pub fn relax_all_modes(&mut self) {
        for p in 0..self.n {
            self.relax_mode(p);
        }
        if self.history.len() >= self.max_history {
            self.history.remove(0);
        }
        self.history.push(self.amplitudes.clone());
    }

    /// Autocorrelation <X_p(t) X_p(0)> at lag τ steps.
    pub fn mode_autocorrelation(&self, p: usize, lag: usize) -> f64 {
        let n = self.history.len();
        if n == 0 || p >= self.n || lag >= n {
            return 0.0;
        }
        let count = n - lag;
        let sum: f64 = (0..count)
            .map(|t| self.history[t][p] * self.history[t + lag][p])
            .sum();
        sum / count as f64
    }

    /// Theoretical autocorrelation: C_p(t) = <X_p²> * exp(-t/τ_p).
    pub fn theory_autocorr(&self, p: usize, t: f64) -> f64 {
        let amps2 = self.equilibrium_amplitude_sq(p);
        let tau = self.tau_p(p);
        amps2 * (-t / tau.max(1e-30)).exp()
    }

    /// Project bead positions onto normal mode p.
    ///
    /// X_p = sqrt(2/N) Σ_j r_j cos(π p (j+0.5) / N)
    pub fn project_positions(&self, positions: &[f64], p: usize) -> f64 {
        let n = self.n as f64;
        let norm = if p == 0 {
            (1.0 / n).sqrt()
        } else {
            (2.0 / n).sqrt()
        };
        norm * positions
            .iter()
            .enumerate()
            .map(|(j, &r)| r * (PI * p as f64 * (j as f64 + 0.5) / n).cos())
            .sum::<f64>()
    }
}

/// Zimm model with hydrodynamic interactions via Oseen tensor.
///
/// The pre-averaged Oseen tensor Tᵢⱼ = (6π η_s |r_ij|)⁻¹ couples the
/// friction of different beads through the solvent.
///
/// References:
/// - Zimm, B.H. (1956). *J. Chem. Phys.* 24, 269.
/// - Doi, M. & Edwards, S.F. (1986). *Theory of Polymer Dynamics*.
#[derive(Debug, Clone)]
pub struct ZimmHydrodynamics {
    /// Number of beads N.
    pub n: usize,
    /// Bond length b.
    pub b: f64,
    /// Solvent viscosity η_s.
    pub eta_s: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Bead positions (1D: x-coordinates only for simplicity).
    pub positions: Vec<f64>,
}

impl ZimmHydrodynamics {
    /// Create a new Zimm hydrodynamics model.
    pub fn new(n: usize, b: f64, eta_s: f64, kt: f64) -> Self {
        let positions: Vec<f64> = (0..n).map(|i| i as f64 * b).collect();
        Self {
            n,
            b,
            eta_s,
            kt,
            positions,
        }
    }

    /// Pre-averaged Oseen tensor element T_ij.
    ///
    /// T_ij = 1 / (6π η_s <|r_ij|>) where <|r_ij|> = b sqrt(|i-j|) (Rouse ideal).
    pub fn oseen_element(&self, i: usize, j: usize) -> f64 {
        if i == j {
            return 0.0;
        }
        let delta_n = (i as isize - j as isize).unsigned_abs();
        let r_avg = self.b * (delta_n as f64).sqrt().max(1e-30);
        1.0 / (6.0 * PI * self.eta_s * r_avg)
    }

    /// Effective friction tensor element Ξ_ij = ζ δ_ij + Σ_k T_ik.
    ///
    /// Using the Kirkwood-Riseman approximation: compute diagonal element.
    pub fn friction_diagonal(&self, i: usize) -> f64 {
        let oseen_sum: f64 = (0..self.n)
            .filter(|&j| j != i)
            .map(|j| self.oseen_element(i, j))
            .sum();
        let zeta_0 = 1.0 / (6.0 * PI * self.eta_s * self.b);
        1.0 / (1.0 / zeta_0 + oseen_sum)
    }

    /// Zimm diffusion coefficient D = kT Σ_i (1/Ξ_ii) / N².
    pub fn diffusion_coefficient(&self) -> f64 {
        let sum: f64 = (0..self.n)
            .map(|i| 1.0 / self.friction_diagonal(i).max(1e-30))
            .sum();
        self.kt * sum / (self.n as f64 * self.n as f64)
    }

    /// Zimm relaxation time of mode p.
    ///
    /// τ_p^Zimm = η_s b³ N^(3ν) / (kT * p^(3ν)), ν = 3/5.
    pub fn zimm_tau_p(&self, p: usize) -> f64 {
        if p == 0 {
            return f64::INFINITY;
        }
        let n = self.n as f64;
        let p_f = p as f64;
        self.eta_s * self.b.powi(3) * n.powf(1.764) / (self.kt * p_f.powf(1.764))
    }

    /// Longest Zimm relaxation time τ_1.
    pub fn zimm_time(&self) -> f64 {
        self.zimm_tau_p(1)
    }

    /// Intrinsic viscosity \[η\] = kT τ_1 / (η_s R³).
    pub fn intrinsic_viscosity(&self) -> f64 {
        let r = self.b * (self.n as f64).powf(0.588);
        self.kt * self.zimm_time() / (self.eta_s * r.powi(3))
    }
}
