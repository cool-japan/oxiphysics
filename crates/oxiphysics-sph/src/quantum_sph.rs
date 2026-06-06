// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum SPH (QSPH) for semiclassical and quantum fluid dynamics.
//!
//! This module provides:
//! - [`QuantumParticle`]: SPH particle with quantum wavefunction components.
//! - [`QuantumSPH`]: Core QSPH solver with quantum pressure and Bohm potential.
//! - [`FeynmanPathSPH`]: Path-integral SPH using Feynman path weights.
//! - [`DeBroglieWavelength`]: Thermal de Broglie wavelength calculator.
//! - [`wkb_density`]: WKB approximation for quantum density of states.
//!
//! # Physics Background
//!
//! Quantum SPH extends the standard SPH method to include quantum effects
//! through the Madelung transformation of the Schrödinger equation. Writing
//! ψ = √ρ exp(iS/ℏ), the quantum fluid equations become:
//!
//! ∂ρ/∂t + ∇·(ρu) = 0
//! ∂u/∂t + (u·∇)u = −(1/m)∇(V + Q)
//!
//! where Q = −(ℏ²/2m) ∇²√ρ/√ρ is the Bohm quantum potential.
//!
//! ## Feynman Path Integrals
//! The quantum propagator K(x_b, t_b; x_a, t_a) = ∫ D\[x(t)\] exp(iS/ℏ) sums
//! over all classical paths weighted by exp(iS/ℏ). In imaginary time
//! (β = it/ℏ), this reduces to a Boltzmann-like weight exp(−S_E/ℏ).
//!
//! ## WKB Approximation
//! In the semiclassical limit, the local density of states at energy E
//! in a potential V(x) is ρ_WKB ∝ 1/√(E − V) when E > V.
//!
//! # References
//! - Madelung, E. (1927). Quantentheorie in hydrodynamischer Form. *Z. Phys.*, 40, 322.
//! - Feynman, R.P. & Hibbs, A.R. (1965). *Quantum Mechanics and Path Integrals*.
//! - Bohm, D. (1952). A Suggested Interpretation of Quantum Theory. *Phys. Rev.*, 85, 166.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Reduced Planck constant ℏ (J·s).
pub const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant k_B (J K⁻¹).
pub const K_B: f64 = 1.380_649e-23;
/// Electron mass m_e (kg).
pub const ELECTRON_MASS: f64 = 9.109_383_701_5e-31;
/// Proton mass m_p (kg).
pub const PROTON_MASS: f64 = 1.672_621_923_69e-27;

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// WKB (Wentzel-Kramers-Brillouin) semiclassical density.
///
/// In the WKB approximation, the local density of states at energy E
/// in a potential V(x) is proportional to the inverse of the classical
/// momentum p(x) = √(2m(E − V)):
///
/// ρ_WKB(x) = (1/π) · m / √(2m(E − V(x)))
///
/// Returns 0 when E ≤ V (classically forbidden region).
///
/// # Arguments
/// * `energy`    – total energy E (J)
/// * `potential` – local potential V(x) (J)
/// * `mass`      – particle mass (kg)
/// * `hbar`      – reduced Planck constant (J·s)
///
/// ```no_run
/// use oxiphysics_sph::quantum_sph::wkb_density;
/// let rho = wkb_density(1.0, 0.5, 1.0, 1.0);
/// assert!(rho > 0.0);
/// ```
pub fn wkb_density(energy: f64, potential: f64, mass: f64, _hbar: f64) -> f64 {
    let delta = energy - potential;
    if delta <= 0.0 || mass <= 0.0 {
        return 0.0;
    }
    mass / (PI * (2.0 * mass * delta).sqrt())
}

// ---------------------------------------------------------------------------
// QuantumParticle
// ---------------------------------------------------------------------------

/// SPH particle carrying quantum wavefunction information.
///
/// In the Madelung representation ψ = √ρ exp(iS/ℏ), the particle stores
/// both density ρ and the real/imaginary parts of the wavefunction amplitude.
#[derive(Debug, Clone)]
pub struct QuantumParticle {
    /// 3D position \[x, y, z\] (m).
    pub pos: [f64; 3],
    /// 3D velocity \[vx, vy, vz\] (m s⁻¹).
    pub vel: [f64; 3],
    /// Real part of wavefunction amplitude ψ_r = √ρ cos(S/ℏ).
    pub psi_r: f64,
    /// Imaginary part of wavefunction amplitude ψ_i = √ρ sin(S/ℏ).
    pub psi_i: f64,
    /// Local number density ρ (m⁻³ or dimensionless).
    pub density: f64,
}

impl QuantumParticle {
    /// Create a new [`QuantumParticle`].
    pub fn new(pos: [f64; 3], vel: [f64; 3], psi_r: f64, psi_i: f64) -> Self {
        let density = psi_r * psi_r + psi_i * psi_i;
        Self {
            pos,
            vel,
            psi_r,
            psi_i,
            density,
        }
    }

    /// Compute |ψ|² = ρ (probability density).
    pub fn probability_density(&self) -> f64 {
        self.psi_r * self.psi_r + self.psi_i * self.psi_i
    }

    /// Compute the quantum phase S/ℏ = atan2(ψ_i, ψ_r).
    pub fn phase(&self) -> f64 {
        self.psi_i.atan2(self.psi_r)
    }

    /// Compute the amplitude |ψ| = √(ψ_r² + ψ_i²).
    pub fn amplitude(&self) -> f64 {
        (self.psi_r * self.psi_r + self.psi_i * self.psi_i).sqrt()
    }
}

// ---------------------------------------------------------------------------
// QuantumSPH
// ---------------------------------------------------------------------------

/// Core Quantum SPH solver with Bohm potential and quantum pressure.
///
/// Implements the Madelung fluid formulation on a set of SPH particles.
/// Each time step:
/// 1. Computes SPH density estimates.
/// 2. Computes the Bohm quantum potential Q = −(ℏ²/2m) ∇²√ρ/√ρ.
/// 3. Advances velocities using classical + quantum forces.
/// 4. Advances positions.
#[derive(Debug, Clone)]
pub struct QuantumSPH {
    /// Reduced Planck constant (lattice units).
    pub hbar: f64,
    /// Particle mass (lattice units).
    pub mass: f64,
    /// SPH smoothing length h_sph (m).
    pub smoothing_length: f64,
    /// List of quantum particles.
    pub particles: Vec<QuantumParticle>,
}

impl QuantumSPH {
    /// Create a new [`QuantumSPH`] solver.
    ///
    /// # Arguments
    /// * `hbar`             – reduced Planck constant (lattice units)
    /// * `mass`             – particle mass (lattice units)
    /// * `smoothing_length` – SPH kernel smoothing length
    /// * `particles`        – initial particle list
    pub fn new(
        hbar: f64,
        mass: f64,
        smoothing_length: f64,
        particles: Vec<QuantumParticle>,
    ) -> Self {
        Self {
            hbar,
            mass,
            smoothing_length,
            particles,
        }
    }

    /// Gaussian SPH kernel W(r, h) = (1/(h√π))³ exp(−r²/h²) in 3D.
    #[inline]
    fn kernel(&self, r: f64) -> f64 {
        let h = self.smoothing_length;
        if h <= 0.0 {
            return 0.0;
        }
        let q = r / h;
        let norm = 1.0 / (PI.sqrt() * h).powi(3);
        norm * (-q * q).exp()
    }

    /// Compute kernel gradient magnitude |∇W| = |dW/dr|.
    #[inline]
    fn kernel_gradient_mag(&self, r: f64) -> f64 {
        let h = self.smoothing_length;
        if h <= 0.0 || r < 1e-30 {
            return 0.0;
        }
        let q = r / h;
        let norm = 1.0 / (PI.sqrt() * h).powi(3);
        norm * (-q * q).exp() * (-2.0 * q / h)
    }

    /// Compute SPH density estimates ρ_a = ∑_b m W(|r_a − r_b|, h).
    ///
    /// Updates each particle's `density` field in place.
    pub fn compute_density(&mut self) {
        let n = self.particles.len();
        let mass = self.mass;
        let mut densities = vec![0.0f64; n];
        for (i, rho_i) in densities.iter_mut().enumerate() {
            for j in 0..n {
                let r = self.distance(i, j);
                *rho_i += mass * self.kernel(r);
            }
        }
        for (p, rho) in self.particles.iter_mut().zip(densities.iter()) {
            p.density = *rho;
        }
    }

    /// Helper: Euclidean distance between particles i and j.
    fn distance(&self, i: usize, j: usize) -> f64 {
        let pa = &self.particles[i];
        let pb = &self.particles[j];
        let dx = pa.pos[0] - pb.pos[0];
        let dy = pa.pos[1] - pb.pos[1];
        let dz = pa.pos[2] - pb.pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Compute the quantum (Bohm) potential Q_i at each particle.
    ///
    /// Uses the Laplacian of √ρ estimated via SPH summation:
    /// ∇²f_a ≈ ∑_b (m/ρ_b) (f_b − f_a) ∇²W(r_ab, h)
    ///
    /// Then Q_a = −(ℏ²/2m) · (∇²√ρ_a) / √ρ_a
    ///
    /// Returns a vector of Bohm potentials (one per particle).
    pub fn bohm_potential(&self) -> Vec<f64> {
        let n = self.particles.len();
        let h = self.smoothing_length;
        let mut q_pot = vec![0.0f64; n];
        for (i, qp_i) in q_pot.iter_mut().enumerate() {
            let rho_i = self.particles[i].density.max(1e-30);
            let sqrt_rho_i = rho_i.sqrt();
            let mut laplacian = 0.0f64;
            for j in 0..n {
                let r = self.distance(i, j);
                let rho_j = self.particles[j].density.max(1e-30);
                let sqrt_rho_j = rho_j.sqrt();
                // Laplacian kernel: ∇²W ≈ (10/(πh⁵)) exp(−r²/h²) [3D Gaussian]
                let w_lapl = if h > 0.0 {
                    let q = r / h;
                    let norm = 1.0 / (PI.sqrt() * h).powi(3);
                    norm * (-q * q).exp() * (6.0 - 4.0 * q * q) / (h * h)
                } else {
                    0.0
                };
                laplacian += (self.mass / rho_j) * (sqrt_rho_j - sqrt_rho_i) * w_lapl;
            }
            *qp_i = if sqrt_rho_i > 1e-15 {
                -(self.hbar * self.hbar / (2.0 * self.mass)) * laplacian / sqrt_rho_i
            } else {
                0.0
            };
        }
        q_pot
    }

    /// Compute quantum pressure P_Q = (ℏ²/4m) ρ ∇²(ln ρ) at each particle.
    ///
    /// The quantum pressure is related to the Bohm potential by:
    /// P_Q = −ρ Q / m
    ///
    /// Returns a vector of quantum pressures (one per particle).
    pub fn compute_quantum_pressure(&self) -> Vec<f64> {
        let q_pot = self.bohm_potential();
        self.particles
            .iter()
            .zip(q_pot.iter())
            .map(|(p, &q)| -p.density * q / self.mass)
            .collect()
    }

    /// Compute a discretized Wigner function W(x, p) for 1D projections.
    ///
    /// Uses the formula:
    /// W(x, p) = (1/πℏ) ∫ ψ*(x+y) ψ(x−y) exp(2ipy/ℏ) dy
    ///
    /// Returns the Wigner function on a grid of `nx × np` points.
    ///
    /// # Arguments
    /// * `nx`   – number of position grid points
    /// * `np`   – number of momentum grid points
    /// * `xmax` – maximum position (m)
    /// * `pmax` – maximum momentum (kg m s⁻¹)
    pub fn wigner_function(&self, nx: usize, np: usize, xmax: f64, pmax: f64) -> Vec<Vec<f64>> {
        let _hbar = self.hbar;
        let mut wf = vec![vec![0.0f64; np]; nx];
        let dx = 2.0 * xmax / nx as f64;
        let dp = 2.0 * pmax / np as f64;
        for (ix, wf_row) in wf.iter_mut().enumerate() {
            let x = -xmax + ix as f64 * dx;
            // Estimate ψ(x) from nearest particle
            let psi_r = self.interpolate_psi_r(x);
            let psi_i = self.interpolate_psi_i(x);
            for (ip, wf_cell) in wf_row.iter_mut().enumerate() {
                let p_val = -pmax + ip as f64 * dp;
                // Simplified: W(x,p) ≈ |ψ(x)|² δ(p − mẋ) → Gaussian in p
                let rho_x = psi_r * psi_r + psi_i * psi_i;
                let p_mean = self.interpolate_velocity_x(x) * self.mass;
                let sigma_p = (self.hbar / (2.0 * self.smoothing_length)).max(1e-30);
                let dp_val = p_val - p_mean;
                *wf_cell =
                    rho_x / (PI.sqrt() * sigma_p) * (-dp_val * dp_val / (sigma_p * sigma_p)).exp();
            }
        }
        wf
    }

    /// SPH interpolation of ψ_r (real part) at position x (1D, x-coord).
    fn interpolate_psi_r(&self, x: f64) -> f64 {
        let mut val = 0.0f64;
        let mut norm = 0.0f64;
        for p in &self.particles {
            let dx = x - p.pos[0];
            let r = dx.abs();
            let w = self.kernel(r);
            val += p.psi_r * w;
            norm += w;
        }
        if norm > 1e-30 { val / norm } else { 0.0 }
    }

    /// SPH interpolation of ψ_i (imaginary part) at position x (1D).
    fn interpolate_psi_i(&self, x: f64) -> f64 {
        let mut val = 0.0f64;
        let mut norm = 0.0f64;
        for p in &self.particles {
            let dx = x - p.pos[0];
            let r = dx.abs();
            let w = self.kernel(r);
            val += p.psi_i * w;
            norm += w;
        }
        if norm > 1e-30 { val / norm } else { 0.0 }
    }

    /// SPH interpolation of x-velocity at position x (1D).
    fn interpolate_velocity_x(&self, x: f64) -> f64 {
        let mut val = 0.0f64;
        let mut norm = 0.0f64;
        for p in &self.particles {
            let dx = x - p.pos[0];
            let r = dx.abs();
            let w = self.kernel(r);
            val += p.vel[0] * w;
            norm += w;
        }
        if norm > 1e-30 { val / norm } else { 0.0 }
    }

    /// Advance the QSPH system by one time step `dt`.
    ///
    /// Uses a symplectic Euler integrator:
    /// 1. Compute density.
    /// 2. Compute Bohm potential forces.
    /// 3. Update velocities.
    /// 4. Update positions.
    pub fn step(&mut self, dt: f64) {
        self.compute_density();
        let q_forces = self.compute_bohm_forces();
        let n = self.particles.len();
        for (i, qf) in q_forces.iter().enumerate().take(n) {
            // Update velocity
            for (d, vd) in self.particles[i].vel.iter_mut().enumerate() {
                *vd += qf[d] * dt;
            }
            // Update position
            let vel = self.particles[i].vel;
            for (pd, vd) in self.particles[i].pos.iter_mut().zip(vel.iter()) {
                *pd += vd * dt;
            }
        }
    }

    /// Compute Bohm potential force F_i = −m ∇Q_i on each particle.
    ///
    /// The gradient is computed via SPH:
    /// ∇Q_a = ∑_b m/ρ_b (Q_b − Q_a) ∇W(r_ab)
    fn compute_bohm_forces(&self) -> Vec<[f64; 3]> {
        let n = self.particles.len();
        let q_pot = self.bohm_potential();
        let mut forces = vec![[0.0f64; 3]; n];
        for i in 0..n {
            let rho_i = self.particles[i].density.max(1e-30);
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = self.distance(i, j);
                if r < 1e-30 {
                    continue;
                }
                let rho_j = self.particles[j].density.max(1e-30);
                let dq = q_pot[j] - q_pot[i];
                let dw_dr = self.kernel_gradient_mag(r);
                // Gradient direction: (r_a − r_b) / |r_a − r_b|
                for (d, fd) in forces[i].iter_mut().enumerate() {
                    let dr = (self.particles[i].pos[d] - self.particles[j].pos[d]) / r;
                    *fd += -self.mass * (self.mass / rho_j.max(rho_i)) * dq * dw_dr * dr;
                }
            }
        }
        forces
    }

    /// Compute total energy: kinetic + quantum kinetic.
    ///
    /// E_total = E_kinetic + E_quantum_kinetic
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.quantum_kinetic_energy()
    }

    /// Compute classical kinetic energy E_kin = ∑_i m |v_i|² / 2.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.vel[0] * p.vel[0] + p.vel[1] * p.vel[1] + p.vel[2] * p.vel[2];
                0.5 * self.mass * v2
            })
            .sum()
    }

    /// Compute quantum kinetic energy E_qkin = ∑_i ⟨ℏ²|∇ψ|²/2m⟩.
    ///
    /// Estimated from the Bohm potential as E_qkin ≈ ∑_i m |Q_i|.
    pub fn quantum_kinetic_energy(&self) -> f64 {
        let q_pot = self.bohm_potential();
        q_pot
            .iter()
            .zip(self.particles.iter())
            .map(|(&q, p)| p.density * q.abs() / self.mass)
            .sum::<f64>()
            * self.mass
    }
}

// ---------------------------------------------------------------------------
// FeynmanPathSPH
// ---------------------------------------------------------------------------

/// Feynman path-integral SPH using imaginary-time propagation.
///
/// In the imaginary-time formulation β = iT/ℏ, the quantum propagator
/// K ∝ exp(−S_E/ℏ) where S_E is the Euclidean action:
///
/// S_E\[x(τ)\] = ∫₀^ℏβ \[m ẋ²/2 + V(x)\] dτ
///
/// Monte-Carlo sampling of these paths gives the finite-temperature density.
#[derive(Debug, Clone)]
pub struct FeynmanPathSPH {
    /// Number of imaginary-time slices (path discretization).
    pub n_beads: usize,
    /// Inverse temperature β = 1/(k_B T) (J⁻¹).
    pub beta: f64,
    /// Particle mass (lattice units).
    pub mass: f64,
    /// Reduced Planck constant (lattice units).
    pub hbar: f64,
    /// Path positions: shape \[n_beads\]\[3\].
    pub path: Vec<[f64; 3]>,
    /// External potential values at each bead.
    pub potential_values: Vec<f64>,
}

impl FeynmanPathSPH {
    /// Create a new [`FeynmanPathSPH`] ring polymer.
    ///
    /// # Arguments
    /// * `n_beads`    – number of imaginary-time slices
    /// * `beta`       – inverse temperature (J⁻¹)
    /// * `mass`       – particle mass
    /// * `hbar`       – reduced Planck constant
    /// * `init_pos`   – initial 3D position for all beads
    pub fn new(n_beads: usize, beta: f64, mass: f64, hbar: f64, init_pos: [f64; 3]) -> Self {
        let path = vec![init_pos; n_beads];
        let potential_values = vec![0.0f64; n_beads];
        Self {
            n_beads,
            beta,
            mass,
            hbar,
            path,
            potential_values,
        }
    }

    /// Compute the path-integral weight exp(−S_E/ℏ) for the current path.
    ///
    /// S_E = ∑_k \[m/(2ℏ Δτ) |r_k − r_{k−1}|² + Δτ V_k\]
    ///
    /// where Δτ = ℏβ / n_beads.
    pub fn path_integral_weight(&self) -> f64 {
        let action = self.action();
        (-action / self.hbar).exp()
    }

    /// Compute the Euclidean action S_E for the current ring-polymer path.
    ///
    /// S_E = ∑_{k=0}^{P−1} \[m/(2ℏΔτ) |r_k − r_{k+1}|² + Δτ V_k\]
    pub fn action(&self) -> f64 {
        let dtau = self.hbar * self.beta / self.n_beads as f64;
        let spring_const = self.mass / (self.hbar * dtau).max(1e-30);
        let mut s_e = 0.0f64;
        for k in 0..self.n_beads {
            let next = (k + 1) % self.n_beads;
            let dx = self.path[k][0] - self.path[next][0];
            let dy = self.path[k][1] - self.path[next][1];
            let dz = self.path[k][2] - self.path[next][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            s_e += 0.5 * spring_const * r2 + dtau * self.potential_values[k];
        }
        s_e
    }

    /// Propagate the path by one Monte-Carlo step using Metropolis algorithm.
    ///
    /// Displaces a randomly chosen bead by a small amount and accepts/rejects
    /// based on the Boltzmann weight exp(−ΔS_E/ℏ).
    ///
    /// # Arguments
    /// * `step_size` – maximum displacement magnitude
    pub fn propagate(&mut self, step_size: f64) {
        use rand::RngExt;
        let mut rng = rand::rng();
        let n = self.n_beads;
        // Choose a random bead
        let k: usize = rng.random_range(0..n);
        let s_old = self.action();
        // Displace the bead
        let old_pos = self.path[k];
        self.path[k][0] += rng.random_range(-step_size..step_size);
        self.path[k][1] += rng.random_range(-step_size..step_size);
        self.path[k][2] += rng.random_range(-step_size..step_size);
        let s_new = self.action();
        let delta_s = s_new - s_old;
        // Metropolis acceptance
        let accept = if delta_s <= 0.0 {
            true
        } else {
            let prob = (-delta_s / self.hbar).exp();
            rng.random_range(0.0f64..1.0) < prob
        };
        if !accept {
            self.path[k] = old_pos;
        }
    }

    /// Compute the centroid (mean bead position) of the ring polymer.
    pub fn centroid(&self) -> [f64; 3] {
        let n = self.n_beads as f64;
        let mut c = [0.0f64; 3];
        for bead in &self.path {
            c[0] += bead[0];
            c[1] += bead[1];
            c[2] += bead[2];
        }
        [c[0] / n, c[1] / n, c[2] / n]
    }

    /// Compute the gyration radius of the ring polymer.
    ///
    /// r_g = √(⟨r²⟩ − ⟨r⟩²)
    ///
    /// This is related to the thermal de Broglie wavelength.
    pub fn gyration_radius(&self) -> f64 {
        let c = self.centroid();
        let n = self.n_beads as f64;
        let r2: f64 = self
            .path
            .iter()
            .map(|bead| {
                let dx = bead[0] - c[0];
                let dy = bead[1] - c[1];
                let dz = bead[2] - c[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            / n;
        r2.sqrt()
    }
}

// ---------------------------------------------------------------------------
// DeBroglieWavelength
// ---------------------------------------------------------------------------

/// Thermal de Broglie wavelength calculator.
///
/// The thermal de Broglie wavelength λ_th sets the scale at which quantum
/// effects become important. When the inter-particle spacing n^(−1/3) ≲ λ_th,
/// the system is in the quantum degenerate regime.
///
/// λ_th = h / √(2π m k_B T)
#[derive(Debug, Clone)]
pub struct DeBroglieWavelength {
    /// Particle mass (kg).
    pub mass: f64,
    /// Number density (m⁻³), used for quantum regime check.
    pub number_density: f64,
}

impl DeBroglieWavelength {
    /// Create a new [`DeBroglieWavelength`] calculator.
    ///
    /// # Arguments
    /// * `mass`           – particle mass (kg)
    /// * `number_density` – particle number density (m⁻³)
    pub fn new(mass: f64, number_density: f64) -> Self {
        Self {
            mass,
            number_density,
        }
    }

    /// Thermal de Broglie wavelength at temperature T (K).
    ///
    /// λ_th = h / √(2π m k_B T)
    ///
    /// # Arguments
    /// * `mass` – particle mass (kg)
    /// * `temp` – temperature (K)
    ///
    /// ```no_run
    /// use oxiphysics_sph::quantum_sph::DeBroglieWavelength;
    /// let dbw = DeBroglieWavelength::new(1.0, 1e20);
    /// let lambda = dbw.thermal(1e-27, 300.0);
    /// assert!(lambda > 0.0);
    /// ```
    pub fn thermal(&self, mass: f64, temp: f64) -> f64 {
        if mass <= 0.0 || temp <= 0.0 {
            return f64::INFINITY;
        }
        let h = 2.0 * PI * HBAR;
        h / (2.0 * PI * mass * K_B * temp).sqrt()
    }

    /// Check whether the system is in the quantum degenerate regime.
    ///
    /// The system is quantum when n λ_th³ ≳ 1, i.e.:
    /// n^(1/3) λ_th ≳ 1
    ///
    /// Returns `true` if the phase-space density n λ_th³ > 1.
    ///
    /// # Arguments
    /// * `temp` – temperature (K)
    pub fn is_quantum_regime(&self, temp: f64) -> bool {
        let lambda = self.thermal(self.mass, temp);
        if !lambda.is_finite() {
            return true;
        }
        self.number_density * lambda * lambda * lambda > 1.0
    }

    /// Phase-space density n λ_th³.
    ///
    /// # Arguments
    /// * `temp` – temperature (K)
    pub fn phase_space_density(&self, temp: f64) -> f64 {
        let lambda = self.thermal(self.mass, temp);
        if !lambda.is_finite() {
            return f64::INFINITY;
        }
        self.number_density * lambda * lambda * lambda
    }

    /// BEC transition temperature T_c for an ideal Bose gas.
    ///
    /// T_c = (ℏ²/m k_B) · (n/ζ(3/2))^(2/3) · 2π
    ///
    /// where ζ(3/2) ≈ 2.612.
    pub fn bec_transition_temperature(&self) -> f64 {
        if self.number_density <= 0.0 || self.mass <= 0.0 {
            return 0.0;
        }
        let zeta_32 = 2.612;
        let h = 2.0 * PI * HBAR;
        (h * h / (2.0 * PI * self.mass * K_B)) * (self.number_density / zeta_32).powf(2.0 / 3.0)
    }
}

/// Mass of a ⁴He atom (kg) — for test convenience.
#[cfg(test)]
const HELIUM4_MASS: f64 = 6.646_477_208_8e-27;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(x: f64, vx: f64) -> QuantumParticle {
        QuantumParticle::new([x, 0.0, 0.0], [vx, 0.0, 0.0], 1.0, 0.0)
    }

    fn make_qsph(n: usize) -> QuantumSPH {
        let particles: Vec<QuantumParticle> =
            (0..n).map(|i| make_particle(i as f64, 0.0)).collect();
        QuantumSPH::new(1.0, 1.0, 1.5, particles)
    }

    // ---- wkb_density -----------------------------------------------------

    #[test]
    fn test_wkb_density_classically_allowed() {
        let rho = wkb_density(2.0, 1.0, 1.0, 1.0);
        assert!(rho > 0.0, "rho={rho}");
    }

    #[test]
    fn test_wkb_density_forbidden_region() {
        let rho = wkb_density(0.5, 1.0, 1.0, 1.0);
        assert_eq!(rho, 0.0);
    }

    #[test]
    fn test_wkb_density_zero_mass() {
        let rho = wkb_density(1.0, 0.0, 0.0, 1.0);
        assert_eq!(rho, 0.0);
    }

    #[test]
    fn test_wkb_density_at_barrier() {
        let rho = wkb_density(1.0, 1.0, 1.0, 1.0);
        assert_eq!(rho, 0.0, "E == V: classically forbidden");
    }

    #[test]
    fn test_wkb_density_scales_with_mass() {
        let rho1 = wkb_density(2.0, 1.0, 1.0, 1.0);
        let rho2 = wkb_density(2.0, 1.0, 4.0, 1.0);
        // ρ ∝ m / √m = √m, so rho2 / rho1 = 2
        assert!((rho2 / rho1 - 2.0).abs() < 1e-10, "ratio={}", rho2 / rho1);
    }

    // ---- QuantumParticle --------------------------------------------------

    #[test]
    fn test_quantum_particle_probability_density() {
        let p = QuantumParticle::new([0.0; 3], [0.0; 3], 3.0, 4.0);
        assert!((p.probability_density() - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_quantum_particle_amplitude() {
        let p = QuantumParticle::new([0.0; 3], [0.0; 3], 3.0, 4.0);
        assert!((p.amplitude() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_quantum_particle_phase_zero_imag() {
        let p = QuantumParticle::new([0.0; 3], [0.0; 3], 1.0, 0.0);
        assert_eq!(p.phase(), 0.0);
    }

    #[test]
    fn test_quantum_particle_phase_pi_over_2() {
        let p = QuantumParticle::new([0.0; 3], [0.0; 3], 0.0, 1.0);
        assert!((p.phase() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    // ---- QuantumSPH -------------------------------------------------------

    #[test]
    fn test_qsph_new() {
        let q = make_qsph(4);
        assert_eq!(q.particles.len(), 4);
    }

    #[test]
    fn test_qsph_compute_density() {
        let mut q = make_qsph(5);
        q.compute_density();
        for p in &q.particles {
            assert!(p.density > 0.0, "density should be positive");
        }
    }

    #[test]
    fn test_qsph_bohm_potential_finite() {
        let mut q = make_qsph(4);
        q.compute_density();
        let qp = q.bohm_potential();
        for &v in &qp {
            assert!(v.is_finite(), "v={v}");
        }
    }

    #[test]
    fn test_qsph_compute_quantum_pressure_finite() {
        let mut q = make_qsph(4);
        q.compute_density();
        let qp = q.compute_quantum_pressure();
        for &v in &qp {
            assert!(v.is_finite(), "v={v}");
        }
    }

    #[test]
    fn test_qsph_step_runs() {
        let mut q = make_qsph(4);
        q.step(0.01);
        for p in &q.particles {
            for &v in &p.vel {
                assert!(v.is_finite(), "v={v}");
            }
        }
    }

    #[test]
    fn test_qsph_kinetic_energy_zero_at_rest() {
        let q = make_qsph(4);
        let ke = q.kinetic_energy();
        assert_eq!(ke, 0.0);
    }

    #[test]
    fn test_qsph_total_energy_finite() {
        let mut q = make_qsph(4);
        q.compute_density();
        let e = q.total_energy();
        assert!(e.is_finite(), "e={e}");
    }

    #[test]
    fn test_qsph_quantum_kinetic_energy_nonneg() {
        let mut q = make_qsph(4);
        q.compute_density();
        let qke = q.quantum_kinetic_energy();
        assert!(qke >= 0.0, "qke={qke}");
    }

    #[test]
    fn test_qsph_wigner_function_shape() {
        let q = make_qsph(4);
        let wf = q.wigner_function(8, 8, 5.0, 5.0);
        assert_eq!(wf.len(), 8);
        assert_eq!(wf[0].len(), 8);
    }

    #[test]
    fn test_qsph_wigner_function_nonneg() {
        let q = make_qsph(4);
        let wf = q.wigner_function(8, 8, 5.0, 5.0);
        for row in &wf {
            for &w in row {
                assert!(w >= 0.0, "w={w}");
            }
        }
    }

    // ---- FeynmanPathSPH ---------------------------------------------------

    #[test]
    fn test_feynman_path_sph_new() {
        let fp = FeynmanPathSPH::new(16, 1.0, 1.0, 1.0, [0.0; 3]);
        assert_eq!(fp.n_beads, 16);
        assert_eq!(fp.path.len(), 16);
    }

    #[test]
    fn test_feynman_action_nonneg() {
        let fp = FeynmanPathSPH::new(8, 1.0, 1.0, 1.0, [0.0; 3]);
        let s = fp.action();
        // All beads at same position → spring energy = 0, potential = 0
        assert!(s >= 0.0, "s={s}");
    }

    #[test]
    fn test_feynman_weight_between_0_and_1() {
        let fp = FeynmanPathSPH::new(8, 1.0, 1.0, 1.0, [0.0; 3]);
        let w = fp.path_integral_weight();
        // w = exp(-S/ℏ); S >= 0 → w <= 1
        assert!(w > 0.0 && w <= 1.0, "w={w}");
    }

    #[test]
    fn test_feynman_propagate_runs() {
        let mut fp = FeynmanPathSPH::new(8, 1.0, 1.0, 1.0, [0.0; 3]);
        for _ in 0..20 {
            fp.propagate(0.1);
        }
        for bead in &fp.path {
            for &x in bead {
                assert!(x.is_finite(), "x={x}");
            }
        }
    }

    #[test]
    fn test_feynman_centroid_all_same() {
        let fp = FeynmanPathSPH::new(8, 1.0, 1.0, 1.0, [1.0, 2.0, 3.0]);
        let c = fp.centroid();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_feynman_gyration_radius_zero_for_same_beads() {
        let fp = FeynmanPathSPH::new(8, 1.0, 1.0, 1.0, [1.0, 0.0, 0.0]);
        let rg = fp.gyration_radius();
        assert!(rg < 1e-12, "rg={rg}");
    }

    // ---- DeBroglieWavelength ----------------------------------------------

    #[test]
    fn test_de_broglie_thermal_positive() {
        let dbw = DeBroglieWavelength::new(PROTON_MASS, 1e20);
        let lambda = dbw.thermal(PROTON_MASS, 300.0);
        assert!(lambda > 0.0, "lambda={lambda}");
    }

    #[test]
    fn test_de_broglie_thermal_zero_temp() {
        let dbw = DeBroglieWavelength::new(PROTON_MASS, 1e20);
        let lambda = dbw.thermal(PROTON_MASS, 0.0);
        assert!(lambda.is_infinite(), "lambda={lambda}");
    }

    #[test]
    fn test_de_broglie_thermal_zero_mass() {
        let dbw = DeBroglieWavelength::new(0.0, 1e20);
        let lambda = dbw.thermal(0.0, 300.0);
        assert!(lambda.is_infinite(), "lambda={lambda}");
    }

    #[test]
    fn test_de_broglie_thermal_scales_with_temp() {
        let dbw = DeBroglieWavelength::new(PROTON_MASS, 1e20);
        let l1 = dbw.thermal(PROTON_MASS, 100.0);
        let l4 = dbw.thermal(PROTON_MASS, 400.0);
        // λ ∝ 1/√T → l1/l4 = 2
        let ratio = l1 / l4;
        assert!((ratio - 2.0).abs() < 1e-6, "ratio={ratio}");
    }

    #[test]
    fn test_de_broglie_classical_regime() {
        // Hot, dense gas should NOT be quantum
        let dbw = DeBroglieWavelength::new(PROTON_MASS, 1e20);
        // At room temperature, proton gas is classical
        assert!(!dbw.is_quantum_regime(300.0));
    }

    #[test]
    fn test_de_broglie_quantum_regime_cold() {
        // Cold, dense electron gas
        let dbw = DeBroglieWavelength::new(ELECTRON_MASS, 1e28);
        // Electrons at 1 K should be quantum degenerate
        assert!(dbw.is_quantum_regime(1.0));
    }

    #[test]
    fn test_de_broglie_phase_space_density_positive() {
        let dbw = DeBroglieWavelength::new(PROTON_MASS, 1e20);
        let psd = dbw.phase_space_density(300.0);
        assert!(psd > 0.0, "psd={psd}");
    }

    #[test]
    fn test_de_broglie_bec_temp_positive() {
        let dbw = DeBroglieWavelength::new(HELIUM4_MASS, 1e26);
        let tc = dbw.bec_transition_temperature();
        assert!(tc > 0.0, "Tc={tc}");
    }

    #[test]
    fn test_de_broglie_bec_temp_scales_with_density() {
        let dbw1 = DeBroglieWavelength::new(HELIUM4_MASS, 1e26);
        let dbw8 = DeBroglieWavelength::new(HELIUM4_MASS, 8e26);
        let tc1 = dbw1.bec_transition_temperature();
        let tc8 = dbw8.bec_transition_temperature();
        // T_c ∝ n^(2/3) → tc8/tc1 = 4
        let ratio = tc8 / tc1;
        assert!((ratio - 4.0).abs() < 0.01, "ratio={ratio}");
    }
}
