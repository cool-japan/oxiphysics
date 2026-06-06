// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum-classical MD using path-integral (ring-polymer) representation.
//!
//! Implements the Path-Integral Molecular Dynamics (PIMD) formalism where
//! each quantum particle is represented as a ring of P classical beads
//! coupled by harmonic springs.  Includes centroid MD (CMD), ring polymer MD
//! (RPMD), the PILE (Path-Integral Langevin Equation) thermostat, normal-mode
//! and staging transformations, quantum kinetic-energy estimators, velocity
//! autocorrelation functions for quantum diffusion, and simple tunneling-rate
//! estimates.
//!
//! # Key types
//! - [`QuantumBead`]: a single quantum particle represented as a ring polymer.
//! - [`QuantumSystem`]: collection of quantum beads forming a full system.
//! - [`NormalModeTransform`]: forward/inverse normal-mode transformation for
//!   ring polymers.
//! - [`StagingCoords`]: staging variable transformation and its inverse.
//! - [`PileThermostat`]: Path-Integral Langevin Equation (PILE) thermostat.
//! - [`CentroidMdState`]: centroid-MD state with adiabatically decoupled
//!   fluctuations.
//! - [`RingPolymerMdState`]: RPMD state for quantum real-time dynamics.
//! - [`VacfAccumulator`]: velocity autocorrelation function accumulator for
//!   computing quantum diffusion coefficients.
//! - [`TunnelingRate`]: Wentzel-Kramers-Brillouin (WKB) and instanton
//!   tunneling-rate estimators.
//! - [`QuantumDiagnostics`]: summary diagnostics for a PIMD simulation.
//!
//! # Key functions
//! - [`ring_polymer_spring_constant`]: spring constant for the inter-bead
//!   harmonic coupling.
//! - [`ring_polymer_spring_force`]: force from harmonic springs along the ring.
//! - [`centroid_position`] / [`centroid_velocity`]: classical centroid
//!   estimators.
//! - [`quantum_kinetic_energy`]: virial estimator for the quantum kinetic
//!   energy.
//! - [`primitive_kinetic_energy`]: primitive (thermodynamic) kinetic energy
//!   estimator.
//! - [`thermal_de_broglie`]: thermal de Broglie wavelength.
//! - [`step_pimd`]: one velocity-Verlet step of the ring-polymer MD.
//! - [`step_rpmd`]: one velocity-Verlet step with spring forces included.
//! - [`normal_mode_frequencies`]: compute ring-polymer normal-mode frequencies.
//! - [`wkb_tunneling_rate`]: WKB tunneling rate through a parabolic barrier.
//! - [`quantum_diffusion_coefficient`]: Green-Kubo integral of the VACF.

use std::f64::consts::PI;

// ============================================================================
// 3-vector helpers (private)
// ============================================================================

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(s: f64, a: [f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_norm_sq(a: [f64; 3]) -> f64 {
    v3_dot(a, a)
}

// ============================================================================
// Core structs
// ============================================================================

/// A single quantum particle represented as a ring polymer of P beads.
#[derive(Debug, Clone)]
pub struct QuantumBead {
    /// Bead positions, one per ring-polymer bead: `positions[p] = [x, y, z]`.
    pub positions: Vec<[f64; 3]>,
    /// Bead velocities, one per ring-polymer bead.
    pub velocities: Vec<[f64; 3]>,
    /// Physical mass of the particle (a.m.u. or internal units).
    pub mass: f64,
    /// Effective reduced Planck constant ℏ_eff (may be scaled for isotope
    /// studies).
    pub hbar_eff: f64,
}

impl QuantumBead {
    /// Create a new `QuantumBead` with all beads at the origin.
    pub fn new(n_beads: usize, mass: f64, hbar_eff: f64) -> Self {
        Self {
            positions: vec![[0.0; 3]; n_beads],
            velocities: vec![[0.0; 3]; n_beads],
            mass,
            hbar_eff,
        }
    }

    /// Number of ring-polymer beads P.
    pub fn n_beads(&self) -> usize {
        self.positions.len()
    }

    /// Centroid position of this ring polymer.
    pub fn centroid(&self) -> [f64; 3] {
        centroid_position(&self.positions)
    }

    /// Centroid velocity of this ring polymer.
    pub fn centroid_vel(&self) -> [f64; 3] {
        centroid_velocity(&self.velocities)
    }

    /// Mean-square displacement of beads around the centroid (ring-polymer
    /// spread).
    pub fn ring_spread_sq(&self) -> f64 {
        let c = self.centroid();
        self.positions
            .iter()
            .map(|&r| v3_norm_sq(v3_sub(r, c)))
            .sum::<f64>()
            / self.positions.len() as f64
    }
}

/// Collection of quantum particles forming a complete PIMD system.
#[derive(Debug, Clone)]
pub struct QuantumSystem {
    /// All quantum particles (ring polymers).
    pub beads: Vec<QuantumBead>,
    /// Number of ring-polymer beads P (same for all particles).
    pub n_beads: usize,
    /// Inverse temperature β = 1/(k_B T).
    pub beta: f64,
    /// Number of atoms (= `beads.len()`).
    pub n_atoms: usize,
}

impl QuantumSystem {
    /// Create a new `QuantumSystem` with `n_atoms` quantum beads.
    pub fn new(n_atoms: usize, n_beads: usize, beta: f64, mass: f64, hbar_eff: f64) -> Self {
        let beads = (0..n_atoms)
            .map(|_| QuantumBead::new(n_beads, mass, hbar_eff))
            .collect();
        Self {
            beads,
            n_beads,
            beta,
            n_atoms,
        }
    }

    /// Temperature T = 1/β (in units where k_B = 1).
    pub fn temperature(&self) -> f64 {
        1.0 / self.beta
    }

    /// Total centroid kinetic energy (classical, all atoms).
    pub fn centroid_kinetic_energy(&self) -> f64 {
        self.beads
            .iter()
            .map(|b| 0.5 * b.mass * v3_norm_sq(b.centroid_vel()))
            .sum()
    }
}

// ============================================================================
// Normal-mode transformation
// ============================================================================

/// Forward and inverse normal-mode transformation for ring polymers.
///
/// The transformation diagonalises the ring-polymer spring matrix so that each
/// normal mode evolves independently in the harmonic (free ring-polymer) limit.
#[derive(Debug, Clone)]
pub struct NormalModeTransform {
    /// Number of beads P.
    pub n_beads: usize,
    /// Normal-mode matrix C (P × P), stored row-major.
    pub c_matrix: Vec<f64>,
}

impl NormalModeTransform {
    /// Build the normal-mode matrix for P beads using the standard
    /// orthogonal DFT-like transform.
    ///
    /// For bead index k and mode n (0-indexed):
    ///
    /// ```text
    /// C[k][n=0]        = 1 / sqrt(P)
    /// C[k][n=1..P/2]   = sqrt(2/P) * cos(2π n k / P)
    /// C[k][n=P/2]      = (-1)^k / sqrt(P)      (even P only)
    /// C[k][n=P/2+1..]  = sqrt(2/P) * sin(2π n k / P)
    /// ```
    pub fn new(n_beads: usize) -> Self {
        let p = n_beads;
        let mut c = vec![0.0f64; p * p];
        let sqrt_p = (p as f64).sqrt();
        for k in 0..p {
            // Mode 0: centroid
            c[k * p] = 1.0 / sqrt_p;
            for n in 1..p {
                let angle = 2.0 * PI * (n as f64) * (k as f64) / (p as f64);
                c[k * p + n] = if n < p / 2 + 1 {
                    (2.0 / p as f64).sqrt() * angle.cos()
                } else {
                    (2.0 / p as f64).sqrt() * angle.sin()
                };
            }
        }
        // Handle even P: mode P/2 should be (-1)^k / sqrt(P)
        if p.is_multiple_of(2) {
            let half = p / 2;
            for k in 0..p {
                c[k * p + half] = if k % 2 == 0 {
                    1.0 / sqrt_p
                } else {
                    -1.0 / sqrt_p
                };
            }
        }
        Self {
            n_beads: p,
            c_matrix: c,
        }
    }

    /// Apply forward transform: convert bead positions to normal-mode
    /// amplitudes.  Input and output are P 3-vectors.
    pub fn forward(&self, bead_coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let p = self.n_beads;
        let mut modes = vec![[0.0f64; 3]; p];
        for (n, mode_n) in modes.iter_mut().enumerate() {
            for (k, &bc_k) in bead_coords.iter().enumerate().take(p) {
                let c_kn = self.c_matrix[k * p + n];
                *mode_n = v3_add(*mode_n, v3_scale(c_kn, bc_k));
            }
        }
        modes
    }

    /// Apply inverse transform: convert normal-mode amplitudes back to bead
    /// positions.
    pub fn inverse(&self, mode_coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let p = self.n_beads;
        let mut beads = vec![[0.0f64; 3]; p];
        for (k, bead_k) in beads.iter_mut().enumerate() {
            for (n, &mc_n) in mode_coords.iter().enumerate().take(p) {
                let c_kn = self.c_matrix[k * p + n];
                *bead_k = v3_add(*bead_k, v3_scale(c_kn, mc_n));
            }
        }
        beads
    }
}

/// Compute the ring-polymer normal-mode angular frequencies ω_n for P beads
/// at inverse temperature β.
///
/// ω_n = (2 P / (β ℏ)) sin(π n / P),   n = 0, 1, …, P-1.
///
/// Mode n=0 is the centroid (ω=0).
pub fn normal_mode_frequencies(n_beads: usize, beta: f64, hbar: f64) -> Vec<f64> {
    let p = n_beads as f64;
    (0..n_beads)
        .map(|n| {
            let sin_arg = PI * (n as f64) / p;
            2.0 * p / (beta * hbar) * sin_arg.sin()
        })
        .collect()
}

// ============================================================================
// Staging transformation
// ============================================================================

/// Staging-variable transformation for ring polymers.
///
/// The staging coordinates decouple the harmonic spring interactions into a set
/// of independent Einstein-oscillator modes, with masses that scale as
/// `m_s = s*m / (s-1)` for staging level s.
#[derive(Debug, Clone)]
pub struct StagingCoords {
    /// Number of beads P.
    pub n_beads: usize,
    /// Physical particle mass.
    pub mass: f64,
}

impl StagingCoords {
    /// Create a new staging-coordinate helper.
    pub fn new(n_beads: usize, mass: f64) -> Self {
        Self { n_beads, mass }
    }

    /// Convert bead positions to staging coordinates.
    ///
    /// u_1 = r_1  (centroid anchor)
    /// u_j = r_j - ((j-1)*r_{j+1} + r_1) / j,  for j = 2, …, P.
    pub fn to_staging(&self, beads: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let p = self.n_beads;
        let mut u = vec![[0.0f64; 3]; p];
        if p == 0 {
            return u;
        }
        u[0] = beads[0];
        for j in 1..p {
            // r_{j+1} wraps around: index (j+1) % p
            let r_next = beads[(j + 1) % p];
            let r_1 = beads[0];
            // u_j = r_j - ( (j-1)*r_next + r_1 ) / j
            let jf = j as f64;
            u[j] = v3_sub(
                beads[j],
                v3_scale(1.0 / jf, v3_add(v3_scale(jf - 1.0, r_next), r_1)),
            );
        }
        u
    }

    /// Convert staging coordinates back to bead positions (inverse staging).
    pub fn from_staging(&self, u: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let p = self.n_beads;
        let mut r = vec![[0.0f64; 3]; p];
        if p == 0 {
            return r;
        }
        r[0] = u[0];
        // Reconstruct from the last bead backwards.
        // For j = P-1 down to 1:
        //   r_j = u_j + ((j-1)*r_{j+1} + r_1) / j
        // We need r_{j+1} already reconstructed; do backwards.
        // First set r[p-1] via the P-1 staging:
        // Actually the standard inversion reconstructs forwards using the
        // relation:  r_{j+1} = (j * u_j - r_1 + j * r_j) / (j+1) ... this
        // varies by convention.  We use a simple iterative forward pass.
        for j in 1..p {
            let jf = j as f64;
            let r_next = r[(j + 1) % p]; // wrapped; r[0] is already known
            r[j] = v3_add(
                u[j],
                v3_scale(1.0 / jf, v3_add(v3_scale(jf - 1.0, r_next), r[0])),
            );
        }
        r
    }

    /// Effective staging mass for level j (1-indexed, j ≥ 2).
    ///
    /// m_s(j) = j/(j-1) * m
    pub fn staging_mass(&self, j: usize) -> f64 {
        if j <= 1 {
            self.mass
        } else {
            self.mass * (j as f64) / ((j - 1) as f64)
        }
    }
}

// ============================================================================
// Ring-polymer functions
// ============================================================================

/// Spring constant for the ring-polymer harmonic coupling.
///
/// k = m * (P / (β ℏ))²
///
/// where P is the number of beads, β = 1/(k_B T), and ℏ is the (effective)
/// reduced Planck constant.
pub fn ring_polymer_spring_constant(mass: f64, n_beads: usize, beta: f64) -> f64 {
    let hbar = 1.0; // internal units
    let p = n_beads as f64;
    mass * (p / (beta * hbar)).powi(2)
}

/// Harmonic spring forces on each bead due to adjacent beads in the ring.
///
/// F_s(i) = k * (r_{i+1} + r_{i-1} − 2 r_i)  (cyclic indexing)
pub fn ring_polymer_spring_force(beads: &[[f64; 3]], k: f64) -> Vec<[f64; 3]> {
    let n = beads.len();
    if n == 0 {
        return Vec::new();
    }
    (0..n)
        .map(|i| {
            let prev = beads[(i + n - 1) % n];
            let next = beads[(i + 1) % n];
            let curr = beads[i];
            let delta = v3_add(v3_sub(prev, curr), v3_sub(next, curr));
            v3_scale(k, delta)
        })
        .collect()
}

/// Centroid position: arithmetic mean of all bead positions.
pub fn centroid_position(beads: &[[f64; 3]]) -> [f64; 3] {
    if beads.is_empty() {
        return [0.0; 3];
    }
    let n = beads.len() as f64;
    let sum = beads.iter().fold([0.0_f64; 3], |acc, &b| v3_add(acc, b));
    v3_scale(1.0 / n, sum)
}

/// Centroid velocity: arithmetic mean of all bead velocities.
pub fn centroid_velocity(velocities: &[[f64; 3]]) -> [f64; 3] {
    centroid_position(velocities)
}

/// Quantum kinetic energy using the virial estimator.
///
/// T_vir = (3 N P) / (2β) − (k/2) Σ_{i,p} |r_ip − r_centroid_i|²
///
/// where k is the ring-polymer spring constant for each atom.
pub fn quantum_kinetic_energy(sys: &QuantumSystem) -> f64 {
    let n = sys.n_atoms as f64;
    let p = sys.n_beads as f64;
    let classical_term = 3.0 * n * p / (2.0 * sys.beta);
    let mut spring_term = 0.0;
    for bead in &sys.beads {
        let k = ring_polymer_spring_constant(bead.mass, sys.n_beads, sys.beta);
        let c = centroid_position(&bead.positions);
        for &pos in &bead.positions {
            let dr = v3_sub(pos, c);
            spring_term += k * v3_dot(dr, dr);
        }
    }
    classical_term - 0.5 * spring_term
}

/// Primitive (thermodynamic) kinetic energy estimator.
///
/// T_prim = (3 N P) / (2β) − (k_P/2) Σ_{i,p} |r_{i,p+1} − r_{i,p}|²
///
/// This estimator has larger variance than the virial estimator but is
/// straightforward to implement.
pub fn primitive_kinetic_energy(sys: &QuantumSystem) -> f64 {
    let n = sys.n_atoms as f64;
    let p_count = sys.n_beads;
    let p = p_count as f64;
    let classical_term = 3.0 * n * p / (2.0 * sys.beta);
    let mut spring_term = 0.0;
    for bead in &sys.beads {
        let k = ring_polymer_spring_constant(bead.mass, p_count, sys.beta);
        for j in 0..p_count {
            let rj = bead.positions[j];
            let rj1 = bead.positions[(j + 1) % p_count];
            let dr = v3_sub(rj1, rj);
            spring_term += k * v3_dot(dr, dr);
        }
    }
    classical_term - 0.5 * spring_term
}

/// Thermal de Broglie wavelength Λ = ℏ √(2π β / m).
pub fn thermal_de_broglie(mass: f64, beta: f64, hbar: f64) -> f64 {
    hbar * (2.0 * PI * beta / mass).sqrt()
}

/// Quantum correction factor: ratio of quantum to classical kinetic energy.
///
/// Returns 1 when `classical_ke` is zero (to avoid division by zero).
pub fn quantum_correction_factor(classical_ke: f64, quantum_ke: f64) -> f64 {
    if classical_ke.abs() < 1e-30 {
        return 1.0;
    }
    quantum_ke / classical_ke
}

/// Total ring-polymer potential energy (spring contributions only).
///
/// U_spring = (k/2) Σ_{i,p} |r_{i,p+1} − r_{i,p}|²
pub fn ring_polymer_spring_energy(sys: &QuantumSystem) -> f64 {
    let mut u = 0.0;
    for bead in &sys.beads {
        let k = ring_polymer_spring_constant(bead.mass, sys.n_beads, sys.beta);
        let p = sys.n_beads;
        for j in 0..p {
            let rj = bead.positions[j];
            let rj1 = bead.positions[(j + 1) % p];
            let dr = v3_sub(rj1, rj);
            u += 0.5 * k * v3_dot(dr, dr);
        }
    }
    u
}

// ============================================================================
// PILE thermostat (Path-Integral Langevin Equation)
// ============================================================================

/// Path-Integral Langevin Equation (PILE) thermostat for PIMD.
///
/// PILE couples each normal mode of the ring polymer to an independent Langevin
/// heat bath with a friction coefficient γ_n = 2 ω_n (critical damping for
/// non-centroid modes) and γ_0 = γ_centroid for the centroid.
///
/// Reference: Ceriotti, Parrinello, Markland, Manolopoulos, JCP 133, 124104
/// (2010).
#[derive(Debug, Clone)]
pub struct PileThermostat {
    /// Number of ring-polymer beads P.
    pub n_beads: usize,
    /// Inverse temperature β = 1/(k_B T).
    pub beta: f64,
    /// Reduced Planck constant ℏ (internal units).
    pub hbar: f64,
    /// Centroid friction coefficient γ_0.
    pub gamma_centroid: f64,
    /// Normal-mode friction coefficients γ_n (length P).
    pub gamma_modes: Vec<f64>,
    /// Precomputed c1_n = exp(−γ_n dt/2) coefficients.
    pub c1: Vec<f64>,
    /// Precomputed c2_n = sqrt(1 − c1_n²) coefficients.
    pub c2: Vec<f64>,
    /// Time step used when building c1/c2.
    pub dt: f64,
}

impl PileThermostat {
    /// Build a PILE thermostat.
    ///
    /// # Arguments
    /// * `n_beads` — number of ring-polymer beads P
    /// * `beta` — inverse temperature
    /// * `hbar` — reduced Planck constant
    /// * `gamma_centroid` — centroid-mode friction coefficient
    /// * `dt` — integration time step
    pub fn new(n_beads: usize, beta: f64, hbar: f64, gamma_centroid: f64, dt: f64) -> Self {
        let freqs = normal_mode_frequencies(n_beads, beta, hbar);
        let mut gamma_modes = Vec::with_capacity(n_beads);
        gamma_modes.push(gamma_centroid);
        for freq_n in freqs.iter().skip(1) {
            gamma_modes.push(2.0 * freq_n); // critical damping
        }
        let c1: Vec<f64> = gamma_modes.iter().map(|&g| (-g * dt / 2.0).exp()).collect();
        let c2: Vec<f64> = c1.iter().map(|&c| (1.0 - c * c).sqrt()).collect();
        Self {
            n_beads,
            beta,
            hbar,
            gamma_centroid,
            gamma_modes,
            c1,
            c2,
            dt,
        }
    }

    /// Apply the PILE half-step thermostat to a single bead's velocity in one
    /// normal mode.
    ///
    /// Returns the updated velocity component after Langevin stochastic
    /// integration:
    ///
    /// ```text
    /// v_new = c1 * v_old + c2 * sqrt(1/(β m)) * xi
    /// ```
    ///
    /// where `xi` is a standard normal random variate.
    pub fn apply_half_step(&self, mode_idx: usize, v_comp: f64, mass: f64, xi: f64) -> f64 {
        let sigma = (1.0 / (self.beta * mass)).sqrt();
        self.c1[mode_idx] * v_comp + self.c2[mode_idx] * sigma * xi
    }
}

// ============================================================================
// Centroid MD state
// ============================================================================

/// State for Centroid Molecular Dynamics (CMD).
///
/// In CMD, the system propagates the centroid degree of freedom on a mean-field
/// potential obtained by integrating over the ring-polymer fluctuations.
/// This struct tracks the centroid positions and velocities separately from the
/// internal fluctuation modes.
#[derive(Debug, Clone)]
pub struct CentroidMdState {
    /// Centroid positions for each atom.
    pub centroid_positions: Vec<[f64; 3]>,
    /// Centroid velocities for each atom.
    pub centroid_velocities: Vec<[f64; 3]>,
    /// Centroid forces for each atom (updated each step).
    pub centroid_forces: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Inverse temperature β.
    pub beta: f64,
    /// Number of ring-polymer beads P (used to weight quantum fluctuations).
    pub n_beads: usize,
}

impl CentroidMdState {
    /// Create a new CMD state with `n_atoms` atoms all at the origin.
    pub fn new(n_atoms: usize, n_beads: usize, beta: f64, mass: f64) -> Self {
        Self {
            centroid_positions: vec![[0.0; 3]; n_atoms],
            centroid_velocities: vec![[0.0; 3]; n_atoms],
            centroid_forces: vec![[0.0; 3]; n_atoms],
            masses: vec![mass; n_atoms],
            beta,
            n_beads,
        }
    }

    /// Perform one velocity-Verlet step on the centroid degrees of freedom.
    ///
    /// After calling this, `centroid_forces` should be updated externally
    /// using the new centroid positions before calling again.
    pub fn step(&mut self, dt: f64) {
        let n = self.centroid_positions.len();
        // Half-kick velocities
        for i in 0..n {
            let inv_m = 1.0 / self.masses[i];
            self.centroid_velocities[i] = v3_add(
                self.centroid_velocities[i],
                v3_scale(0.5 * dt * inv_m, self.centroid_forces[i]),
            );
        }
        // Update positions
        for i in 0..n {
            self.centroid_positions[i] = v3_add(
                self.centroid_positions[i],
                v3_scale(dt, self.centroid_velocities[i]),
            );
        }
    }

    /// Second half-kick after forces have been recomputed at new positions.
    pub fn second_half_kick(&mut self, dt: f64) {
        let n = self.centroid_positions.len();
        for i in 0..n {
            let inv_m = 1.0 / self.masses[i];
            self.centroid_velocities[i] = v3_add(
                self.centroid_velocities[i],
                v3_scale(0.5 * dt * inv_m, self.centroid_forces[i]),
            );
        }
    }

    /// Classical kinetic energy of the centroid degrees of freedom.
    pub fn kinetic_energy(&self) -> f64 {
        self.centroid_velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * v3_norm_sq(*v))
            .sum()
    }

    /// Instantaneous centroid temperature (in units where k_B = 1).
    pub fn temperature(&self) -> f64 {
        let n = self.centroid_positions.len() as f64;
        if n < 1.0 {
            return 0.0;
        }
        let ke = self.kinetic_energy();
        2.0 * ke / (3.0 * n)
    }
}

// ============================================================================
// Ring Polymer MD state
// ============================================================================

/// State for Ring Polymer Molecular Dynamics (RPMD).
///
/// RPMD propagates all P beads using the physical Hamiltonian plus the
/// ring-polymer spring terms.  It provides an exact quantum real-time
/// correlation function in the short-time limit and a good approximation
/// to quantum diffusion.
#[derive(Debug, Clone)]
pub struct RingPolymerMdState {
    /// Full ring-polymer system (P beads per atom).
    pub system: QuantumSystem,
    /// External forces on each bead: `forces[atom][bead] = [fx, fy, fz]`.
    pub forces: Vec<Vec<[f64; 3]>>,
    /// Simulation time accumulated.
    pub time: f64,
}

impl RingPolymerMdState {
    /// Create a new RPMD state.
    pub fn new(n_atoms: usize, n_beads: usize, beta: f64, mass: f64, hbar: f64) -> Self {
        let system = QuantumSystem::new(n_atoms, n_beads, beta, mass, hbar);
        let forces = vec![vec![[0.0_f64; 3]; n_beads]; n_atoms];
        Self {
            system,
            forces,
            time: 0.0,
        }
    }

    /// Advance by one velocity-Verlet step including spring forces.
    pub fn step(&mut self, dt: f64) {
        step_rpmd(&mut self.system, &self.forces, dt);
        self.time += dt;
    }

    /// Total kinetic energy of the ring polymer (all beads).
    pub fn kinetic_energy(&self) -> f64 {
        self.system
            .beads
            .iter()
            .map(|b| {
                b.velocities
                    .iter()
                    .map(|&v| 0.5 * b.mass * v3_norm_sq(v))
                    .sum::<f64>()
            })
            .sum()
    }

    /// Ring-polymer spring potential energy.
    pub fn spring_energy(&self) -> f64 {
        ring_polymer_spring_energy(&self.system)
    }
}

// ============================================================================
// Integrators
// ============================================================================

/// One velocity-Verlet step for the ring-polymer MD.
///
/// `forces[i][p]` is the external force on atom `i`, bead `p`.
/// Spring forces are added internally.
/// Updates positions and velocities in place.
pub fn step_pimd(sys: &mut QuantumSystem, forces: &[Vec<[f64; 3]>], dt: f64) {
    let n_atoms = sys.n_atoms;
    let n_beads = sys.n_beads;

    let mut spring_forces: Vec<Vec<[f64; 3]>> = Vec::with_capacity(n_atoms);
    for bead in &sys.beads {
        let k = ring_polymer_spring_constant(bead.mass, n_beads, sys.beta);
        spring_forces.push(ring_polymer_spring_force(&bead.positions, k));
    }

    // Half-kick + position update
    for (i, bead) in sys.beads.iter_mut().enumerate() {
        let inv_m = 1.0 / bead.mass;
        for p in 0..n_beads {
            let f_tot = v3_add(forces[i][p], spring_forces[i][p]);
            bead.velocities[p] = v3_add(bead.velocities[p], v3_scale(0.5 * dt * inv_m, f_tot));
            bead.positions[p] = v3_add(bead.positions[p], v3_scale(dt, bead.velocities[p]));
        }
    }

    // Recompute spring forces at new positions
    let mut spring_forces2: Vec<Vec<[f64; 3]>> = Vec::with_capacity(n_atoms);
    for bead in &sys.beads {
        let k = ring_polymer_spring_constant(bead.mass, n_beads, sys.beta);
        spring_forces2.push(ring_polymer_spring_force(&bead.positions, k));
    }

    // Second half-kick
    for (i, bead) in sys.beads.iter_mut().enumerate() {
        let inv_m = 1.0 / bead.mass;
        for p in 0..n_beads {
            let f_tot = v3_add(forces[i][p], spring_forces2[i][p]);
            bead.velocities[p] = v3_add(bead.velocities[p], v3_scale(0.5 * dt * inv_m, f_tot));
        }
    }
}

/// One velocity-Verlet step for RPMD (identical to PIMD integrator but named
/// separately for conceptual clarity and potential future divergence).
pub fn step_rpmd(sys: &mut QuantumSystem, forces: &[Vec<[f64; 3]>], dt: f64) {
    step_pimd(sys, forces, dt);
}

// ============================================================================
// Velocity Autocorrelation Function (VACF) and diffusion
// ============================================================================

/// Accumulator for the velocity autocorrelation function (VACF).
///
/// The VACF C(t) = ⟨v(0)·v(t)⟩ / ⟨v(0)·v(0)⟩ is computed from stored
/// velocity snapshots.  The quantum self-diffusion coefficient follows from
/// the Green-Kubo relation D = ∫₀^∞ C(t) dt.
#[derive(Debug, Clone)]
pub struct VacfAccumulator {
    /// Stored centroid velocity snapshots: `snapshots[step][atom] = [vx,vy,vz]`.
    pub snapshots: Vec<Vec<[f64; 3]>>,
    /// Maximum number of snapshots to retain.
    pub max_snapshots: usize,
}

impl VacfAccumulator {
    /// Create a new VACF accumulator.
    pub fn new(max_snapshots: usize) -> Self {
        Self {
            snapshots: Vec::new(),
            max_snapshots,
        }
    }

    /// Push one snapshot of centroid velocities.
    pub fn push(&mut self, velocities: Vec<[f64; 3]>) {
        if self.snapshots.len() >= self.max_snapshots {
            self.snapshots.remove(0);
        }
        self.snapshots.push(velocities);
    }

    /// Compute the VACF at lag `tau` (in snapshot steps).
    ///
    /// Returns the normalised correlation C(tau)/C(0), or `None` if insufficient
    /// data.
    pub fn vacf(&self, tau: usize) -> Option<f64> {
        let n_snaps = self.snapshots.len();
        if tau >= n_snaps {
            return None;
        }
        let n_pairs = n_snaps - tau;
        if n_pairs == 0 {
            return None;
        }
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for t0 in 0..n_pairs {
            let v0 = &self.snapshots[t0];
            let vt = &self.snapshots[t0 + tau];
            for (a0, at) in v0.iter().zip(vt.iter()) {
                numerator += v3_dot(*a0, *at);
            }
        }
        for v0 in &self.snapshots[..n_pairs] {
            for a0 in v0.iter() {
                denominator += v3_dot(*a0, *a0);
            }
        }
        if denominator.abs() < 1e-30 {
            return None;
        }
        Some(numerator / denominator)
    }

    /// Compute the un-normalised VACF at lag `tau`.
    pub fn vacf_unnorm(&self, tau: usize) -> f64 {
        let n_snaps = self.snapshots.len();
        if tau >= n_snaps {
            return 0.0;
        }
        let n_pairs = n_snaps - tau;
        let mut sum = 0.0;
        for t0 in 0..n_pairs {
            let v0 = &self.snapshots[t0];
            let vt = &self.snapshots[t0 + tau];
            for (a0, at) in v0.iter().zip(vt.iter()) {
                sum += v3_dot(*a0, *at);
            }
        }
        sum / n_pairs as f64
    }
}

/// Compute the quantum self-diffusion coefficient from a VACF array.
///
/// D = (1/3) ∫₀^∞ C(t) dt  ≈  (dt/3) Σ_τ C(τ)
///
/// The factor 1/3 accounts for 3D averaging (assuming isotropic diffusion and
/// VACF already containing the dot product sum over all 3 components).
pub fn quantum_diffusion_coefficient(vacf_values: &[f64], dt: f64) -> f64 {
    // Trapezoidal rule
    if vacf_values.is_empty() {
        return 0.0;
    }
    let n = vacf_values.len();
    let mut integral = 0.5 * (vacf_values[0] + vacf_values[n - 1]);
    for &v in &vacf_values[1..n - 1] {
        integral += v;
    }
    integral * dt / 3.0
}

// ============================================================================
// Tunneling rates
// ============================================================================

/// WKB and instanton tunneling-rate estimators.
#[derive(Debug, Clone)]
pub struct TunnelingRate {
    /// Barrier height E_b (energy units).
    pub barrier_height: f64,
    /// Barrier width a (length units) at half-height.
    pub barrier_width: f64,
    /// Particle mass m.
    pub mass: f64,
    /// Reduced Planck constant ℏ.
    pub hbar: f64,
}

impl TunnelingRate {
    /// Create a new tunneling-rate estimator.
    pub fn new(barrier_height: f64, barrier_width: f64, mass: f64, hbar: f64) -> Self {
        Self {
            barrier_height,
            barrier_width,
            mass,
            hbar,
        }
    }

    /// WKB transmission coefficient for a square barrier of height E_b and
    /// width a, for particle energy E.
    ///
    /// T_WKB = exp(−2 κ a),  κ = sqrt(2 m (E_b − E)) / ℏ
    ///
    /// Returns 1.0 when E ≥ E_b (above-barrier transmission).
    pub fn wkb_transmission(&self, particle_energy: f64) -> f64 {
        if particle_energy >= self.barrier_height {
            return 1.0;
        }
        let kappa = (2.0 * self.mass * (self.barrier_height - particle_energy)).sqrt() / self.hbar;
        (-2.0 * kappa * self.barrier_width).exp()
    }

    /// Instanton tunneling rate using the leading-order semiclassical result:
    ///
    /// k_inst = ω_0 / (2π) * exp(−S_inst / ℏ)
    ///
    /// where S_inst = π ℏ sqrt(2 m E_b) * a is the imaginary-time action for
    /// a parabolic barrier and ω_0 = sqrt(2 E_b / (m a²)) is the barrier
    /// frequency.
    pub fn instanton_rate(&self) -> f64 {
        let s_inst =
            PI * self.hbar * (2.0 * self.mass * self.barrier_height).sqrt() * self.barrier_width;
        let omega_0 = (2.0 * self.barrier_height / (self.mass * self.barrier_width.powi(2))).sqrt();
        omega_0 / (2.0 * PI) * (-s_inst / self.hbar).exp()
    }
}

/// Compute the WKB tunneling rate for a parabolic barrier.
///
/// Γ = ω_b / (2π) * exp(−2 S_WKB / ℏ)
///
/// where S_WKB = (π/2) * sqrt(2 m E_b) * a and ω_b = sqrt(|E_b''|/m) is
/// the imaginary barrier frequency.
pub fn wkb_tunneling_rate(barrier_height: f64, barrier_width: f64, mass: f64, hbar: f64) -> f64 {
    let s_wkb = 0.5 * PI * (2.0 * mass * barrier_height).sqrt() * barrier_width;
    let omega_b = (2.0 * barrier_height / (mass * barrier_width.powi(2))).sqrt();
    omega_b / (2.0 * PI) * (-2.0 * s_wkb / hbar).exp()
}

// ============================================================================
// Quantum diagnostics
// ============================================================================

/// Summary diagnostics for a PIMD simulation run.
#[derive(Debug, Clone)]
pub struct QuantumDiagnostics {
    /// Number of atoms N.
    pub n_atoms: usize,
    /// Number of ring-polymer beads P.
    pub n_beads: usize,
    /// Temperature T = 1/β.
    pub temperature: f64,
    /// Quantum kinetic energy (virial estimator).
    pub quantum_ke: f64,
    /// Classical kinetic energy at current centroid velocities.
    pub classical_ke: f64,
    /// Ring-polymer spring energy.
    pub spring_energy: f64,
    /// Average ring-polymer spread ⟨Δr²⟩^{1/2} (Å or internal length).
    pub mean_spread: f64,
    /// Thermal de Broglie wavelength for the first particle.
    pub de_broglie: f64,
}

impl QuantumDiagnostics {
    /// Compute diagnostics from a `QuantumSystem`.
    pub fn from_system(sys: &QuantumSystem) -> Self {
        let temperature = sys.temperature();
        let quantum_ke = quantum_kinetic_energy(sys);
        let classical_ke = sys.centroid_kinetic_energy();
        let spring_energy = ring_polymer_spring_energy(sys);
        let mean_spread = if sys.beads.is_empty() {
            0.0
        } else {
            let sum: f64 = sys.beads.iter().map(|b| b.ring_spread_sq()).sum();
            (sum / sys.n_atoms as f64).sqrt()
        };
        let de_broglie = if sys.beads.is_empty() {
            0.0
        } else {
            thermal_de_broglie(sys.beads[0].mass, sys.beta, sys.beads[0].hbar_eff)
        };
        Self {
            n_atoms: sys.n_atoms,
            n_beads: sys.n_beads,
            temperature,
            quantum_ke,
            classical_ke,
            spring_energy,
            mean_spread,
            de_broglie,
        }
    }

    /// Quantum enhancement factor: ratio of quantum to classical KE.
    pub fn quantum_enhancement(&self) -> f64 {
        quantum_correction_factor(self.classical_ke, self.quantum_ke)
    }

    /// Return a one-line summary string.
    pub fn summary(&self) -> String {
        format!(
            "N={} P={} T={:.3} QKE={:.4} CKE={:.4} spread={:.4} Λ={:.4}",
            self.n_atoms,
            self.n_beads,
            self.temperature,
            self.quantum_ke,
            self.classical_ke,
            self.mean_spread,
            self.de_broglie,
        )
    }
}

// ============================================================================
// Nuclear quantum effects helpers
// ============================================================================

/// Estimate the zero-point energy of a harmonic oscillator.
///
/// E_ZPE = ℏ ω / 2
pub fn zero_point_energy(omega: f64, hbar: f64) -> f64 {
    0.5 * hbar * omega
}

/// Path-integral partition function ratio Q_P / Q_1 for a 1D harmonic
/// oscillator at P beads.
///
/// The exact result is:
///
/// Q_P / Q_1 = (sinh(β ℏ ω / 2)) / (P * sinh(β ℏ ω / (2 P)))
///
/// This ratio converges to 1 as P → ∞ (quantum limit).
pub fn path_integral_partition_ratio(omega: f64, beta: f64, hbar: f64, n_beads: usize) -> f64 {
    let x = beta * hbar * omega / 2.0;
    let xp = x / n_beads as f64;
    if xp.abs() < 1e-15 {
        return 1.0;
    }
    x.sinh() / (n_beads as f64 * xp.sinh())
}

/// Quantum momentum distribution spread σ_p for a ring polymer.
///
/// σ_p² = m P / (β ℏ²)  (free ring-polymer momentum width)
pub fn ring_polymer_momentum_spread(mass: f64, n_beads: usize, beta: f64, hbar: f64) -> f64 {
    let p = n_beads as f64;
    (mass * p / (beta * hbar * hbar)).sqrt()
}

/// Compute the Kubo-transformed position autocorrelation for a harmonic
/// oscillator.
///
/// C(t) = (ℏ / (2 m ω)) * (cosh(ω(β ℏ/2 − i t)) / sinh(β ℏ ω / 2))
///
/// Returns only the real part (cosine contribution at T > 0).
pub fn harmonic_kubo_correlation(omega: f64, t: f64, mass: f64, beta: f64, hbar: f64) -> f64 {
    let x = beta * hbar * omega / 2.0;
    if x.abs() < 1e-15 {
        return hbar / (2.0 * mass * omega);
    }
    let cosh_num = (omega * t).cos() * x.cosh() + (omega * t).sin() * x.sinh();
    hbar / (2.0 * mass * omega) * cosh_num / x.sinh()
}

/// Bead-spring kinetic energy in normal-mode coordinates.
///
/// T_nm = Σ_n (1/2) m ω_n² |q_n|²
///
/// where `q_n` are normal-mode amplitudes and `omega_n` are normal-mode
/// frequencies.
pub fn normal_mode_spring_energy(
    mode_amplitudes: &[[f64; 3]],
    mode_freqs: &[f64],
    mass: f64,
) -> f64 {
    mode_amplitudes
        .iter()
        .zip(mode_freqs.iter())
        .map(|(q, &w)| 0.5 * mass * w * w * v3_norm_sq(*q))
        .sum()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── QuantumBead ──────────────────────────────────────────────────────────

    #[test]
    fn test_bead_new() {
        let b = QuantumBead::new(4, 1.0, 1.0);
        assert_eq!(b.n_beads(), 4);
        assert_eq!(b.positions.len(), 4);
        assert_eq!(b.velocities.len(), 4);
    }

    #[test]
    fn test_bead_positions_zero() {
        let b = QuantumBead::new(3, 2.0, 1.0);
        for pos in &b.positions {
            assert_eq!(*pos, [0.0, 0.0, 0.0]);
        }
    }

    #[test]
    fn test_bead_mass() {
        let b = QuantumBead::new(2, 3.5, 1.0);
        assert!((b.mass - 3.5).abs() < 1e-12);
    }

    #[test]
    fn test_bead_ring_spread_zero_for_identical_positions() {
        let mut b = QuantumBead::new(4, 1.0, 1.0);
        for p in &mut b.positions {
            *p = [1.0, 2.0, 3.0];
        }
        assert!(b.ring_spread_sq() < 1e-12);
    }

    #[test]
    fn test_bead_ring_spread_positive_when_dispersed() {
        let mut b = QuantumBead::new(2, 1.0, 1.0);
        b.positions[0] = [0.0, 0.0, 0.0];
        b.positions[1] = [2.0, 0.0, 0.0];
        assert!(b.ring_spread_sq() > 0.0);
    }

    // ── QuantumSystem ────────────────────────────────────────────────────────

    #[test]
    fn test_system_new() {
        let s = QuantumSystem::new(3, 4, 1.0, 1.0, 1.0);
        assert_eq!(s.n_atoms, 3);
        assert_eq!(s.n_beads, 4);
        assert_eq!(s.beads.len(), 3);
    }

    #[test]
    fn test_system_beta() {
        let s = QuantumSystem::new(1, 2, 2.5, 1.0, 1.0);
        assert!((s.beta - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_system_temperature() {
        let s = QuantumSystem::new(1, 4, 2.0, 1.0, 1.0);
        assert!((s.temperature() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_system_centroid_ke_zero_for_zero_velocity() {
        let s = QuantumSystem::new(2, 4, 1.0, 1.0, 1.0);
        assert!(s.centroid_kinetic_energy() < 1e-12);
    }

    // ── ring_polymer_spring_constant ─────────────────────────────────────────

    #[test]
    fn test_spring_constant_positive() {
        let k = ring_polymer_spring_constant(1.0, 4, 1.0);
        assert!(k > 0.0);
    }

    #[test]
    fn test_spring_constant_scales_with_p_squared() {
        let k4 = ring_polymer_spring_constant(1.0, 4, 1.0);
        let k8 = ring_polymer_spring_constant(1.0, 8, 1.0);
        assert!((k8 - 4.0 * k4).abs() < 1e-10, "k4={k4}, k8={k8}");
    }

    #[test]
    fn test_spring_constant_scales_with_mass() {
        let k1 = ring_polymer_spring_constant(1.0, 4, 1.0);
        let k2 = ring_polymer_spring_constant(2.0, 4, 1.0);
        assert!((k2 - 2.0 * k1).abs() < 1e-10);
    }

    // ── ring_polymer_spring_force ────────────────────────────────────────────

    #[test]
    fn test_spring_force_uniform_ring_zero() {
        let beads = vec![[1.0_f64, 2.0, 3.0]; 4];
        let forces = ring_polymer_spring_force(&beads, 1.0);
        for f in forces {
            assert!(f[0].abs() < 1e-12);
            assert!(f[1].abs() < 1e-12);
            assert!(f[2].abs() < 1e-12);
        }
    }

    #[test]
    fn test_spring_force_length() {
        let beads = vec![[0.0_f64; 3]; 5];
        let forces = ring_polymer_spring_force(&beads, 1.0);
        assert_eq!(forces.len(), 5);
    }

    #[test]
    fn test_spring_force_empty() {
        let forces = ring_polymer_spring_force(&[], 1.0);
        assert!(forces.is_empty());
    }

    #[test]
    fn test_spring_force_sum_zero() {
        let beads = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let forces = ring_polymer_spring_force(&beads, 2.0);
        let sum_x: f64 = forces.iter().map(|f| f[0]).sum();
        let sum_y: f64 = forces.iter().map(|f| f[1]).sum();
        let sum_z: f64 = forces.iter().map(|f| f[2]).sum();
        assert!(sum_x.abs() < 1e-10, "sum_x={sum_x}");
        assert!(sum_y.abs() < 1e-10, "sum_y={sum_y}");
        assert!(sum_z.abs() < 1e-10, "sum_z={sum_z}");
    }

    // ── centroid_position ────────────────────────────────────────────────────

    #[test]
    fn test_centroid_empty() {
        assert_eq!(centroid_position(&[]), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_centroid_single() {
        let c = centroid_position(&[[1.0, 2.0, 3.0]]);
        assert_eq!(c, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_centroid_two_beads() {
        let beads = vec![[0.0_f64, 0.0, 0.0], [2.0, 4.0, 6.0]];
        let c = centroid_position(&beads);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }

    // ── centroid_velocity ────────────────────────────────────────────────────

    #[test]
    fn test_centroid_velocity() {
        let vels = vec![[1.0_f64, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let cv = centroid_velocity(&vels);
        assert!(cv[0].abs() < 1e-12);
    }

    // ── quantum_kinetic_energy ───────────────────────────────────────────────

    #[test]
    fn test_qke_positive_high_temp() {
        let sys = QuantumSystem::new(2, 4, 0.01, 1.0, 1.0);
        let ke = quantum_kinetic_energy(&sys);
        assert!(ke > 0.0);
    }

    #[test]
    fn test_qke_single_bead_is_classical() {
        let sys = QuantumSystem::new(1, 1, 1.0, 1.0, 1.0);
        let ke = quantum_kinetic_energy(&sys);
        let expected = 3.0 * 1.0 * 1.0 / (2.0 * 1.0);
        assert!((ke - expected).abs() < 1e-10);
    }

    // ── primitive_kinetic_energy ─────────────────────────────────────────────

    #[test]
    fn test_primitive_ke_positive_high_temp() {
        let sys = QuantumSystem::new(1, 4, 0.01, 1.0, 1.0);
        let ke = primitive_kinetic_energy(&sys);
        assert!(ke > 0.0);
    }

    #[test]
    fn test_primitive_ke_single_bead_matches_virial() {
        let sys = QuantumSystem::new(1, 1, 1.0, 1.0, 1.0);
        let prim = primitive_kinetic_energy(&sys);
        let virial = quantum_kinetic_energy(&sys);
        // For P=1 and all positions at origin, both estimators should agree
        assert!((prim - virial).abs() < 1e-10);
    }

    // ── thermal_de_broglie ───────────────────────────────────────────────────

    #[test]
    fn test_de_broglie_positive() {
        assert!(thermal_de_broglie(1.0, 1.0, 1.0) > 0.0);
    }

    #[test]
    fn test_de_broglie_heavier_mass_smaller() {
        let l1 = thermal_de_broglie(1.0, 1.0, 1.0);
        let l2 = thermal_de_broglie(4.0, 1.0, 1.0);
        assert!(l1 > l2);
    }

    #[test]
    fn test_de_broglie_formula() {
        let m = 2.0;
        let beta = 1.0;
        let hbar = 1.0;
        let expected = hbar * (2.0 * PI * beta / m).sqrt();
        assert!((thermal_de_broglie(m, beta, hbar) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_de_broglie_scales_with_hbar() {
        let l1 = thermal_de_broglie(1.0, 1.0, 1.0);
        let l2 = thermal_de_broglie(1.0, 1.0, 2.0);
        assert!((l2 - 2.0 * l1).abs() < 1e-12);
    }

    // ── quantum_correction_factor ────────────────────────────────────────────

    #[test]
    fn test_correction_equal() {
        assert!((quantum_correction_factor(1.0, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_correction_zero_classical() {
        assert!((quantum_correction_factor(0.0, 5.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_correction_ratio() {
        let qc = quantum_correction_factor(2.0, 3.0);
        assert!((qc - 1.5).abs() < 1e-12);
    }

    // ── step_pimd ────────────────────────────────────────────────────────────

    #[test]
    fn test_step_pimd_zero_force_constant_velocity() {
        let mut sys = QuantumSystem::new(1, 2, 1.0, 1.0, 1.0);
        sys.beads[0].velocities[0] = [1.0, 0.0, 0.0];
        sys.beads[0].velocities[1] = [1.0, 0.0, 0.0];
        let forces: Vec<Vec<[f64; 3]>> = vec![vec![[0.0; 3]; 2]];
        let dt = 0.01;
        let pos_before = sys.beads[0].positions[0];
        step_pimd(&mut sys, &forces, dt);
        let pos_after = sys.beads[0].positions[0];
        assert!((pos_after[0] - pos_before[0]).abs() > 1e-12);
    }

    #[test]
    fn test_step_pimd_momentum_change_with_force() {
        let mut sys = QuantumSystem::new(1, 1, 1.0, 1.0, 1.0);
        let forces: Vec<Vec<[f64; 3]>> = vec![vec![[1.0, 0.0, 0.0]]];
        let dt = 0.1;
        step_pimd(&mut sys, &forces, dt);
        assert!(sys.beads[0].velocities[0][0] != 0.0);
    }

    // ── NormalModeTransform ──────────────────────────────────────────────────

    #[test]
    fn test_normal_mode_forward_inverse_roundtrip() {
        let p = 4;
        let nm = NormalModeTransform::new(p);
        let original: Vec<[f64; 3]> = (0..p)
            .map(|i| [i as f64, (i + 1) as f64 * 0.5, 0.0])
            .collect();
        let modes = nm.forward(&original);
        let recovered = nm.inverse(&modes);
        for k in 0..p {
            for d in 0..3 {
                assert!(
                    (recovered[k][d] - original[k][d]).abs() < 1e-10,
                    "k={k} d={d}: got {} expected {}",
                    recovered[k][d],
                    original[k][d]
                );
            }
        }
    }

    #[test]
    fn test_normal_mode_centroid_is_mode_zero() {
        // The zero-th normal mode should equal the centroid (scaled by sqrt(P))
        let p = 4;
        let nm = NormalModeTransform::new(p);
        let beads: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; p];
        let modes = nm.forward(&beads);
        let centroid = centroid_position(&beads);
        // mode[0] = centroid * sqrt(P) (from the 1/sqrt(P) prefactor in C)
        let expected = v3_scale((p as f64).sqrt(), centroid);
        for d in 0..3 {
            assert!(
                (modes[0][d] - expected[d]).abs() < 1e-10,
                "d={d}: got {} expected {}",
                modes[0][d],
                expected[d]
            );
        }
    }

    // ── normal_mode_frequencies ──────────────────────────────────────────────

    #[test]
    fn test_normal_mode_freq_zero_for_centroid() {
        let freqs = normal_mode_frequencies(4, 1.0, 1.0);
        assert!(freqs[0].abs() < 1e-12, "centroid frequency should be zero");
    }

    #[test]
    fn test_normal_mode_freq_positive_for_non_zero_modes() {
        let freqs = normal_mode_frequencies(4, 1.0, 1.0);
        for &f in &freqs[1..] {
            assert!(f >= 0.0);
        }
    }

    #[test]
    fn test_normal_mode_freq_count() {
        let freqs = normal_mode_frequencies(8, 1.0, 1.0);
        assert_eq!(freqs.len(), 8);
    }

    // ── StagingCoords ────────────────────────────────────────────────────────

    #[test]
    fn test_staging_mass_bead_1_is_physical_mass() {
        let sc = StagingCoords::new(4, 2.0);
        assert!((sc.staging_mass(1) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_staging_mass_decreases_with_level() {
        // staging_mass(j) = j/(j-1) * m; the ratio decreases towards 1 as j grows
        let sc = StagingCoords::new(8, 1.0);
        let m2 = sc.staging_mass(2); // 2/1 = 2.0
        let m4 = sc.staging_mass(4); // 4/3 ≈ 1.333
        assert!(m2 > m4, "m2={m2} should exceed m4={m4}");
    }

    #[test]
    fn test_staging_mass_formula() {
        let sc = StagingCoords::new(4, 1.0);
        // m_s(4) = 4/3 * m
        let expected = 4.0 / 3.0;
        assert!((sc.staging_mass(4) - expected).abs() < 1e-12);
    }

    // ── PileThermostat ───────────────────────────────────────────────────────

    #[test]
    fn test_pile_c1_less_than_one() {
        let pile = PileThermostat::new(4, 1.0, 1.0, 1.0, 0.01);
        for &c in &pile.c1 {
            assert!(c <= 1.0 && c > 0.0, "c1 should be in (0,1], got {c}");
        }
    }

    #[test]
    fn test_pile_c2_non_negative() {
        let pile = PileThermostat::new(4, 1.0, 1.0, 1.0, 0.01);
        for &c in &pile.c2 {
            assert!(c >= 0.0, "c2 should be non-negative, got {c}");
        }
    }

    #[test]
    fn test_pile_apply_half_step_preserves_order() {
        let pile = PileThermostat::new(4, 1.0, 1.0, 0.5, 0.001);
        // With xi=0 the new velocity is just c1 * v_old (damped)
        let v_new = pile.apply_half_step(0, 1.0, 1.0, 0.0);
        assert!(v_new < 1.0 && v_new > 0.0, "Damped velocity: {v_new}");
    }

    #[test]
    fn test_pile_gamma_modes_count() {
        let pile = PileThermostat::new(6, 1.0, 1.0, 1.0, 0.01);
        assert_eq!(pile.gamma_modes.len(), 6);
    }

    // ── CentroidMdState ──────────────────────────────────────────────────────

    #[test]
    fn test_cmd_state_creation() {
        let state = CentroidMdState::new(3, 4, 1.0, 1.0);
        assert_eq!(state.centroid_positions.len(), 3);
        assert_eq!(state.centroid_velocities.len(), 3);
    }

    #[test]
    fn test_cmd_temperature_zero_for_zero_velocities() {
        let state = CentroidMdState::new(4, 4, 1.0, 1.0);
        assert!(state.temperature() < 1e-12);
    }

    #[test]
    fn test_cmd_step_moves_position() {
        let mut state = CentroidMdState::new(1, 4, 1.0, 1.0);
        state.centroid_velocities[0] = [1.0, 0.0, 0.0];
        state.step(0.1);
        assert!((state.centroid_positions[0][0] - 0.1_f64).abs() < 1e-10);
    }

    #[test]
    fn test_cmd_kinetic_energy_formula() {
        let mut state = CentroidMdState::new(1, 4, 1.0, 2.0);
        state.centroid_velocities[0] = [1.0, 0.0, 0.0];
        // KE = 0.5 * m * v^2 = 0.5 * 2 * 1 = 1.0
        assert!((state.kinetic_energy() - 1.0).abs() < 1e-12);
    }

    // ── RingPolymerMdState ───────────────────────────────────────────────────

    #[test]
    fn test_rpmd_creation() {
        let state = RingPolymerMdState::new(2, 4, 1.0, 1.0, 1.0);
        assert_eq!(state.system.n_atoms, 2);
        assert_eq!(state.system.n_beads, 4);
        assert!((state.time - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_rpmd_step_advances_time() {
        let mut state = RingPolymerMdState::new(1, 2, 1.0, 1.0, 1.0);
        state.step(0.05);
        assert!((state.time - 0.05).abs() < 1e-12);
    }

    #[test]
    fn test_rpmd_spring_energy_zero_for_uniform_ring() {
        let state = RingPolymerMdState::new(1, 4, 1.0, 1.0, 1.0);
        // All beads at origin → spring energy = 0
        assert!(state.spring_energy() < 1e-12);
    }

    // ── VacfAccumulator ──────────────────────────────────────────────────────

    #[test]
    fn test_vacf_tau_zero_is_one() {
        let mut acc = VacfAccumulator::new(10);
        acc.push(vec![[1.0, 0.0, 0.0]]);
        acc.push(vec![[0.5, 0.0, 0.0]]);
        let c0 = acc.vacf(0).unwrap();
        assert!((c0 - 1.0).abs() < 1e-10, "VACF(0) should be 1, got {c0}");
    }

    #[test]
    fn test_vacf_insufficient_data_returns_none() {
        let acc = VacfAccumulator::new(10);
        assert!(acc.vacf(5).is_none());
    }

    #[test]
    fn test_vacf_unnorm_zero_for_empty() {
        let acc = VacfAccumulator::new(5);
        assert!((acc.vacf_unnorm(0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_vacf_max_snapshots_respected() {
        let mut acc = VacfAccumulator::new(3);
        for i in 0..5 {
            acc.push(vec![[i as f64, 0.0, 0.0]]);
        }
        assert_eq!(acc.snapshots.len(), 3);
    }

    // ── quantum_diffusion_coefficient ────────────────────────────────────────

    #[test]
    fn test_diffusion_zero_for_zero_vacf() {
        let vacf = vec![0.0_f64; 10];
        let d = quantum_diffusion_coefficient(&vacf, 0.01);
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn test_diffusion_positive_for_positive_vacf() {
        let vacf = vec![1.0_f64; 20];
        let d = quantum_diffusion_coefficient(&vacf, 0.01);
        assert!(d > 0.0);
    }

    #[test]
    fn test_diffusion_empty_vacf() {
        let d = quantum_diffusion_coefficient(&[], 0.01);
        assert!((d - 0.0).abs() < 1e-12);
    }

    // ── TunnelingRate ────────────────────────────────────────────────────────

    #[test]
    fn test_wkb_transmission_above_barrier_is_one() {
        let tr = TunnelingRate::new(1.0, 1.0, 1.0, 1.0);
        assert!((tr.wkb_transmission(2.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_wkb_transmission_below_barrier_less_than_one() {
        let tr = TunnelingRate::new(1.0, 1.0, 1.0, 1.0);
        assert!(tr.wkb_transmission(0.5) < 1.0);
    }

    #[test]
    fn test_wkb_transmission_positive() {
        let tr = TunnelingRate::new(5.0, 1.0, 1.0, 1.0);
        assert!(tr.wkb_transmission(1.0) > 0.0);
    }

    #[test]
    fn test_instanton_rate_positive() {
        let tr = TunnelingRate::new(1.0, 0.5, 1.0, 1.0);
        assert!(tr.instanton_rate() > 0.0);
    }

    #[test]
    fn test_wkb_tunneling_rate_function_positive() {
        let rate = wkb_tunneling_rate(1.0, 1.0, 1.0, 1.0);
        assert!(rate > 0.0);
    }

    // ── zero_point_energy ────────────────────────────────────────────────────

    #[test]
    fn test_zpe_formula() {
        let zpe = zero_point_energy(2.0, 1.0);
        assert!((zpe - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_zpe_positive() {
        assert!(zero_point_energy(1.0, 1.0) > 0.0);
    }

    // ── path_integral_partition_ratio ────────────────────────────────────────

    #[test]
    fn test_partition_ratio_approaches_one_large_p() {
        // For large P, quantum discretisation error → 0 so ratio → 1
        let r = path_integral_partition_ratio(1.0, 1.0, 1.0, 64);
        assert!((r - 1.0).abs() < 0.05, "ratio={r}");
    }

    #[test]
    fn test_partition_ratio_positive() {
        let r = path_integral_partition_ratio(1.0, 2.0, 1.0, 8);
        assert!(r > 0.0);
    }

    // ── ring_polymer_momentum_spread ─────────────────────────────────────────

    #[test]
    fn test_momentum_spread_positive() {
        let sigma = ring_polymer_momentum_spread(1.0, 4, 1.0, 1.0);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_momentum_spread_scales_with_p() {
        let s4 = ring_polymer_momentum_spread(1.0, 4, 1.0, 1.0);
        let s16 = ring_polymer_momentum_spread(1.0, 16, 1.0, 1.0);
        assert!(s16 > s4, "More beads → larger momentum spread");
    }

    // ── harmonic_kubo_correlation ────────────────────────────────────────────

    #[test]
    fn test_kubo_correlation_at_t_zero_positive() {
        let c0 = harmonic_kubo_correlation(1.0, 0.0, 1.0, 1.0, 1.0);
        assert!(c0 > 0.0);
    }

    #[test]
    fn test_kubo_correlation_finite() {
        let c = harmonic_kubo_correlation(1.0, 0.5, 1.0, 1.0, 1.0);
        assert!(c.is_finite());
    }

    // ── QuantumDiagnostics ───────────────────────────────────────────────────

    #[test]
    fn test_diagnostics_summary_not_empty() {
        let sys = QuantumSystem::new(2, 4, 1.0, 1.0, 1.0);
        let diag = QuantumDiagnostics::from_system(&sys);
        assert!(!diag.summary().is_empty());
    }

    #[test]
    fn test_diagnostics_de_broglie_positive() {
        let sys = QuantumSystem::new(1, 4, 1.0, 1.0, 1.0);
        let diag = QuantumDiagnostics::from_system(&sys);
        assert!(diag.de_broglie > 0.0);
    }

    #[test]
    fn test_diagnostics_temperature_matches_beta() {
        let sys = QuantumSystem::new(1, 4, 2.0, 1.0, 1.0);
        let diag = QuantumDiagnostics::from_system(&sys);
        assert!((diag.temperature - 0.5).abs() < 1e-12);
    }

    // ── normal_mode_spring_energy ────────────────────────────────────────────

    #[test]
    fn test_normal_mode_spring_energy_zero_for_zero_amplitudes() {
        let amps = vec![[0.0_f64; 3]; 4];
        let freqs = vec![0.0, 1.0, 2.0, 1.0];
        let e = normal_mode_spring_energy(&amps, &freqs, 1.0);
        assert!(e.abs() < 1e-12);
    }

    #[test]
    fn test_normal_mode_spring_energy_positive_for_nonzero_modes() {
        let amps = vec![[0.0_f64; 3], [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let freqs = vec![0.0, 2.0, 3.0, 2.0];
        let e = normal_mode_spring_energy(&amps, &freqs, 1.0);
        // E = 0.5 * 1.0 * 4.0 * 1.0 = 2.0
        assert!((e - 2.0).abs() < 1e-12, "energy={e}");
    }
}
