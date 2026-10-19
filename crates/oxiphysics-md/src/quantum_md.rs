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
//! # Conventions
//!
//! A particle of mass `m` and effective Planck constant `ℏ` is a ring of `P`
//! beads obeying the ring-polymer Hamiltonian (Craig-Manolopoulos RPMD, i-PI)
//!
//! ```text
//! H_P = Σ_p [ m |v_p|² / 2 + ½ k |r_{p+1} − r_p|² + V(r_p) ],
//! k = m ω_P²,   ω_P = P / (β ℏ),   r_P ≡ r_0,
//! ```
//!
//! i.e. the *full* physical potential `V` acts on every bead and the external
//! force of bead `p` is `forces[i][p] = −∇V(r_ip)` (never the spring force).
//! The classical Boltzmann distribution of `H_P` at the inverse temperature
//! `β / P` reproduces the quantum canonical ensemble at `β`; the free ring
//! polymer decouples into the normal modes of [`NormalModeTransform`] with the
//! frequencies of [`normal_mode_frequencies`].  Units are those of the caller
//! with `k_B = 1`.
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
//! - [`ring_polymer_spring_constant`]: spring constant `m (P / (β ℏ))²` for the
//!   inter-bead harmonic coupling.
//! - [`ring_polymer_spring_force`]: force from harmonic springs along the ring.
//! - [`centroid_position`] / [`centroid_velocity`]: classical centroid
//!   estimators.
//! - [`quantum_kinetic_energy`]: centroid-virial estimator for the quantum
//!   kinetic energy (needs the external forces on every bead).
//! - [`primitive_kinetic_energy`]: primitive (thermodynamic) kinetic energy
//!   estimator.
//! - [`thermal_de_broglie`]: thermal de Broglie wavelength.
//! - [`step_pimd`]: one velocity-Verlet step of the ring-polymer MD: half-kick,
//!   exact normal-mode propagation of the free ring polymer, force refresh,
//!   half-kick.
//! - [`step_rpmd`]: the same integrator under its RPMD name.
//! - [`normal_mode_frequencies`]: compute ring-polymer normal-mode frequencies.
//! - [`wkb_tunneling_rate`]: WKB tunneling rate through a parabolic barrier.
//! - [`quantum_diffusion_coefficient`]: Green-Kubo integral of the VACF.

use crate::MdError;
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
    /// Normal-mode transform and mode buffers reused by [`step_pimd`]; rebuilt
    /// whenever the number of beads changes.
    mode_scratch: Option<ModeScratch>,
}

/// Normal-mode transform and the mode-space buffers of one ring polymer,
/// cached in [`QuantumSystem`] so that [`step_pimd`] does not rebuild the
/// `P × P` transform matrix or reallocate on every step.
#[derive(Debug, Clone)]
struct ModeScratch {
    transform: NormalModeTransform,
    position_modes: Vec<[f64; 3]>,
    velocity_modes: Vec<[f64; 3]>,
}

impl ModeScratch {
    fn new(n_beads: usize) -> Self {
        Self {
            transform: NormalModeTransform::new(n_beads),
            position_modes: vec![[0.0; 3]; n_beads],
            velocity_modes: vec![[0.0; 3]; n_beads],
        }
    }
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
            mode_scratch: None,
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
        let mut modes = vec![[0.0f64; 3]; self.n_beads];
        self.forward_into(bead_coords, &mut modes);
        modes
    }

    /// Forward transform into a caller-provided buffer (no allocation).
    ///
    /// `modes[n] = Σ_k C[k][n] · bead_coords[k]`; at most `P` entries of both
    /// slices are used.
    pub fn forward_into(&self, bead_coords: &[[f64; 3]], modes: &mut [[f64; 3]]) {
        let p = self.n_beads;
        for (n, mode_n) in modes.iter_mut().enumerate().take(p) {
            let mut amplitude = [0.0f64; 3];
            for (k, &bc_k) in bead_coords.iter().enumerate().take(p) {
                let c_kn = self.c_matrix[k * p + n];
                amplitude = v3_add(amplitude, v3_scale(c_kn, bc_k));
            }
            *mode_n = amplitude;
        }
    }

    /// Apply inverse transform: convert normal-mode amplitudes back to bead
    /// positions.
    pub fn inverse(&self, mode_coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut beads = vec![[0.0f64; 3]; self.n_beads];
        self.inverse_into(mode_coords, &mut beads);
        beads
    }

    /// Inverse transform into a caller-provided buffer (no allocation).
    ///
    /// `beads[k] = Σ_n C[k][n] · mode_coords[n]`; at most `P` entries of both
    /// slices are used.
    pub fn inverse_into(&self, mode_coords: &[[f64; 3]], beads: &mut [[f64; 3]]) {
        let p = self.n_beads;
        for (k, bead_k) in beads.iter_mut().enumerate().take(p) {
            let mut position = [0.0f64; 3];
            for (n, &mc_n) in mode_coords.iter().enumerate().take(p) {
                let c_kn = self.c_matrix[k * p + n];
                position = v3_add(position, v3_scale(c_kn, mc_n));
            }
            *bead_k = position;
        }
    }
}

/// Compute the ring-polymer normal-mode angular frequencies ω_n for P beads
/// at inverse temperature β.
///
/// ω_n = (2 P / (β ℏ)) sin(π n / P),   n = 0, 1, …, P-1.
///
/// Mode n=0 is the centroid (ω=0).  The index `n` is the column index of
/// [`NormalModeTransform`]; the sine modes (`n > P/2`) carry the frequency of
/// their cosine partner `P − n`, which `sin(π n / P)` reproduces.
pub fn normal_mode_frequencies(n_beads: usize, beta: f64, hbar: f64) -> Vec<f64> {
    (0..n_beads)
        .map(|n| mode_frequency(n, n_beads, beta, hbar))
        .collect()
}

/// Angular frequency `ω_n = (2P / (β ℏ)) sin(π n / P)` of normal mode `n`.
#[inline]
fn mode_frequency(n: usize, n_beads: usize, beta: f64, hbar: f64) -> f64 {
    let p = n_beads as f64;
    2.0 * p / (beta * hbar) * (PI * (n as f64) / p).sin()
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
/// reduced Planck constant of the particle (`QuantumBead::hbar_eff`).  This is
/// `k = m ω_P²` with `ω_P = P / (β ℏ)`; the springs act between neighbouring
/// beads with the *full* weight of the physical potential on each bead (see the
/// module-level conventions).
pub fn ring_polymer_spring_constant(mass: f64, n_beads: usize, beta: f64, hbar: f64) -> f64 {
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

/// Check that `sys` is internally consistent and that `forces` has one entry per
/// bead of every ring polymer.
fn check_bead_forces(sys: &QuantumSystem, forces: &[Vec<[f64; 3]>]) -> Result<(), MdError> {
    if sys.beads.len() != sys.n_atoms {
        return Err(MdError::General(format!(
            "QuantumSystem declares {} atoms but holds {} ring polymers",
            sys.n_atoms,
            sys.beads.len()
        )));
    }
    if forces.len() != sys.beads.len() {
        return Err(MdError::General(format!(
            "bead forces given for {} atoms but the system has {}",
            forces.len(),
            sys.beads.len()
        )));
    }
    for (i, (bead, bead_forces)) in sys.beads.iter().zip(forces.iter()).enumerate() {
        if bead.positions.len() != sys.n_beads || bead.velocities.len() != sys.n_beads {
            return Err(MdError::General(format!(
                "atom {i} has {} positions and {} velocities but the system has {} beads",
                bead.positions.len(),
                bead.velocities.len(),
                sys.n_beads
            )));
        }
        if bead_forces.len() != sys.n_beads {
            return Err(MdError::General(format!(
                "atom {i} has {} bead forces but the system has {} beads",
                bead_forces.len(),
                sys.n_beads
            )));
        }
    }
    Ok(())
}

/// Quantum kinetic energy using the centroid-virial estimator.
///
/// T_cv = 3N/(2β) − (1/(2P)) Σ_{i,p} (r_ip − c_i) · F_ip
///
/// where `c_i` is the centroid of atom `i` and `forces[i][p] = −∇V(r_ip)` is the
/// external (physical) force on bead `p` of atom `i`, exactly the array that is
/// passed to [`step_pimd`].  For a harmonic potential `F = −mω²r` the correction
/// is non-negative, so `T_cv ≥ 3N/(2β)`.  Unlike the primitive estimator its
/// variance stays bounded as the number of beads grows.
///
/// For one bead the ring has no internal structure and the estimator reduces to
/// the classical `3N/(2β)`.
///
/// # Errors
/// [`MdError::General`] if `forces` does not hold one force per bead of every
/// ring polymer.
pub fn quantum_kinetic_energy(
    sys: &QuantumSystem,
    forces: &[Vec<[f64; 3]>],
) -> Result<f64, MdError> {
    check_bead_forces(sys, forces)?;
    if sys.n_beads == 0 {
        return Ok(0.0);
    }
    let n = sys.n_atoms as f64;
    let p = sys.n_beads as f64;
    let classical_term = 3.0 * n / (2.0 * sys.beta);
    let mut virial_sum = 0.0;
    for (bead, bead_forces) in sys.beads.iter().zip(forces.iter()) {
        let c = centroid_position(&bead.positions);
        for (&pos, &force) in bead.positions.iter().zip(bead_forces.iter()) {
            virial_sum += v3_dot(v3_sub(pos, c), force);
        }
    }
    Ok(classical_term - virial_sum / (2.0 * p))
}

/// Primitive (thermodynamic) kinetic energy estimator.
///
/// T_prim = (3 N P) / (2β) − Σ_{i} (k_i / (2P)) Σ_p |r_{i,p+1} − r_{i,p}|²
///
/// with the spring constant `k_i = m_i (P / (β ℏ_i))²` of
/// [`ring_polymer_spring_constant`] evaluated with the particle's own
/// `hbar_eff`, i.e. `k_i / (2P) = m_i P / (2 β² ℏ_i²)` per link.  The
/// estimator needs no forces but has a variance that grows with the number of
/// beads (compare [`quantum_kinetic_energy`]).
pub fn primitive_kinetic_energy(sys: &QuantumSystem) -> f64 {
    let n = sys.n_atoms as f64;
    let p_count = sys.n_beads;
    if p_count == 0 {
        return 0.0;
    }
    let p = p_count as f64;
    let classical_term = 3.0 * n * p / (2.0 * sys.beta);
    let mut spring_term = 0.0;
    for bead in &sys.beads {
        let k = ring_polymer_spring_constant(bead.mass, p_count, sys.beta, bead.hbar_eff);
        for j in 0..p_count {
            let rj = bead.positions[j];
            let rj1 = bead.positions[(j + 1) % p_count];
            let dr = v3_sub(rj1, rj);
            spring_term += k / (2.0 * p) * v3_dot(dr, dr);
        }
    }
    classical_term - spring_term
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
        let k = ring_polymer_spring_constant(bead.mass, sys.n_beads, sys.beta, bead.hbar_eff);
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
///
/// The external forces are cached between steps (velocity Verlet): they must be
/// valid for the current bead positions before the first call of
/// [`RingPolymerMdState::step`], which [`RingPolymerMdState::refresh_forces`]
/// takes care of; `step` itself refreshes them after every drift.
#[derive(Debug, Clone)]
pub struct RingPolymerMdState {
    /// Full ring-polymer system (P beads per atom).
    pub system: QuantumSystem,
    /// External forces `−∇V` on each bead for the *current* positions:
    /// `forces[atom][bead] = [fx, fy, fz]`.
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

    /// Recompute the cached external forces for the current bead positions
    /// with `force_fn` (same contract as in [`step_pimd`]).  Call this once
    /// after setting the initial positions: [`RingPolymerMdState::new`] starts
    /// with zero forces.
    pub fn refresh_forces<F>(&mut self, mut force_fn: F)
    where
        F: FnMut(&QuantumSystem, &mut [Vec<[f64; 3]>]),
    {
        force_fn(&self.system, &mut self.forces);
    }

    /// Advance by one velocity-Verlet step including the ring-polymer springs.
    ///
    /// `force_fn` recomputes the external forces after the drift (see
    /// [`step_pimd`]); the cached `forces` are updated in place and the time
    /// advances by `dt`.
    ///
    /// # Errors
    /// [`MdError::General`] if the cached forces do not match the system shape.
    pub fn step<F>(&mut self, dt: f64, force_fn: F) -> Result<(), MdError>
    where
        F: FnMut(&QuantumSystem, &mut [Vec<[f64; 3]>]),
    {
        step_rpmd(&mut self.system, &mut self.forces, dt, force_fn)?;
        self.time += dt;
        Ok(())
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

/// Rotate one normal mode `(q, v)` of frequency `omega` exactly over `dt`:
/// `q' = q cos ωdt + (v/ω) sin ωdt`, `v' = −ω q sin ωdt + v cos ωdt`; the
/// centroid (`ω = 0`) drifts freely.
#[inline]
fn rotate_mode(q: &mut [f64; 3], v: &mut [f64; 3], omega: f64, dt: f64) {
    if omega == 0.0 {
        *q = v3_add(*q, v3_scale(dt, *v));
        return;
    }
    let (sin_wt, cos_wt) = (omega * dt).sin_cos();
    let q0 = *q;
    let v0 = *v;
    *q = v3_add(v3_scale(cos_wt, q0), v3_scale(sin_wt / omega, v0));
    *v = v3_add(v3_scale(-omega * sin_wt, q0), v3_scale(cos_wt, v0));
}

/// Advance every ring polymer of `sys` by `dt` under its springs alone.
///
/// The spring Hamiltonian is diagonal in the normal modes of
/// [`NormalModeTransform`] (all of mass `m`), so each mode is propagated by the
/// exact harmonic rotation of [`rotate_mode`] with `ω_n` from
/// [`normal_mode_frequencies`] and the particle's own `ℏ`.  The result is exact
/// for any `dt` (the integrator is not limited by the stiff spring frequencies)
/// and conserves the free ring-polymer energy to rounding error.
fn propagate_free_ring_polymers(sys: &mut QuantumSystem, dt: f64) {
    let n_beads = sys.n_beads;
    if n_beads == 0 {
        return;
    }
    let beta = sys.beta;
    let mut scratch = match sys.mode_scratch.take() {
        Some(cached) if cached.transform.n_beads == n_beads => cached,
        _ => ModeScratch::new(n_beads),
    };
    for bead in sys.beads.iter_mut() {
        scratch
            .transform
            .forward_into(&bead.positions, &mut scratch.position_modes);
        scratch
            .transform
            .forward_into(&bead.velocities, &mut scratch.velocity_modes);
        for (n, (q_mode, v_mode)) in scratch
            .position_modes
            .iter_mut()
            .zip(scratch.velocity_modes.iter_mut())
            .enumerate()
        {
            let omega = mode_frequency(n, n_beads, beta, bead.hbar_eff);
            rotate_mode(q_mode, v_mode, omega, dt);
        }
        scratch
            .transform
            .inverse_into(&scratch.position_modes, &mut bead.positions);
        scratch
            .transform
            .inverse_into(&scratch.velocity_modes, &mut bead.velocities);
    }
    sys.mode_scratch = Some(scratch);
}

/// Kick every bead velocity with `forces` over the time `half_dt`
/// (`v += half_dt F / m`).
fn kick_beads(sys: &mut QuantumSystem, forces: &[Vec<[f64; 3]>], half_dt: f64) {
    for (bead, bead_forces) in sys.beads.iter_mut().zip(forces.iter()) {
        let scale = half_dt / bead.mass;
        for (velocity, &force) in bead.velocities.iter_mut().zip(bead_forces.iter()) {
            *velocity = v3_add(*velocity, v3_scale(scale, force));
        }
    }
}

/// One velocity-Verlet step for the ring-polymer MD.
///
/// `forces[i][p]` is the external force `−∇V(r_ip)` on atom `i`, bead `p`, and
/// must be valid for the bead positions on entry (the previous step leaves it in
/// that state; [`RingPolymerMdState::refresh_forces`] initialises it).  The step
/// is the symmetric splitting
///
/// 1. half-kick `v += (dt / 2m) F` with the cached forces,
/// 2. exact propagation of the free ring polymer over `dt` in normal-mode space
///    (harmonic rotation of every mode with `ω_n`, free drift of the centroid),
/// 3. `force_fn(sys, forces)` recomputes the external forces at the **new**
///    positions and stores them in `forces`,
/// 4. half-kick with the refreshed forces.
///
/// The map is symplectic and time-reversible; it conserves
/// `H_P = Σ m v²/2 + Σ_p ½ k |r_{p+1} − r_p|² + Σ_p V(r_p)` up to a bounded
/// `O(dt²)` fluctuation and conserves the free ring-polymer energy to rounding
/// error for any `dt`.  The previous implementation kicked the second half with
/// the *old* external forces and only re-evaluated the springs, which is not
/// velocity Verlet for any position-dependent potential.
///
/// The normal-mode transform is cached in `sys` and rebuilt only when the number
/// of beads changes.
///
/// # Errors
/// [`MdError::General`] if `forces` does not hold one force per bead of every
/// ring polymer.  A wrongly shaped `forces` on entry leaves `sys` untouched; if
/// `force_fn` leaves it wrongly shaped, the beads have already been advanced and
/// the second half-kick is not applied.
pub fn step_pimd<F>(
    sys: &mut QuantumSystem,
    forces: &mut [Vec<[f64; 3]>],
    dt: f64,
    mut force_fn: F,
) -> Result<(), MdError>
where
    F: FnMut(&QuantumSystem, &mut [Vec<[f64; 3]>]),
{
    check_bead_forces(sys, forces)?;
    kick_beads(sys, forces, 0.5 * dt);
    propagate_free_ring_polymers(sys, dt);
    force_fn(sys, forces);
    check_bead_forces(sys, forces)?;
    kick_beads(sys, forces, 0.5 * dt);
    Ok(())
}

/// One velocity-Verlet step for RPMD (identical to PIMD integrator but named
/// separately for conceptual clarity and potential future divergence).
///
/// # Errors
/// See [`step_pimd`].
pub fn step_rpmd<F>(
    sys: &mut QuantumSystem,
    forces: &mut [Vec<[f64; 3]>],
    dt: f64,
    force_fn: F,
) -> Result<(), MdError>
where
    F: FnMut(&QuantumSystem, &mut [Vec<[f64; 3]>]),
{
    step_pimd(sys, forces, dt, force_fn)
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
    /// Compute diagnostics from a `QuantumSystem` and the external bead forces
    /// `forces[i][p] = −∇V(r_ip)` needed by the virial kinetic-energy estimator.
    ///
    /// # Errors
    /// [`MdError::General`] if `forces` does not match the shape of `sys`
    /// (see [`quantum_kinetic_energy`]).
    pub fn from_system(sys: &QuantumSystem, forces: &[Vec<[f64; 3]>]) -> Result<Self, MdError> {
        let temperature = sys.temperature();
        let quantum_ke = quantum_kinetic_energy(sys, forces)?;
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
        Ok(Self {
            n_atoms: sys.n_atoms,
            n_beads: sys.n_beads,
            temperature,
            quantum_ke,
            classical_ke,
            spring_energy,
            mean_spread,
            de_broglie,
        })
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
mod pimd_tests;

#[cfg(test)]
mod tests;
