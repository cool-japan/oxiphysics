// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polymer chain simulation module.
//!
//! This module provides:
//! - [`PolymerBead`]: Single bead (position, velocity, mass, charge).
//! - [`FeneSpring`]: Finitely Extensible Nonlinear Elastic (FENE) spring.
//! - [`KremerGrestChain`]: Full bead-spring chain (Kremer–Grest model).
//! - [`WormlikeChain`]: WLC force-extension via the Marko–Siggia formula.
//! - [`RadiusOfGyration`]: Rg computation from bead positions.
//! - [`RouseModel`]: Normal-mode relaxation times and diffusion coefficient.
//! - [`flory_exponent_good_solvent`]: ν ≈ 0.588 (Flory exponent, good solvent).
//! - [`end_to_end_ideal`]: ⟨R²⟩^½ for an ideal chain.
//! - [`rg_ideal`]: Radius of gyration for an ideal (Gaussian) chain.

// ---------------------------------------------------------------------------
// Helper: 3-D vector arithmetic
// ---------------------------------------------------------------------------

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn mag3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Flory exponent ν for a polymer in a good solvent (≈ 0.588).
///
/// Gives the scaling of the end-to-end distance: R ~ N^ν.
///
/// ```no_run
/// use oxiphysics_md::polymer_sim::flory_exponent_good_solvent;
/// assert!((flory_exponent_good_solvent() - 0.588).abs() < 1e-12);
/// ```
pub fn flory_exponent_good_solvent() -> f64 {
    0.588
}

/// Root-mean-square end-to-end distance of an ideal (Gaussian) chain.
///
/// `R = sqrt(N) * b`
///
/// # Arguments
/// * `n` – Number of Kuhn segments.
/// * `b` – Kuhn segment length (m).
///
/// ```no_run
/// use oxiphysics_md::polymer_sim::end_to_end_ideal;
/// let r = end_to_end_ideal(100, 1.0);
/// assert!((r - 10.0).abs() < 1e-10);
/// ```
pub fn end_to_end_ideal(n: usize, b: f64) -> f64 {
    (n as f64).sqrt() * b
}

/// Radius of gyration of an ideal (Gaussian) chain.
///
/// `R_g = b * sqrt(N / 6)`
///
/// # Arguments
/// * `n` – Number of Kuhn segments.
/// * `b` – Kuhn segment length (m).
///
/// ```no_run
/// use oxiphysics_md::polymer_sim::rg_ideal;
/// let rg = rg_ideal(6, 1.0);
/// assert!((rg - 1.0).abs() < 1e-10);
/// ```
pub fn rg_ideal(n: usize, b: f64) -> f64 {
    b * (n as f64 / 6.0).sqrt()
}

// ---------------------------------------------------------------------------
// PolymerBead
// ---------------------------------------------------------------------------

/// A single bead in a bead-spring polymer model.
///
/// Stores Newtonian degrees of freedom plus physical parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct PolymerBead {
    /// Position vector (m).
    pub position: [f64; 3],
    /// Velocity vector (m s⁻¹).
    pub velocity: [f64; 3],
    /// Bead mass (kg or reduced units).
    pub mass: f64,
    /// Partial charge (e or reduced units).
    pub charge: f64,
}

impl PolymerBead {
    /// Create a new polymer bead at rest.
    pub fn new(position: [f64; 3], mass: f64, charge: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            charge,
        }
    }

    /// Kinetic energy ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        0.5 * self.mass * v2
    }
}

// ---------------------------------------------------------------------------
// FeneSpring
// ---------------------------------------------------------------------------

/// Finitely Extensible Nonlinear Elastic (FENE) spring.
///
/// The FENE potential diverges as r → R_max, preventing chain crossing.
/// Force magnitude: `F(r) = -k r / (1 - (r/R_max)²)`.
#[derive(Debug, Clone)]
pub struct FeneSpring {
    /// Maximum extension R_max (reduced units or m).
    pub r_max: f64,
    /// Spring constant k (reduced units or N m⁻¹).
    pub k_spring: f64,
}

impl FeneSpring {
    /// Create a new FENE spring.
    pub fn new(r_max: f64, k_spring: f64) -> Self {
        Self { r_max, k_spring }
    }

    /// Compute the FENE restoring force magnitude at extension `r`.
    ///
    /// `F = -k r / (1 - (r/R_max)²)`
    ///
    /// Returns a *negative* value (restoring), saturating toward ±∞ as
    /// r → R_max.  Returns 0 when r ≤ 0 and the simulation clamps r < R_max.
    ///
    /// # Arguments
    /// * `r` – Bond extension (must be < `r_max` for physical results).
    pub fn compute_force(&self, r: f64) -> f64 {
        if r <= 0.0 {
            return 0.0;
        }
        let ratio = r / self.r_max;
        let denom = 1.0 - ratio * ratio;
        if denom <= 0.0 {
            // Bond stretched beyond R_max – return a very large restoring force
            return -1e30 * r.signum();
        }
        -self.k_spring * r / denom
    }

    /// FENE potential energy U(r) = -½ k R_max² ln(1 - (r/R_max)²).
    ///
    /// Returns a very large value when r ≥ R_max.
    pub fn potential_energy(&self, r: f64) -> f64 {
        if r <= 0.0 {
            return 0.0;
        }
        let ratio = r / self.r_max;
        let arg = 1.0 - ratio * ratio;
        if arg <= 0.0 {
            return 1e30;
        }
        -0.5 * self.k_spring * self.r_max * self.r_max * arg.ln()
    }
}

// ---------------------------------------------------------------------------
// KremerGrestChain
// ---------------------------------------------------------------------------

/// Kremer–Grest bead-spring polymer chain.
///
/// Combines FENE bonds between adjacent beads with WCA (purely repulsive
/// Lennard-Jones) excluded-volume interactions between all bead pairs.
#[derive(Debug, Clone)]
pub struct KremerGrestChain {
    /// Bead degrees of freedom.
    pub beads: Vec<PolymerBead>,
    /// FENE springs connecting adjacent beads.
    pub fene_springs: Vec<FeneSpring>,
    /// LJ σ for bead–bead excluded volume.
    pub lj_sigma: f64,
    /// LJ ε for bead–bead excluded volume.
    pub lj_epsilon: f64,
}

impl KremerGrestChain {
    /// Create a new Kremer–Grest chain with `n` beads, stretched along x.
    ///
    /// # Arguments
    /// * `n`           – Number of beads.
    /// * `r_max`       – FENE maximum extension (typically 1.5 σ).
    /// * `k_fene`      – FENE spring constant.
    /// * `lj_sigma`    – LJ σ.
    /// * `lj_epsilon`  – LJ ε.
    pub fn new(n: usize, r_max: f64, k_fene: f64, lj_sigma: f64, lj_epsilon: f64) -> Self {
        let beads = (0..n)
            .map(|i| PolymerBead::new([i as f64 * lj_sigma, 0.0, 0.0], 1.0, 0.0))
            .collect();
        let fene_springs = (0..n.saturating_sub(1))
            .map(|_| FeneSpring::new(r_max, k_fene))
            .collect();
        Self {
            beads,
            fene_springs,
            lj_sigma,
            lj_epsilon,
        }
    }

    /// Number of beads.
    pub fn n_beads(&self) -> usize {
        self.beads.len()
    }

    /// Perform one velocity-Verlet integration step.
    ///
    /// # Arguments
    /// * `dt` – Time step.
    pub fn step(&mut self, dt: f64) {
        let forces = self.compute_forces();
        // Update positions and half-step velocities
        for (bead, f) in self.beads.iter_mut().zip(forces.iter()) {
            let inv_m = 1.0 / bead.mass;
            for (v, (p, &fk)) in bead
                .velocity
                .iter_mut()
                .zip(bead.position.iter_mut().zip(f.iter()))
            {
                *v += 0.5 * fk * inv_m * dt;
                *p += *v * dt;
            }
        }
        // Second force evaluation
        let forces2 = self.compute_forces();
        for (bead, f2) in self.beads.iter_mut().zip(forces2.iter()) {
            let inv_m = 1.0 / bead.mass;
            for (v, &fk) in bead.velocity.iter_mut().zip(f2.iter()) {
                *v += 0.5 * fk * inv_m * dt;
            }
        }
    }

    /// Compute total force on each bead (FENE bonds + WCA repulsion).
    fn compute_forces(&self) -> Vec<[f64; 3]> {
        let n = self.beads.len();
        let mut forces = vec![[0.0_f64; 3]; n];

        // FENE bonds between adjacent beads
        for i in 0..n.saturating_sub(1) {
            let ri = self.beads[i].position;
            let rj = self.beads[i + 1].position;
            let dr = sub3(rj, ri);
            let r = mag3(dr).max(1e-12);
            let f_mag = self.fene_springs[i].compute_force(r);
            for k in 0..3 {
                let fk = f_mag * dr[k] / r;
                forces[i][k] += fk;
                forces[i + 1][k] -= fk;
            }
        }

        // WCA excluded-volume (truncated LJ at r_cut = 2^(1/6) σ)
        let r_cut = 2.0_f64.powf(1.0 / 6.0) * self.lj_sigma;
        let eps = self.lj_epsilon;
        let sig = self.lj_sigma;
        for i in 0..n {
            for j in (i + 2)..n {
                let dr = sub3(self.beads[j].position, self.beads[i].position);
                let r = mag3(dr).max(1e-12);
                if r < r_cut {
                    let s6 = (sig / r).powi(6);
                    let f_mag = 24.0 * eps / (r * r) * (2.0 * s6 * s6 - s6);
                    for k in 0..3 {
                        let fk = f_mag * dr[k];
                        forces[i][k] -= fk;
                        forces[j][k] += fk;
                    }
                }
            }
        }

        forces
    }

    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.beads.iter().map(|b| b.kinetic_energy()).sum()
    }
}

// ---------------------------------------------------------------------------
// WormlikeChain
// ---------------------------------------------------------------------------

/// Wormlike Chain (WLC) model for semi-flexible polymers.
///
/// Uses the Marko–Siggia interpolation formula for the force-extension relation.
#[derive(Debug, Clone)]
pub struct WormlikeChain {
    /// Contour length L_c (m).
    pub contour_length: f64,
    /// Persistence length L_p (m).
    pub persistence_length: f64,
    /// Boltzmann thermal energy k_B T (J).
    pub kbt: f64,
}

impl WormlikeChain {
    /// Create a new WLC.
    ///
    /// # Arguments
    /// * `contour_length`    – Contour length L_c (m).
    /// * `persistence_length`– Persistence length L_p (m).
    /// * `kbt`               – Thermal energy k_B T (J).
    pub fn new(contour_length: f64, persistence_length: f64, kbt: f64) -> Self {
        Self {
            contour_length,
            persistence_length,
            kbt,
        }
    }

    /// End-to-end distance for a freely-jointed WLC at zero force.
    ///
    /// `R = sqrt(2 L_p L_c)` (valid for L_c >> L_p).
    pub fn end_to_end_distance(&self) -> f64 {
        (2.0 * self.persistence_length * self.contour_length).sqrt()
    }

    /// Marko–Siggia force-extension formula.
    ///
    /// `F = (k_B T / L_p) [ 1/(4(1 - x)²) - 1/4 + x ]`
    ///
    /// where x = extension / L_c ∈ \[0, 1).
    ///
    /// Returns 0 for non-positive extension and a large value when
    /// extension ≥ L_c.
    ///
    /// # Arguments
    /// * `extension` – End-to-end distance z (m), must be < L_c.
    pub fn force_extension(&self, extension: f64) -> f64 {
        if extension <= 0.0 {
            return 0.0;
        }
        let lc = self.contour_length;
        let lp = self.persistence_length;
        if extension >= lc {
            return 1e30; // Diverges as x → 1
        }
        let x = extension / lc;
        let denom = 1.0 - x;
        self.kbt / lp * (0.25 / (denom * denom) - 0.25 + x)
    }
}

// ---------------------------------------------------------------------------
// RadiusOfGyration
// ---------------------------------------------------------------------------

/// Utility for computing the radius of gyration of a set of beads.
///
/// R_g² = (1/N) Σ |r_i - r_cm|²
#[derive(Debug, Clone, Default)]
pub struct RadiusOfGyration;

impl RadiusOfGyration {
    /// Compute R_g from a slice of [`PolymerBead`].
    ///
    /// Returns 0 for empty or single-bead systems.
    pub fn compute(beads: &[PolymerBead]) -> f64 {
        let n = beads.len();
        if n <= 1 {
            return 0.0;
        }
        // Centre of mass (equal mass assumed for simplicity)
        let inv_n = 1.0 / n as f64;
        let cm: [f64; 3] = {
            let mut s = [0.0_f64; 3];
            for b in beads {
                s[0] += b.position[0];
                s[1] += b.position[1];
                s[2] += b.position[2];
            }
            [s[0] * inv_n, s[1] * inv_n, s[2] * inv_n]
        };
        let sum_sq: f64 = beads
            .iter()
            .map(|b| {
                let dr = sub3(b.position, cm);
                dr[0].powi(2) + dr[1].powi(2) + dr[2].powi(2)
            })
            .sum();
        (sum_sq * inv_n).sqrt()
    }
}

// ---------------------------------------------------------------------------
// RouseModel
// ---------------------------------------------------------------------------

/// Rouse bead-spring model for polymer dynamics in a viscous solvent.
///
/// Describes the normal-mode (Rouse) relaxation spectrum of an ideal chain.
#[derive(Debug, Clone)]
pub struct RouseModel {
    /// Number of Kuhn segments N.
    pub n_segments: usize,
    /// Kuhn segment length b (m).
    pub kuhn_length: f64,
    /// Friction coefficient per bead ζ (kg s⁻¹).
    pub friction_coeff: f64,
    /// Thermal energy k_B T (J).
    pub kbt: f64,
}

impl RouseModel {
    /// Create a new Rouse-model polymer.
    pub fn new(n_segments: usize, kuhn_length: f64, friction_coeff: f64, kbt: f64) -> Self {
        Self {
            n_segments,
            kuhn_length,
            friction_coeff,
            kbt,
        }
    }

    /// Relaxation time of Rouse mode p.
    ///
    /// `τ_p = (ζ N² b²) / (3 π² k_B T p²)`
    ///
    /// # Arguments
    /// * `p` – Mode index (1 = longest mode).
    ///
    /// Returns the Rouse relaxation time τ_p (s).
    pub fn rouse_mode_relaxation(&self, p: usize) -> f64 {
        if p == 0 {
            return f64::INFINITY;
        }
        let n = self.n_segments as f64;
        let b = self.kuhn_length;
        let zeta = self.friction_coeff;
        let kbt = self.kbt;
        let pp = p as f64;
        (zeta * n * n * b * b) / (3.0 * std::f64::consts::PI * std::f64::consts::PI * kbt * pp * pp)
    }

    /// Centre-of-mass diffusion coefficient D = k_B T / (N ζ).
    pub fn diffusion_coefficient(&self) -> f64 {
        let n = self.n_segments as f64;
        self.kbt / (n * self.friction_coeff)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- free functions ------------------------------------------------------

    #[test]
    fn test_flory_exponent_value() {
        assert!((flory_exponent_good_solvent() - 0.588).abs() < 1e-12);
    }

    #[test]
    fn test_end_to_end_ideal_100_beads() {
        // R = sqrt(100) * 1.0 = 10.0
        assert!((end_to_end_ideal(100, 1.0) - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_end_to_end_ideal_scales_with_sqrt_n() {
        let r1 = end_to_end_ideal(4, 1.0);
        let r4 = end_to_end_ideal(16, 1.0);
        assert!((r4 - 2.0 * r1).abs() < 1e-10);
    }

    #[test]
    fn test_end_to_end_ideal_scales_with_b() {
        let r1 = end_to_end_ideal(9, 1.0);
        let r2 = end_to_end_ideal(9, 2.0);
        assert!((r2 - 2.0 * r1).abs() < 1e-10);
    }

    #[test]
    fn test_rg_ideal_6_segments() {
        // Rg = 1 * sqrt(6/6) = 1.0
        assert!((rg_ideal(6, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rg_ideal_scales_with_sqrt_n() {
        let rg1 = rg_ideal(4, 1.0);
        let rg4 = rg_ideal(16, 1.0);
        assert!((rg4 - 2.0 * rg1).abs() < 1e-10);
    }

    #[test]
    fn test_rg_ideal_ratio_to_end_to_end() {
        // R_g = R_ee / sqrt(6)
        let n = 100;
        let b = 1.0;
        let ree = end_to_end_ideal(n, b);
        let rg = rg_ideal(n, b);
        assert!((rg - ree / 6.0_f64.sqrt()).abs() < 1e-10);
    }

    // -- PolymerBead ---------------------------------------------------------

    #[test]
    fn test_polymer_bead_at_rest() {
        let bead = PolymerBead::new([1.0, 2.0, 3.0], 1.0, 0.0);
        assert_eq!(bead.velocity, [0.0; 3]);
    }

    #[test]
    fn test_polymer_bead_kinetic_energy_zero() {
        let bead = PolymerBead::new([0.0; 3], 2.0, 0.0);
        assert_eq!(bead.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_polymer_bead_kinetic_energy_nonzero() {
        let mut bead = PolymerBead::new([0.0; 3], 2.0, 0.0);
        bead.velocity = [1.0, 0.0, 0.0];
        assert!((bead.kinetic_energy() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_polymer_bead_charge_stored() {
        let bead = PolymerBead::new([0.0; 3], 1.0, -1.6e-19);
        assert!((bead.charge + 1.6e-19).abs() < 1e-30);
    }

    // -- FeneSpring ----------------------------------------------------------

    #[test]
    fn test_fene_force_zero_extension() {
        let spring = FeneSpring::new(1.5, 30.0);
        assert_eq!(spring.compute_force(0.0), 0.0);
    }

    #[test]
    fn test_fene_force_negative_restoring() {
        let spring = FeneSpring::new(1.5, 30.0);
        let f = spring.compute_force(0.5);
        assert!(f < 0.0, "FENE force should be restoring (negative)");
    }

    #[test]
    fn test_fene_force_magnitude_increases_with_r() {
        let spring = FeneSpring::new(1.5, 30.0);
        let f1 = spring.compute_force(0.3).abs();
        let f2 = spring.compute_force(0.8).abs();
        assert!(f2 > f1);
    }

    #[test]
    fn test_fene_force_diverges_near_rmax() {
        let spring = FeneSpring::new(1.0, 30.0);
        let f = spring.compute_force(0.999).abs();
        assert!(f > 1000.0);
    }

    #[test]
    fn test_fene_potential_zero_extension() {
        let spring = FeneSpring::new(1.5, 30.0);
        assert_eq!(spring.potential_energy(0.0), 0.0);
    }

    #[test]
    fn test_fene_potential_positive() {
        let spring = FeneSpring::new(1.5, 30.0);
        let u = spring.potential_energy(0.5);
        assert!(u > 0.0);
    }

    #[test]
    fn test_fene_potential_increases_with_r() {
        let spring = FeneSpring::new(1.5, 30.0);
        let u1 = spring.potential_energy(0.3);
        let u2 = spring.potential_energy(0.8);
        assert!(u2 > u1);
    }

    // -- WormlikeChain -------------------------------------------------------

    #[test]
    fn test_wlc_end_to_end_formula() {
        let wlc = WormlikeChain::new(1000e-9, 50e-9, 4.1e-21);
        let r = wlc.end_to_end_distance();
        let expected = (2.0 * 50e-9 * 1000e-9_f64).sqrt();
        assert!((r - expected).abs() < 1e-20);
    }

    #[test]
    fn test_wlc_force_zero_extension() {
        let wlc = WormlikeChain::new(1.0, 0.05, 4.1e-21);
        assert_eq!(wlc.force_extension(0.0), 0.0);
    }

    #[test]
    fn test_wlc_force_small_extension_positive() {
        let wlc = WormlikeChain::new(1.0, 0.05, 4.1e-21);
        let f = wlc.force_extension(0.1);
        assert!(f > 0.0);
    }

    #[test]
    fn test_wlc_force_increases_with_extension() {
        let wlc = WormlikeChain::new(1.0, 0.05, 4.1e-21);
        let f1 = wlc.force_extension(0.3);
        let f2 = wlc.force_extension(0.7);
        assert!(f2 > f1);
    }

    #[test]
    fn test_wlc_force_near_contour_large() {
        // Use kBT = 1.0 (reduced units) so force is easily > 1.0 near contour
        let wlc = WormlikeChain::new(1.0, 0.05, 1.0);
        let f = wlc.force_extension(0.999);
        assert!(f > 1.0);
    }

    #[test]
    fn test_wlc_force_at_contour_length_saturates() {
        let wlc = WormlikeChain::new(1.0, 0.05, 4.1e-21);
        let f = wlc.force_extension(1.0);
        assert!(f >= 1e30);
    }

    // -- RadiusOfGyration ----------------------------------------------------

    #[test]
    fn test_rg_single_bead_zero() {
        let beads = vec![PolymerBead::new([0.0; 3], 1.0, 0.0)];
        assert_eq!(RadiusOfGyration::compute(&beads), 0.0);
    }

    #[test]
    fn test_rg_two_beads_symmetric() {
        let beads = vec![
            PolymerBead::new([-1.0, 0.0, 0.0], 1.0, 0.0),
            PolymerBead::new([1.0, 0.0, 0.0], 1.0, 0.0),
        ];
        // cm = 0, sum_sq = 2, Rg = sqrt(2/2) = 1.0
        assert!((RadiusOfGyration::compute(&beads) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rg_four_beads_square() {
        // Beads at corners of a unit square in xy-plane
        let beads = vec![
            PolymerBead::new([1.0, 1.0, 0.0], 1.0, 0.0),
            PolymerBead::new([-1.0, 1.0, 0.0], 1.0, 0.0),
            PolymerBead::new([1.0, -1.0, 0.0], 1.0, 0.0),
            PolymerBead::new([-1.0, -1.0, 0.0], 1.0, 0.0),
        ];
        // cm = 0, each |dr|² = 2, Rg = sqrt(4*2/4) = sqrt(2)
        let rg = RadiusOfGyration::compute(&beads);
        assert!((rg - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_rg_empty_zero() {
        let beads: Vec<PolymerBead> = vec![];
        assert_eq!(RadiusOfGyration::compute(&beads), 0.0);
    }

    // -- RouseModel ----------------------------------------------------------

    #[test]
    fn test_rouse_relaxation_mode1() {
        let rm = RouseModel::new(100, 1.0, 1.0, 1.0);
        let tau = rm.rouse_mode_relaxation(1);
        let expected = (1.0 * 100.0 * 100.0 * 1.0 * 1.0)
            / (3.0 * std::f64::consts::PI * std::f64::consts::PI * 1.0 * 1.0);
        assert!((tau - expected).abs() < 1e-10);
    }

    #[test]
    fn test_rouse_relaxation_mode_p_scales_inverse_p2() {
        let rm = RouseModel::new(100, 1.0, 1.0, 1.0);
        let tau1 = rm.rouse_mode_relaxation(1);
        let tau2 = rm.rouse_mode_relaxation(2);
        assert!((tau2 - tau1 / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_rouse_relaxation_mode0_infinity() {
        let rm = RouseModel::new(100, 1.0, 1.0, 1.0);
        assert_eq!(rm.rouse_mode_relaxation(0), f64::INFINITY);
    }

    #[test]
    fn test_rouse_diffusion_coefficient() {
        let rm = RouseModel::new(10, 1.0, 1.0, 1.0);
        let d = rm.diffusion_coefficient();
        // D = kBT / (N * zeta) = 1 / (10 * 1) = 0.1
        assert!((d - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_rouse_diffusion_scales_inverse_n() {
        let rm1 = RouseModel::new(10, 1.0, 1.0, 1.0);
        let rm2 = RouseModel::new(20, 1.0, 1.0, 1.0);
        assert!((rm2.diffusion_coefficient() - 0.5 * rm1.diffusion_coefficient()).abs() < 1e-12);
    }

    #[test]
    fn test_rouse_relaxation_scales_with_n_squared() {
        let rm1 = RouseModel::new(10, 1.0, 1.0, 1.0);
        let rm2 = RouseModel::new(20, 1.0, 1.0, 1.0);
        let ratio = rm2.rouse_mode_relaxation(1) / rm1.rouse_mode_relaxation(1);
        assert!((ratio - 4.0).abs() < 1e-10);
    }

    // -- KremerGrestChain ----------------------------------------------------

    #[test]
    fn test_kremer_grest_chain_creation() {
        let chain = KremerGrestChain::new(5, 1.5, 30.0, 1.0, 1.0);
        assert_eq!(chain.n_beads(), 5);
        assert_eq!(chain.fene_springs.len(), 4);
    }

    #[test]
    fn test_kremer_grest_step_runs() {
        let mut chain = KremerGrestChain::new(4, 1.5, 30.0, 1.0, 1.0);
        chain.step(1e-4);
        assert_eq!(chain.n_beads(), 4);
    }

    #[test]
    fn test_kremer_grest_kinetic_energy_zero_initially() {
        let chain = KremerGrestChain::new(5, 1.5, 30.0, 1.0, 1.0);
        assert!((chain.kinetic_energy() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_kremer_grest_beads_move_after_step() {
        // Use an asymmetric initial chain so the end beads experience a net force
        let mut chain = KremerGrestChain::new(3, 1.5, 30.0, 1.0, 1.0);
        // Displace the first bead to break symmetry
        chain.beads[0].position[1] = 0.3;
        let pos_before = chain.beads[0].position;
        chain.step(1e-3);
        let pos_after = chain.beads[0].position;
        let moved = pos_before
            .iter()
            .zip(pos_after.iter())
            .any(|(a, b)| (a - b).abs() > 1e-15);
        assert!(moved);
    }
}
