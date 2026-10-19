// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum chemistry molecular dynamics (Born-Oppenheimer MD, Car-Parrinello).
//!
//! This module provides:
//! - [`HartreeFockMd`]: RHF/STO-3G Born-Oppenheimer MD (SCF energies, finite-difference forces)
//! - [`DensityFunctionalMd`]: real-space grid Kohn-Sham DFT Born-Oppenheimer MD
//!   (LDA, PBE, B3LYP) with free-space FFT Poisson and LOBPCG diagonalisation
//! - [`CarParrinelloMd`]: Car-Parrinello extended Lagrangian MD
//! - [`SemiempiricalMd`]: model core-core pair potentials flavoured after AM1, PM3, xTB and DFTB
//! - [`NebMethod`]: Nudged elastic band for reaction path finding

mod dft_md;
#[cfg(test)]
mod dft_md_tests;
mod error;
mod hartree_fock_md;
mod ks_eigensolver;
mod neb;
mod real_space_grid;
mod semiempirical;
mod xc;

pub use dft_md::{DensityFunctionalMd, KsEnergies, KsScfReport, KsStepReport};
pub use error::QuantumChemistryMdError;
pub use hartree_fock_md::{AMU_TO_ELECTRON_MASS, HartreeFockMd};
pub use neb::NebMethod;
pub use real_space_grid::UniformGrid;
pub use semiempirical::SemiempiricalMd;

// ---------------------------------------------------------------------------
// Enumerations
// ---------------------------------------------------------------------------

/// Exchange-correlation functional type for [`DensityFunctionalMd`].
#[derive(Debug, Clone, PartialEq)]
pub enum XcFunctional {
    /// Local density approximation: Slater exchange + VWN5 correlation.
    Lda,
    /// Generalized gradient approximation — PBE exchange and correlation.
    Pbe,
    /// Hybrid functional B3LYP (20 % exact exchange, B88, VWN-RPA, LYP).
    B3lyp,
}

/// Family of model pair potential used by [`SemiempiricalMd`].
///
/// The variants select fixed model parameters named after the semiempirical
/// method they imitate; they are not the parameterised NDDO/tight-binding
/// Hamiltonians themselves.
#[derive(Debug, Clone, PartialEq)]
pub enum SemiempiricalMethod {
    /// Austin Model 1.
    Am1,
    /// Parameterized Model 3.
    Pm3,
    /// Extended tight-binding method (GFN2-xTB).
    Xtb,
    /// Density functional tight binding.
    Dftb,
}

// ---------------------------------------------------------------------------
// CarParrinelloMd
// ---------------------------------------------------------------------------

/// Car-Parrinello extended-Lagrangian MD driver.
///
/// Propagates both ionic positions and electronic degrees of freedom
/// (Kohn-Sham orbitals ψ) simultaneously using a fictitious electron mass μ.
/// Orthonormality of orbitals is enforced via Gram-Schmidt.
#[derive(Debug, Clone)]
pub struct CarParrinelloMd {
    /// Fictitious electron mass μ in atomic units.
    pub fictitious_mass: f64,
    /// Kohn-Sham orbitals: `psi[i]` is orbital i sampled on a grid.
    pub psi: Vec<Vec<f64>>,
    /// Ionic positions \[Bohr\].
    pub ionic_positions: Vec<[f64; 3]>,
    /// Ionic velocities \[Bohr/a.u.\].
    pub ionic_velocities: Vec<[f64; 3]>,
    /// Orbital velocities (time-derivative of ψ coefficients).
    pub psi_dot: Vec<Vec<f64>>,
    /// Ionic masses \[a.u.\].
    pub ionic_masses: Vec<f64>,
}

impl CarParrinelloMd {
    /// Create a Car-Parrinello MD driver.
    ///
    /// # Arguments
    /// * `fictitious_mass` – μ, fictitious electron mass (typical: 400–900 a.u.).
    /// * `psi` – Initial orbital coefficients.
    /// * `ionic_positions` – Initial ionic positions.
    /// * `ionic_masses` – Ionic masses.
    pub fn new(
        fictitious_mass: f64,
        psi: Vec<Vec<f64>>,
        ionic_positions: Vec<[f64; 3]>,
        ionic_masses: Vec<f64>,
    ) -> Self {
        let n_orb = psi.len();
        let n_coeff = psi.first().map(|v| v.len()).unwrap_or(0);
        let n_ions = ionic_positions.len();
        Self {
            fictitious_mass,
            psi,
            ionic_positions,
            ionic_velocities: vec![[0.0; 3]; n_ions],
            psi_dot: vec![vec![0.0; n_coeff]; n_orb],
            ionic_masses,
        }
    }

    /// Advance ionic and orbital degrees of freedom by `dt` (velocity Verlet).
    ///
    /// # Arguments
    /// * `dt` – Time step in atomic units.
    pub fn step(&mut self, dt: f64) {
        // Update ionic positions
        let forces = self.ionic_forces();
        // Update ionic positions
        for ((pos, vel), (f, &m)) in self
            .ionic_positions
            .iter_mut()
            .zip(self.ionic_velocities.iter())
            .zip(forces.iter().zip(self.ionic_masses.iter()))
        {
            for k in 0..3 {
                pos[k] += dt * vel[k] + 0.5 * dt * dt * f[k] / m;
            }
        }
        // Update orbital coefficients
        for (orb, orb_dot) in self.psi.iter_mut().zip(self.psi_dot.iter()) {
            for (c, cd) in orb.iter_mut().zip(orb_dot.iter()) {
                *c += dt * cd;
            }
        }
        self.orthogonalize_orbitals();
        // Update ionic velocities (simple Euler for now)
        let forces2 = self.ionic_forces();
        for ((vel, &m), (f, f2)) in self
            .ionic_velocities
            .iter_mut()
            .zip(self.ionic_masses.iter())
            .zip(forces.iter().zip(forces2.iter()))
        {
            for k in 0..3 {
                vel[k] += 0.5 * dt * (f[k] + f2[k]) / m;
            }
        }
    }

    /// Kinetic energy of the fictitious electron degrees of freedom.
    ///
    /// T_e = (μ/2) Σ_i ∫ |ψ̇_i|² dΩ
    pub fn kinetic_energy_electrons(&self) -> f64 {
        let mut ke = 0.0;
        for orb_dot in &self.psi_dot {
            ke += orb_dot.iter().map(|&x| x * x).sum::<f64>();
        }
        0.5 * self.fictitious_mass * ke
    }

    /// Orthonormalize Kohn-Sham orbitals using Gram-Schmidt.
    pub fn orthogonalize_orbitals(&mut self) {
        let n = self.psi.len();
        for i in 0..n {
            // Subtract projections onto already-orthogonalized orbitals
            for j in 0..i {
                let dot: f64 = self.psi[i]
                    .iter()
                    .zip(self.psi[j].iter())
                    .map(|(a, b)| a * b)
                    .sum();
                let psi_j: Vec<f64> = self.psi[j].clone();
                for (c, pj) in self.psi[i].iter_mut().zip(psi_j.iter()) {
                    *c -= dot * pj;
                }
            }
            // Normalize
            let norm: f64 = self.psi[i].iter().map(|&x| x * x).sum::<f64>().sqrt();
            if norm > 1e-12 {
                for c in &mut self.psi[i] {
                    *c /= norm;
                }
            }
        }
    }

    /// Pairwise nuclear-repulsion forces on the ions.
    ///
    /// Returns the force vector `Vec<[f64; 3]>` (one 3-vector per ion). This is
    /// a deliberately **simplified** model: it sums a pairwise repulsive ~1/r²
    /// interaction between the ionic centres,
    ///
    /// F_i = Σ_{j≠i} (R_i − R_j) / |R_i − R_j|³,
    ///
    /// the gradient of a 1/r repulsion with unit coupling. It is *not* a full
    /// Hellmann-Feynman force: the electronic contribution from the orbitals
    /// `psi` and the true nuclear charges are not included (they cannot be
    /// recovered from the stored coefficient data), so callers needing
    /// ab-initio forces must supply them from a real electronic-structure
    /// evaluation.
    pub fn ionic_forces(&self) -> Vec<[f64; 3]> {
        let n = self.ionic_positions.len();
        let mut f = vec![[0.0f64; 3]; n];
        for (i, f_i) in f.iter_mut().enumerate() {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = dist(&self.ionic_positions[i], &self.ionic_positions[j]).max(1e-6);
                // Repulsive gradient of a unit 1/r interaction: +(R_i - R_j)/r^3.
                let coeff = 1.0 / (r * r * r);
                for (fk, (&pi_k, &pj_k)) in f_i.iter_mut().zip(
                    self.ionic_positions[i]
                        .iter()
                        .zip(self.ionic_positions[j].iter()),
                ) {
                    *fk += coeff * (pi_k - pj_k);
                }
            }
        }
        f
    }
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

/// Euclidean distance between two 3-D points.
fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CarParrinelloMd ---

    #[test]
    fn test_cp_orthogonalize_single_orbital() {
        let mut cp = CarParrinelloMd::new(
            500.0,
            vec![vec![3.0, 4.0]],
            vec![[0.0, 0.0, 0.0]],
            vec![12.0],
        );
        cp.orthogonalize_orbitals();
        let norm: f64 = cp.psi[0].iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cp_orthogonalize_two_orbitals() {
        let mut cp = CarParrinelloMd::new(
            500.0,
            vec![vec![1.0, 0.0, 0.0], vec![1.0, 1.0, 0.0]],
            vec![[0.0, 0.0, 0.0]],
            vec![12.0],
        );
        cp.orthogonalize_orbitals();
        let dot: f64 = cp.psi[0]
            .iter()
            .zip(cp.psi[1].iter())
            .map(|(a, b)| a * b)
            .sum();
        assert!(dot.abs() < 1e-10);
    }

    #[test]
    fn test_cp_electron_kinetic_energy_zero_velocity() {
        let cp = CarParrinelloMd::new(
            500.0,
            vec![vec![1.0, 0.0]],
            vec![[0.0, 0.0, 0.0]],
            vec![12.0],
        );
        assert_eq!(cp.kinetic_energy_electrons(), 0.0);
    }

    #[test]
    fn test_cp_step_position_changes() {
        let mut cp = CarParrinelloMd::new(
            500.0,
            vec![vec![0.0, 1.0]],
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        cp.ionic_velocities[0] = [0.01, 0.0, 0.0];
        let pos_before = cp.ionic_positions[0][0];
        cp.step(0.1);
        let pos_after = cp.ionic_positions[0][0];
        assert!((pos_after - pos_before).abs() > 1e-6);
    }

    #[test]
    fn test_cp_ionic_forces_public() {
        let cp = CarParrinelloMd::new(
            500.0,
            vec![vec![1.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        let f = cp.ionic_forces();
        assert_eq!(f.len(), 2);
        // Two ions on the x-axis: force must be non-zero, purely along x,
        // repulsive (ion 0 pushed toward -x, ion 1 toward +x) and obey
        // Newton's third law f0 = -f1.
        assert!(f[0][0].abs() > 1e-9);
        assert!(f[0][0] < 0.0);
        assert!(f[1][0] > 0.0);
        for (a, b) in f[0].iter().zip(f[1].iter()) {
            assert!((a + b).abs() < 1e-12);
        }
        assert!(f[0][1].abs() < 1e-12 && f[0][2].abs() < 1e-12);
    }

    #[test]
    fn test_dist_utility() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 4.0, 0.0];
        assert!((dist(&a, &b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_cp_kinetic_energy_nonzero() {
        let mut cp = CarParrinelloMd::new(
            400.0,
            vec![vec![1.0, 0.0]],
            vec![[0.0, 0.0, 0.0]],
            vec![12.0],
        );
        cp.psi_dot[0] = vec![1.0, 0.0];
        let ke = cp.kinetic_energy_electrons();
        assert!((ke - 200.0).abs() < 1e-10);
    }
}
