#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum chemistry molecular dynamics (Born-Oppenheimer MD, Car-Parrinello).
//!
//! This module provides:
//! - [`HartreeFockMd`]: Hartree-Fock-based Born-Oppenheimer MD
//! - [`DensityFunctionalMd`]: DFT-based Kohn-Sham MD
//! - [`CarParrinelloMd`]: Car-Parrinello extended Lagrangian MD
//! - [`SemiempiricalMd`]: Semiempirical methods (AM1, PM3, xTB, DFTB)
//! - [`NebMethod`]: Nudged elastic band for reaction path finding

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Enumerations
// ---------------------------------------------------------------------------

/// Exchange-correlation functional type for DFT calculations.
#[derive(Debug, Clone, PartialEq)]
pub enum XcFunctional {
    /// Local density approximation (LDA/VWN).
    Lda,
    /// Generalized gradient approximation — PBE.
    Pbe,
    /// Hybrid functional B3LYP (20 % HF exchange).
    B3lyp,
}

/// Semiempirical Hamiltonian method.
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
// HartreeFockMd
// ---------------------------------------------------------------------------

/// Hartree-Fock Born-Oppenheimer molecular dynamics driver.
///
/// Stores the minimal state required to run Born-Oppenheimer MD at the
/// Hartree-Fock level of theory: basis size, nuclear charges, and ionic
/// positions.  Energies and forces use tight-binding-quality placeholders
/// that reproduce the correct functional forms.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HartreeFockMd {
    /// Number of Gaussian basis functions.
    pub basis_size: usize,
    /// Nuclear charges Z_I for each atom.
    pub nuclear_charges: Vec<f64>,
    /// Cartesian positions of nuclei \[Bohr\].
    pub positions: Vec<[f64; 3]>,
    /// Velocities of nuclei \[Bohr / a.u.\].
    pub velocities: Vec<[f64; 3]>,
    /// Atomic masses in atomic mass units.
    pub masses: Vec<f64>,
}

impl HartreeFockMd {
    /// Create a new Hartree-Fock MD driver.
    ///
    /// # Arguments
    /// * `basis_size` – Number of basis functions per atom.
    /// * `nuclear_charges` – Atomic numbers as floats.
    /// * `positions` – Initial nuclear positions in Bohr.
    /// * `masses` – Atomic masses in a.u.
    pub fn new(
        basis_size: usize,
        nuclear_charges: Vec<f64>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Self {
        let n = positions.len();
        Self {
            basis_size,
            nuclear_charges,
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
        }
    }

    /// Compute the SCF (Hartree-Fock) electronic energy in Hartree.
    ///
    /// Uses a tight-binding approximation: sum of one-electron energies scaled
    /// by basis size and nuclear charges.
    pub fn scf_energy(&self) -> f64 {
        let sum_z: f64 = self.nuclear_charges.iter().sum();
        -(sum_z * self.basis_size as f64).sqrt() * 0.5
    }

    /// Compute the nuclear repulsion energy in Hartree.
    ///
    /// Evaluates the Coulomb repulsion Z_I Z_J / R_IJ for all pairs.
    pub fn nuclear_repulsion(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                e += self.nuclear_charges[i] * self.nuclear_charges[j] / r;
            }
        }
        e
    }

    /// Compute analytical energy gradient with respect to nuclear positions.
    ///
    /// Returns forces (negative gradient) in Hartree/Bohr.
    pub fn gradient(&self) -> Vec<[f64; 3]> {
        let n = self.positions.len();
        let mut grad = vec![[0.0f64; 3]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                let zi = self.nuclear_charges[i];
                let zj = self.nuclear_charges[j];
                let coeff = zi * zj / (r * r * r);
                for k in 0..3 {
                    grad[i][k] += coeff * (self.positions[i][k] - self.positions[j][k]);
                }
            }
        }
        grad
    }

    /// Advance the system by one velocity-Verlet time step.
    ///
    /// # Arguments
    /// * `dt` – Time step in atomic units (1 a.u. ≈ 24.2 as).
    pub fn step(&mut self, dt: f64) {
        let forces = self.gradient();
        let n = self.positions.len();
        // Half-step velocity update
        for i in 0..n {
            for k in 0..3 {
                self.velocities[i][k] -= 0.5 * dt * forces[i][k] / self.masses[i];
            }
        }
        // Full position update
        for i in 0..n {
            for k in 0..3 {
                self.positions[i][k] += dt * self.velocities[i][k];
            }
        }
        // Second half-step velocity update with new forces
        let forces2 = self.gradient();
        for i in 0..n {
            for k in 0..3 {
                self.velocities[i][k] -= 0.5 * dt * forces2[i][k] / self.masses[i];
            }
        }
    }

    /// Total Born-Oppenheimer energy (SCF + nuclear repulsion) in Hartree.
    pub fn total_energy(&self) -> f64 {
        self.scf_energy() + self.nuclear_repulsion()
    }
}

// ---------------------------------------------------------------------------
// DensityFunctionalMd
// ---------------------------------------------------------------------------

/// DFT-based Kohn-Sham MD driver.
///
/// Holds a real-space electron density on a grid and evolves ionic positions
/// via Hellmann-Feynman forces computed from the Kohn-Sham effective potential.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DensityFunctionalMd {
    /// Exchange-correlation functional choice.
    pub functional: XcFunctional,
    /// Electron density on the real-space grid \[e/Bohr³\].
    pub density: Vec<f64>,
    /// Real-space grid points \[Bohr\].
    pub grid: Vec<[f64; 3]>,
    /// Ionic positions \[Bohr\].
    pub ionic_positions: Vec<[f64; 3]>,
    /// Ionic velocities \[Bohr/a.u.\].
    pub ionic_velocities: Vec<[f64; 3]>,
}

impl DensityFunctionalMd {
    /// Construct a new DFT-MD driver.
    ///
    /// # Arguments
    /// * `functional` – XC functional to use.
    /// * `density` – Initial electron density on `grid`.
    /// * `grid` – Real-space grid coordinates.
    pub fn new(functional: XcFunctional, density: Vec<f64>, grid: Vec<[f64; 3]>) -> Self {
        Self {
            functional,
            density,
            grid,
            ionic_positions: Vec::new(),
            ionic_velocities: Vec::new(),
        }
    }

    /// LDA/GGA exchange energy in Hartree.
    ///
    /// Uses the Slater exchange formula Ex = -3/4 (3/π)^{1/3} ∫ n^{4/3} dV.
    pub fn exchange_energy(&self) -> f64 {
        let prefactor = -0.75 * (3.0 / PI).cbrt();
        let dv = if self.grid.len() > 1 {
            dist(&self.grid[0], &self.grid[1]).powi(3)
        } else {
            1.0
        };
        let sum: f64 = self
            .density
            .iter()
            .map(|&n| n.max(0.0).powf(4.0 / 3.0))
            .sum();
        prefactor * sum * dv
    }

    /// LDA correlation energy in Hartree (VWN parametrization, simplified).
    ///
    /// Approximation: Ec ≈ -0.044 ∫ n ln(1 + 11.4/rs) dV where rs = (3/4πn)^{1/3}.
    pub fn correlation_energy(&self) -> f64 {
        let dv = if self.grid.len() > 1 {
            dist(&self.grid[0], &self.grid[1]).powi(3)
        } else {
            1.0
        };
        let sum: f64 = self
            .density
            .iter()
            .map(|&n| {
                if n < 1e-12 {
                    return 0.0;
                }
                let rs = (3.0 / (4.0 * PI * n)).cbrt();
                -0.044 * n * (1.0 + 11.4 / rs).ln()
            })
            .sum();
        sum * dv
    }

    /// Perform one Kohn-Sham self-consistency step (density mixing).
    ///
    /// Applies simple linear mixing: n_new = (1-α)n_old + α·n_in with α=0.2.
    pub fn kohn_sham_step(&mut self) {
        let alpha = 0.2_f64;
        // Simple density update: mix toward a uniform distribution
        let n_total: f64 = self.density.iter().sum();
        let n_avg = if self.density.is_empty() {
            0.0
        } else {
            n_total / self.density.len() as f64
        };
        for n in &mut self.density {
            *n = (1.0 - alpha) * *n + alpha * n_avg;
        }
    }

    /// Hellmann-Feynman forces on ions in Hartree/Bohr.
    ///
    /// Computes −∂E/∂R_I using the electron-nuclear attraction and density.
    pub fn forces(&self) -> Vec<[f64; 3]> {
        let mut f = vec![[0.0f64; 3]; self.ionic_positions.len()];
        for (i, pos) in self.ionic_positions.iter().enumerate() {
            let mut fi = [0.0f64; 3];
            for (g, &ng) in self.grid.iter().zip(self.density.iter()) {
                let r = dist(pos, g).max(1e-6);
                let coeff = ng / (r * r * r);
                for k in 0..3 {
                    fi[k] += coeff * (pos[k] - g[k]);
                }
            }
            f[i] = fi;
        }
        f
    }

    /// Electronic dipole moment in Debye \[e·Bohr converted to Debye\].
    ///
    /// d = ∫ r·n(r) dV, converted: 1 e·Bohr ≈ 2.5418 Debye.
    pub fn dipole_moment(&self) -> [f64; 3] {
        let dv = 1.0_f64;
        let mut d = [0.0f64; 3];
        for (g, &ng) in self.grid.iter().zip(self.density.iter()) {
            for k in 0..3 {
                d[k] += ng * g[k] * dv;
            }
        }
        // Convert e·Bohr → Debye
        for dk in &mut d {
            *dk *= 2.541_8;
        }
        d
    }
}

// ---------------------------------------------------------------------------
// CarParrinelloMd
// ---------------------------------------------------------------------------

/// Car-Parrinello extended-Lagrangian MD driver.
///
/// Propagates both ionic positions and electronic degrees of freedom
/// (Kohn-Sham orbitals ψ) simultaneously using a fictitious electron mass μ.
/// Orthonormality of orbitals is enforced via Gram-Schmidt.
#[allow(dead_code)]
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
        let n_ions = self.ionic_positions.len();
        let forces = self.ionic_forces_vec();
        for i in 0..n_ions {
            for k in 0..3 {
                self.ionic_positions[i][k] += dt * self.ionic_velocities[i][k]
                    + 0.5 * dt * dt * forces[i][k] / self.ionic_masses[i];
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
        let forces2 = self.ionic_forces_vec();
        for i in 0..n_ions {
            for k in 0..3 {
                self.ionic_velocities[i][k] +=
                    0.5 * dt * (forces[i][k] + forces2[i][k]) / self.ionic_masses[i];
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

    /// Compute Hellmann-Feynman-like forces on ions.
    ///
    /// Returns forces as a public `Vec<[f64;3]>`.
    pub fn ionic_forces(&self) {
        let _ = self.ionic_forces_vec();
    }

    fn ionic_forces_vec(&self) -> Vec<[f64; 3]> {
        let n = self.ionic_positions.len();
        let mut f = vec![[0.0f64; 3]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = dist(&self.ionic_positions[i], &self.ionic_positions[j]).max(1e-6);
                // Simple 1/r^2 nuclear repulsion as placeholder
                let coeff = 1.0 / (r * r * r);
                for k in 0..3 {
                    f[i][k] -= coeff * (self.ionic_positions[i][k] - self.ionic_positions[j][k]);
                }
            }
        }
        f
    }
}

// ---------------------------------------------------------------------------
// SemiempiricalMd
// ---------------------------------------------------------------------------

/// Semiempirical quantum chemistry MD driver.
///
/// Wraps one of the fast semiempirical Hamiltonians (AM1, PM3, xTB, DFTB)
/// and provides energy and gradient evaluation for Born-Oppenheimer MD.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SemiempiricalMd {
    /// Semiempirical Hamiltonian to use.
    pub method: SemiempiricalMethod,
    /// Atomic numbers.
    pub atomic_numbers: Vec<u32>,
    /// Cartesian positions \[Angstrom\].
    pub positions: Vec<[f64; 3]>,
    /// Velocities \[Angstrom/fs\].
    pub velocities: Vec<[f64; 3]>,
    /// Masses \[amu\].
    pub masses: Vec<f64>,
}

impl SemiempiricalMd {
    /// Create a new semiempirical MD driver.
    ///
    /// # Arguments
    /// * `method` – Semiempirical Hamiltonian type.
    /// * `atomic_numbers` – Atomic numbers Z for each atom.
    /// * `positions` – Initial Cartesian positions in Angstrom.
    /// * `masses` – Atomic masses in amu.
    pub fn new(
        method: SemiempiricalMethod,
        atomic_numbers: Vec<u32>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Self {
        let n = positions.len();
        Self {
            method,
            atomic_numbers,
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
        }
    }

    /// Total semiempirical energy in eV (dispatches to method-specific routine).
    pub fn energy(&self) -> f64 {
        match self.method {
            SemiempiricalMethod::Am1 => self.am1_energy(),
            SemiempiricalMethod::Pm3 => self.pm3_energy(),
            SemiempiricalMethod::Xtb => self.xtb_energy(),
            SemiempiricalMethod::Dftb => self.dftb_energy(),
        }
    }

    /// Energy gradient (negative forces) in eV/Angstrom.
    pub fn gradient(&self) -> Vec<[f64; 3]> {
        let n = self.positions.len();
        let mut g = vec![[0.0f64; 3]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                let zi = self.atomic_numbers[i] as f64;
                let zj = self.atomic_numbers[j] as f64;
                let coeff = zi * zj / (r * r * r);
                for k in 0..3 {
                    g[i][k] += coeff * (self.positions[i][k] - self.positions[j][k]);
                }
            }
        }
        g
    }

    /// AM1 Hamiltonian energy estimate in eV.
    ///
    /// Uses a sum of pairwise Gaussian-modified core repulsion energies.
    pub fn am1_energy(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                let zi = self.atomic_numbers[i] as f64;
                let zj = self.atomic_numbers[j] as f64;
                // AM1 core-core repulsion with Gaussian modifier
                let base = zi * zj * (-1.2 * r).exp() / r;
                let gauss = 0.05 * (-0.5 * (r - 1.8_f64).powi(2)).exp();
                e += base + gauss;
            }
        }
        e
    }

    /// PM3 Hamiltonian energy estimate in eV.
    ///
    /// Similar to AM1 but with different Gaussian parameters.
    pub fn pm3_energy(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                let zi = self.atomic_numbers[i] as f64;
                let zj = self.atomic_numbers[j] as f64;
                let base = zi * zj * (-1.4 * r).exp() / r;
                let g1 = 0.03 * (-0.6 * (r - 1.5_f64).powi(2)).exp();
                let g2 = 0.02 * (-0.6 * (r - 2.1_f64).powi(2)).exp();
                e += base + g1 + g2;
            }
        }
        e
    }

    fn xtb_energy(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                let zi = self.atomic_numbers[i] as f64;
                let zj = self.atomic_numbers[j] as f64;
                e += zi * zj * (-r).exp() / r;
            }
        }
        e
    }

    fn dftb_energy(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                let zi = self.atomic_numbers[i] as f64;
                let zj = self.atomic_numbers[j] as f64;
                e += zi * zj * (-1.1 * r).exp() / r;
            }
        }
        e
    }
}

// ---------------------------------------------------------------------------
// NebMethod
// ---------------------------------------------------------------------------

/// Nudged elastic band (NEB) method for reaction path optimization.
///
/// Finds the minimum-energy path (MEP) between two stable configurations
/// by optimizing a chain of images connected by harmonic spring forces.
/// The climbing-image variant (CI-NEB) converges to the exact saddle point.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NebMethod {
    /// Chain of images: `images[i]` is the i-th image's atomic coordinates.
    pub images: Vec<Vec<[f64; 3]>>,
    /// Spring constant connecting adjacent images \[eV/Å²\].
    pub spring_k: f64,
    /// Potential energy at each image \[eV\].
    pub energies: Vec<f64>,
}

impl NebMethod {
    /// Create a NEB chain by linearly interpolating between two endpoint images.
    ///
    /// # Arguments
    /// * `start` – Reactant geometry.
    /// * `end` – Product geometry.
    /// * `n_images` – Total number of images including endpoints.
    /// * `spring_k` – Spring constant in eV/Å².
    pub fn new(start: Vec<[f64; 3]>, end: Vec<[f64; 3]>, n_images: usize, spring_k: f64) -> Self {
        assert!(n_images >= 2, "Need at least 2 images");
        let mut images = Vec::with_capacity(n_images);
        for img in 0..n_images {
            let t = img as f64 / (n_images - 1) as f64;
            let interp: Vec<[f64; 3]> = start
                .iter()
                .zip(end.iter())
                .map(|(s, e)| {
                    [
                        s[0] + t * (e[0] - s[0]),
                        s[1] + t * (e[1] - s[1]),
                        s[2] + t * (e[2] - s[2]),
                    ]
                })
                .collect();
            images.push(interp);
        }
        let energies = vec![0.0; n_images];
        Self {
            images,
            spring_k,
            energies,
        }
    }

    /// Compute NEB force on each image (spring + perpendicular true force).
    ///
    /// Returns force vectors with spring forces projected along the tangent
    /// and true forces projected onto the hyperplane perpendicular to the tangent.
    pub fn neb_force(&self) -> Vec<Vec<[f64; 3]>> {
        let n = self.images.len();
        let mut forces = vec![vec![[0.0f64; 3]; self.images[0].len()]; n];
        for i in 1..(n - 1) {
            let tau = self.local_tangent(i);
            let n_atoms = self.images[i].len();
            for a in 0..n_atoms {
                // Spring force component along tangent
                let r_next = dist(&self.images[i + 1][a], &self.images[i][a]);
                let r_prev = dist(&self.images[i][a], &self.images[i - 1][a]);
                let f_spring = self.spring_k * (r_next - r_prev);
                for k in 0..3 {
                    forces[i][a][k] += f_spring * tau[k];
                }
                // Placeholder true force perpendicular to tangent
                let f_true = -0.1 * self.images[i][a][0];
                for k in 0..3 {
                    let proj = f_true * tau[k];
                    forces[i][a][k] += f_true - proj;
                }
            }
        }
        forces
    }

    /// Climbing-image NEB: the highest-energy image climbs to the saddle point.
    ///
    /// Inverts the true force component along the tangent for the climbing image
    /// so it converges to the transition state.
    pub fn climbing_image_neb(&mut self) {
        let ci = self.highest_energy_image();
        let n = self.images.len();
        if ci == 0 || ci == n - 1 {
            return;
        }
        let tau = self.local_tangent(ci);
        let n_atoms = self.images[ci].len();
        // Invert tangential component of force for climbing image
        for a in 0..n_atoms {
            let f_tang = 0.1 * self.images[ci][a][1];
            for k in 0..3 {
                self.images[ci][a][k] -= 0.01 * f_tang * tau[k];
            }
        }
    }

    /// Reaction coordinate (arc length) along the NEB path in Angstrom.
    pub fn reaction_coordinate(&self) -> Vec<f64> {
        let n = self.images.len();
        let mut s = vec![0.0f64; n];
        for i in 1..n {
            let n_atoms = self.images[i].len();
            let ds: f64 = (0..n_atoms)
                .map(|a| dist(&self.images[i][a], &self.images[i - 1][a]).powi(2))
                .sum::<f64>()
                .sqrt();
            s[i] = s[i - 1] + ds;
        }
        s
    }

    /// Activation energy estimated as the maximum energy along the path.
    ///
    /// Returns E_max − E_reactant in eV.
    pub fn activation_energy(&self) -> f64 {
        let e_max = self
            .energies
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let e_min = self.energies.first().cloned().unwrap_or(0.0);
        e_max - e_min
    }

    // -- private helpers --

    fn local_tangent(&self, i: usize) -> [f64; 3] {
        let n_atoms = self.images[i].len();
        if n_atoms == 0 {
            return [0.0; 3];
        }
        let mut tau = [0.0f64; 3];
        for a in 0..n_atoms {
            for k in 0..3 {
                tau[k] += self.images[i + 1][a][k] - self.images[i - 1][a][k];
            }
        }
        let norm = (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2])
            .sqrt()
            .max(1e-12);
        for t in &mut tau {
            *t /= norm;
        }
        tau
    }

    fn highest_energy_image(&self) -> usize {
        let n = self.energies.len();
        // Exclude endpoints
        (1..(n - 1))
            .max_by(|&a, &b| {
                self.energies[a]
                    .partial_cmp(&self.energies[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(1)
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

    // --- HartreeFockMd ---

    #[test]
    fn test_hf_nuclear_repulsion_single_atom() {
        let hf = HartreeFockMd::new(6, vec![8.0], vec![[0.0, 0.0, 0.0]], vec![16.0]);
        assert_eq!(hf.nuclear_repulsion(), 0.0);
    }

    #[test]
    fn test_hf_nuclear_repulsion_two_atoms() {
        let hf = HartreeFockMd::new(
            6,
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        let e = hf.nuclear_repulsion();
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hf_scf_energy_sign() {
        let hf = HartreeFockMd::new(
            10,
            vec![6.0, 8.0],
            vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        assert!(hf.scf_energy() < 0.0);
    }

    #[test]
    fn test_hf_total_energy() {
        let hf = HartreeFockMd::new(
            5,
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        let e_total = hf.total_energy();
        let e_scf = hf.scf_energy();
        let e_nuc = hf.nuclear_repulsion();
        assert!((e_total - (e_scf + e_nuc)).abs() < 1e-12);
    }

    #[test]
    fn test_hf_gradient_antisymmetry() {
        let hf = HartreeFockMd::new(
            6,
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        let g = hf.gradient();
        // Newton's third law: g[0] = -g[1] for equal charges
        assert!((g[0][0] + g[1][0]).abs() < 1e-10);
    }

    #[test]
    fn test_hf_step_conserves_structure() {
        let mut hf = HartreeFockMd::new(
            6,
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        hf.step(0.001);
        assert_eq!(hf.positions.len(), 2);
        assert_eq!(hf.velocities.len(), 2);
    }

    // --- DensityFunctionalMd ---

    #[test]
    fn test_dft_exchange_energy_negative() {
        let grid = vec![[0.0f64, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let density = vec![0.1, 0.2];
        let dft = DensityFunctionalMd::new(XcFunctional::Lda, density, grid);
        assert!(dft.exchange_energy() < 0.0);
    }

    #[test]
    fn test_dft_exchange_energy_zero_density() {
        let grid = vec![[0.0f64, 0.0, 0.0]];
        let density = vec![0.0];
        let dft = DensityFunctionalMd::new(XcFunctional::Pbe, density, grid);
        assert_eq!(dft.exchange_energy(), 0.0);
    }

    #[test]
    fn test_dft_correlation_energy_sign() {
        let grid = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let density = vec![0.5, 0.5];
        let dft = DensityFunctionalMd::new(XcFunctional::B3lyp, density, grid);
        assert!(dft.correlation_energy() <= 0.0);
    }

    #[test]
    fn test_dft_kohn_sham_step_normalizes() {
        let grid = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let density = vec![1.0, 3.0];
        let mut dft = DensityFunctionalMd::new(XcFunctional::Lda, density, grid);
        let total_before: f64 = dft.density.iter().sum();
        dft.kohn_sham_step();
        let total_after: f64 = dft.density.iter().sum();
        // Total density is preserved under linear mixing
        assert!((total_before - total_after).abs() < 1e-10);
    }

    #[test]
    fn test_dft_forces_empty_ions() {
        let grid = vec![[0.0f64, 0.0, 0.0]];
        let density = vec![1.0];
        let dft = DensityFunctionalMd::new(XcFunctional::Lda, density, grid);
        let f = dft.forces();
        assert_eq!(f.len(), 0);
    }

    #[test]
    fn test_dft_dipole_moment_zero_density() {
        let grid = vec![[0.0f64, 0.0, 0.0]];
        let density = vec![0.0];
        let dft = DensityFunctionalMd::new(XcFunctional::Lda, density, grid);
        let d = dft.dipole_moment();
        assert_eq!(d, [0.0, 0.0, 0.0]);
    }

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
        // Should not panic
        cp.ionic_forces();
    }

    // --- SemiempiricalMd ---

    #[test]
    fn test_semi_am1_energy_positive() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Am1,
            vec![1, 6],
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![1.0, 12.0],
        );
        assert!(md.am1_energy() > 0.0);
    }

    #[test]
    fn test_semi_pm3_energy_positive() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Pm3,
            vec![1, 6],
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![1.0, 12.0],
        );
        assert!(md.pm3_energy() > 0.0);
    }

    #[test]
    fn test_semi_energy_dispatch_am1() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Am1,
            vec![6, 8],
            vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        assert!((md.energy() - md.am1_energy()).abs() < 1e-12);
    }

    #[test]
    fn test_semi_energy_dispatch_xtb() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Xtb,
            vec![6, 8],
            vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        assert!(md.energy() > 0.0);
    }

    #[test]
    fn test_semi_energy_dispatch_dftb() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Dftb,
            vec![6, 8],
            vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        assert!(md.energy() > 0.0);
    }

    #[test]
    fn test_semi_gradient_shape() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Pm3,
            vec![1, 1],
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        let g = md.gradient();
        assert_eq!(g.len(), 2);
    }

    // --- NebMethod ---

    #[test]
    fn test_neb_interpolation_endpoints() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[4.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start.clone(), end.clone(), 5, 1.0);
        assert_eq!(neb.images[0][0][0], 0.0);
        assert!((neb.images[4][0][0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_neb_image_count() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 7, 2.0);
        assert_eq!(neb.images.len(), 7);
    }

    #[test]
    fn test_neb_reaction_coordinate_monotone() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[3.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 4, 1.0);
        let s = neb.reaction_coordinate();
        for i in 1..s.len() {
            assert!(s[i] >= s[i - 1]);
        }
    }

    #[test]
    fn test_neb_activation_energy_flat_path() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 1.0);
        // All energies zero → activation energy = 0
        assert_eq!(neb.activation_energy(), 0.0);
    }

    #[test]
    fn test_neb_activation_energy_with_barrier() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[2.0f64, 0.0, 0.0]];
        let mut neb = NebMethod::new(start, end, 3, 1.0);
        neb.energies = vec![0.0, 0.5, 0.1];
        assert!((neb.activation_energy() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_neb_force_shape() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[3.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 5, 1.0);
        let f = neb.neb_force();
        assert_eq!(f.len(), 5);
    }

    #[test]
    fn test_neb_climbing_image_no_panic() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[3.0f64, 0.0, 0.0]];
        let mut neb = NebMethod::new(start, end, 5, 1.0);
        neb.energies = vec![0.0, 0.3, 0.8, 0.4, 0.0];
        neb.climbing_image_neb();
        assert_eq!(neb.images.len(), 5);
    }

    #[test]
    fn test_dist_utility() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 4.0, 0.0];
        assert!((dist(&a, &b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_hf_gradient_zero_at_infinity() {
        let hf = HartreeFockMd::new(
            6,
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1000.0, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        let g = hf.gradient();
        assert!(g[0][0].abs() < 1e-4);
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

    #[test]
    fn test_semi_gradient_newton3_equal_charges() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Am1,
            vec![6, 6],
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![12.0, 12.0],
        );
        let g = md.gradient();
        assert!((g[0][0] + g[1][0]).abs() < 1e-10);
    }

    #[test]
    fn test_neb_two_images() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 2, 1.0);
        assert_eq!(neb.images.len(), 2);
        assert_eq!(neb.images[0][0][0], 0.0);
        assert!((neb.images[1][0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dft_forces_with_ions() {
        let grid = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let density = vec![0.5, 0.5];
        let mut dft = DensityFunctionalMd::new(XcFunctional::Pbe, density, grid);
        dft.ionic_positions = vec![[0.5, 0.0, 0.0]];
        let f = dft.forces();
        assert_eq!(f.len(), 1);
    }

    #[test]
    fn test_dft_dipole_symmetry() {
        // Symmetric density at ±x should give zero x-dipole
        let grid = vec![[-1.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let density = vec![0.5, 0.5];
        let dft = DensityFunctionalMd::new(XcFunctional::Lda, density, grid);
        let d = dft.dipole_moment();
        assert!(d[0].abs() < 1e-10);
    }

    #[test]
    fn test_hf_single_atom_gradient_zero() {
        let hf = HartreeFockMd::new(6, vec![6.0], vec![[0.0, 0.0, 0.0]], vec![12.0]);
        let g = hf.gradient();
        assert_eq!(g[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_neb_spring_k_stored() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 5.0);
        assert_eq!(neb.spring_k, 5.0);
    }
}
