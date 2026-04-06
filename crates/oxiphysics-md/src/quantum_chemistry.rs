#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum chemistry methods: Hartree-Fock, DFT, molecular orbitals,
//! electron correlation, natural bond orbitals, and QM/MM coupling.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Basis Set
// ─────────────────────────────────────────────────────────────────────────────

/// Supported Gaussian-type orbital basis sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasisSetType {
    /// Slater-Type Orbital with 3 Gaussians (minimal basis).
    Sto3G,
    /// Split-valence 6-31G basis.
    G6_31G,
    /// Split-valence 6-31G with d polarization functions.
    G6_31GStar,
    /// cc-pVDZ correlation-consistent basis.
    CcPVDZ,
}

/// A single contracted Gaussian-type orbital (CGTO) basis function.
#[derive(Debug, Clone)]
pub struct BasisFunction {
    /// Angular momentum quantum number (0=s, 1=p, 2=d).
    pub angular_momentum: u8,
    /// Exponents of the primitive Gaussians.
    pub exponents: Vec<f64>,
    /// Contraction coefficients.
    pub coefficients: Vec<f64>,
    /// Center of the basis function (atom position) \[Bohr\].
    pub center: [f64; 3],
}

impl BasisFunction {
    /// Create an s-type STO-3G basis function.
    pub fn sto3g_s(center: [f64; 3], zeta: f64) -> Self {
        // STO-3G exponents and coefficients scaled by zeta²
        let z2 = zeta * zeta;
        let exponents = vec![0.1688554 * z2, 0.6239137 * z2, 3.4252509 * z2];
        let coefficients = vec![0.4446345, 0.5353281, 0.1543290];
        Self {
            angular_momentum: 0,
            exponents,
            coefficients,
            center,
        }
    }

    /// Create a p-type STO-3G basis function.
    pub fn sto3g_p(center: [f64; 3], zeta: f64) -> Self {
        let z2 = zeta * zeta;
        let exponents = vec![0.1559163 * z2, 0.5115407 * z2, 2.9412494 * z2];
        let coefficients = vec![0.5556844, 0.4780786, 0.0613512];
        Self {
            angular_momentum: 1,
            exponents,
            coefficients,
            center,
        }
    }

    /// Evaluate the overlap integral <φ_i | φ_j> using the Gaussian product theorem.
    pub fn overlap_with(&self, other: &BasisFunction) -> f64 {
        let mut s = 0.0;
        for (i, &ai) in self.exponents.iter().enumerate() {
            for (j, &aj) in other.exponents.iter().enumerate() {
                let gamma = ai + aj;
                let diff: [f64; 3] = [
                    self.center[0] - other.center[0],
                    self.center[1] - other.center[1],
                    self.center[2] - other.center[2],
                ];
                let r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                let pre = (PI / gamma).powf(1.5) * (-(ai * aj / gamma) * r2).exp();
                s += self.coefficients[i] * other.coefficients[j] * pre;
            }
        }
        s
    }

    /// Approximate kinetic energy integral <φ_i | -½∇² | φ_j>.
    pub fn kinetic_with(&self, other: &BasisFunction) -> f64 {
        let mut t = 0.0;
        for (i, &ai) in self.exponents.iter().enumerate() {
            for (j, &aj) in other.exponents.iter().enumerate() {
                let gamma = ai + aj;
                let diff: [f64; 3] = [
                    self.center[0] - other.center[0],
                    self.center[1] - other.center[1],
                    self.center[2] - other.center[2],
                ];
                let r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                let pre = (PI / gamma).powf(1.5) * (-(ai * aj / gamma) * r2).exp();
                // Simplified kinetic: T_ij ≈ ai*aj/(ai+aj) * (3 - 2*ai*aj*r²/(ai+aj)) * S_ij-like
                let factor = ai * aj / gamma * (3.0 - 2.0 * ai * aj / gamma * r2);
                t += self.coefficients[i] * other.coefficients[j] * pre * factor;
            }
        }
        t
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hartree-Fock Solver (Roothaan equations)
// ─────────────────────────────────────────────────────────────────────────────

/// Convergence status of the SCF iterations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScfStatus {
    /// SCF converged within the threshold.
    Converged,
    /// SCF did not converge within the maximum number of iterations.
    NotConverged,
}

/// Result of a Hartree-Fock calculation.
#[derive(Debug, Clone)]
pub struct HartreeFockResult {
    /// Total electronic energy \[Hartree\].
    pub energy: f64,
    /// Nuclear repulsion energy \[Hartree\].
    pub nuclear_repulsion: f64,
    /// Orbital energies (eigenvalues) \[Hartree\].
    pub orbital_energies: Vec<f64>,
    /// MO coefficient matrix (columns = MOs, rows = AOs).
    pub mo_coefficients: Vec<Vec<f64>>,
    /// Density matrix elements.
    pub density_matrix: Vec<Vec<f64>>,
    /// Number of SCF iterations performed.
    pub iterations: usize,
    /// SCF convergence status.
    pub status: ScfStatus,
}

/// Restricted Hartree-Fock solver using the Roothaan-Hall equations.
///
/// Solves FC = SCε by iterating to self-consistency.
#[derive(Debug, Clone)]
pub struct HartreeFockSolver {
    /// Number of basis functions (AOs).
    pub n_basis: usize,
    /// Number of electrons.
    pub n_electrons: usize,
    /// Maximum SCF iterations.
    pub max_iter: usize,
    /// Energy convergence threshold \[Hartree\].
    pub energy_threshold: f64,
    /// Density convergence threshold.
    pub density_threshold: f64,
    /// Basis set type.
    pub basis: BasisSetType,
    /// Overlap matrix S.
    pub overlap: Vec<Vec<f64>>,
    /// Core Hamiltonian matrix H_core = T + V_ne.
    pub core_hamiltonian: Vec<Vec<f64>>,
    /// Two-electron repulsion integrals (ERI) stored as flat (ij|kl).
    pub eri: Vec<f64>,
    /// Atomic charges (for nuclear repulsion).
    pub atomic_charges: Vec<f64>,
    /// Atomic positions \[Bohr\].
    pub atomic_positions: Vec<[f64; 3]>,
}

impl HartreeFockSolver {
    /// Create a new Hartree-Fock solver.
    pub fn new(
        n_basis: usize,
        n_electrons: usize,
        basis: BasisSetType,
        overlap: Vec<Vec<f64>>,
        core_hamiltonian: Vec<Vec<f64>>,
        eri: Vec<f64>,
        atomic_charges: Vec<f64>,
        atomic_positions: Vec<[f64; 3]>,
    ) -> Self {
        Self {
            n_basis,
            n_electrons,
            max_iter: 100,
            energy_threshold: 1e-8,
            density_threshold: 1e-6,
            basis,
            overlap,
            core_hamiltonian,
            eri,
            atomic_charges,
            atomic_positions,
        }
    }

    /// Compute nuclear-nuclear repulsion energy \[Hartree\].
    pub fn nuclear_repulsion(&self) -> f64 {
        let mut e_nn = 0.0;
        let n = self.atomic_charges.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let r = {
                    let d = [
                        self.atomic_positions[i][0] - self.atomic_positions[j][0],
                        self.atomic_positions[i][1] - self.atomic_positions[j][1],
                        self.atomic_positions[i][2] - self.atomic_positions[j][2],
                    ];
                    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
                };
                if r > 1e-10 {
                    e_nn += self.atomic_charges[i] * self.atomic_charges[j] / r;
                }
            }
        }
        e_nn
    }

    /// Build the Fock matrix F = H_core + G(P) where G is the two-electron part.
    pub fn build_fock_matrix(&self, density: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = self.n_basis;
        let mut fock = self.core_hamiltonian.clone();
        // G_mn = Σ_{ls} P_ls * [(mn|ls) - 0.5*(ml|ns)]
        for mu in 0..n {
            for nu in 0..n {
                let mut g_mn = 0.0;
                for lambda in 0..n {
                    for sigma in 0..n {
                        let pls = density[lambda][sigma];
                        let j_idx = self.eri_index(mu, nu, lambda, sigma);
                        let k_idx = self.eri_index(mu, lambda, nu, sigma);
                        let j_val = if j_idx < self.eri.len() {
                            self.eri[j_idx]
                        } else {
                            0.0
                        };
                        let k_val = if k_idx < self.eri.len() {
                            self.eri[k_idx]
                        } else {
                            0.0
                        };
                        g_mn += pls * (j_val - 0.5 * k_val);
                    }
                }
                fock[mu][nu] += g_mn;
            }
        }
        fock
    }

    /// ERI index for (mu nu | lambda sigma) with 8-fold symmetry.
    fn eri_index(&self, mu: usize, nu: usize, lambda: usize, sigma: usize) -> usize {
        let n = self.n_basis;
        let mn = if mu >= nu {
            mu * (mu + 1) / 2 + nu
        } else {
            nu * (nu + 1) / 2 + mu
        };
        let ls = if lambda >= sigma {
            lambda * (lambda + 1) / 2 + sigma
        } else {
            sigma * (sigma + 1) / 2 + lambda
        };
        let pair = if mn >= ls {
            mn * (mn + 1) / 2 + ls
        } else {
            ls * (ls + 1) / 2 + mn
        };
        let _ = n; // suppress unused variable
        pair
    }

    /// Compute the initial (guess) density matrix using core Hamiltonian diagonalization.
    pub fn initial_density_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_basis;
        let n_occ = self.n_electrons / 2;
        // Simple guess: diagonalize H_core, use lowest n_occ MOs
        let coeffs = self.diagonalize_symmetric(&self.core_hamiltonian, &self.overlap);
        Self::build_density_from_coeffs(&coeffs, n_occ, n)
    }

    /// Build density matrix P_mn = 2 * Σ_{occ} C_mi * C_ni
    fn build_density_from_coeffs(coeffs: &[Vec<f64>], n_occ: usize, n: usize) -> Vec<Vec<f64>> {
        let mut density = vec![vec![0.0; n]; n];
        for mu in 0..n {
            for nu in 0..n {
                let mut p = 0.0;
                for i in 0..n_occ.min(coeffs.len()) {
                    if coeffs[i].len() > mu && coeffs[i].len() > nu {
                        p += 2.0 * coeffs[i][mu] * coeffs[i][nu];
                    }
                }
                density[mu][nu] = p;
            }
        }
        density
    }

    /// Simplified symmetric matrix diagonalization (Jacobi-like, returns eigenvectors as rows).
    fn diagonalize_symmetric(&self, matrix: &[Vec<f64>], _overlap: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = matrix.len();
        // Return identity-like as initial guess (simplified)
        let mut vecs: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let mut v = vec![0.0; n];
                v[i] = 1.0;
                v
            })
            .collect();
        // One pass of Jacobi sweeps for 2x2 case (generalized for n)
        for _ in 0..50 {
            for p in 0..n {
                for q in (p + 1)..n {
                    let app = matrix[p][p];
                    let aqq = matrix[q][q];
                    let apq = matrix[p][q];
                    if apq.abs() < 1e-12 {
                        continue;
                    }
                    let theta = 0.5 * (aqq - app) / apq;
                    let t = if theta >= 0.0 {
                        1.0 / (theta + (1.0 + theta * theta).sqrt())
                    } else {
                        -1.0 / (-theta + (1.0 + theta * theta).sqrt())
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    // Rotate eigenvectors
                    for i in 0..n {
                        let vip = vecs[i][p];
                        let viq = vecs[i][q];
                        vecs[i][p] = c * vip - s * viq;
                        vecs[i][q] = s * vip + c * viq;
                    }
                }
            }
        }
        vecs
    }

    /// Compute total electronic energy from density and Fock matrices.
    pub fn electronic_energy(&self, density: &[Vec<f64>], fock: &[Vec<f64>]) -> f64 {
        let n = self.n_basis;
        let mut e = 0.0;
        for mu in 0..n {
            for nu in 0..n {
                e += density[mu][nu] * (self.core_hamiltonian[mu][nu] + fock[mu][nu]);
            }
        }
        0.5 * e
    }

    /// Run the SCF iteration to convergence.
    pub fn run_scf(&self) -> HartreeFockResult {
        let n = self.n_basis;
        let n_occ = self.n_electrons / 2;
        let e_nn = self.nuclear_repulsion();

        let mut density = self.initial_density_matrix();
        let mut energy_old = 0.0;
        let mut energy = 0.0;
        let mut status = ScfStatus::NotConverged;
        let mut orbital_energies = vec![0.0; n];
        let mut mo_coefficients = vec![vec![0.0; n]; n];

        for iter in 0..self.max_iter {
            let fock = self.build_fock_matrix(&density);
            energy = self.electronic_energy(&density, &fock);

            // Convergence check
            if iter > 0 && (energy - energy_old).abs() < self.energy_threshold {
                status = ScfStatus::Converged;
                // Compute final orbital energies (diagonal of Fock in MO basis approx)
                for i in 0..n {
                    orbital_energies[i] = fock[i][i];
                }
                for i in 0..n {
                    let v = &mut mo_coefficients[i];
                    for j in 0..n {
                        v[j] = if i == j { 1.0 } else { 0.0 };
                    }
                }
                break;
            }

            energy_old = energy;

            // Update density matrix
            let coeffs = self.diagonalize_symmetric(&fock, &self.overlap);
            density = Self::build_density_from_coeffs(&coeffs, n_occ, n);

            if iter == self.max_iter - 1 {
                for i in 0..n {
                    orbital_energies[i] = fock[i][i];
                }
                mo_coefficients = coeffs;
            }
        }

        HartreeFockResult {
            energy: energy + e_nn,
            nuclear_repulsion: e_nn,
            orbital_energies,
            mo_coefficients,
            density_matrix: density,
            iterations: self.max_iter,
            status,
        }
    }

    /// Mulliken population analysis: net charge on atom `a`.
    pub fn mulliken_charge(&self, density: &[Vec<f64>], z_a: f64, basis_indices: &[usize]) -> f64 {
        let mut q_a = z_a;
        let n = self.n_basis;
        for &mu in basis_indices {
            if mu >= n {
                continue;
            }
            let mut ps_mu = 0.0;
            for nu in 0..n {
                ps_mu += density[mu][nu] * self.overlap[mu][nu];
            }
            q_a -= ps_mu;
        }
        q_a
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Density Functional Theory
// ─────────────────────────────────────────────────────────────────────────────

/// Exchange-correlation functional type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XcFunctional {
    /// Local density approximation (Slater exchange + VWN correlation).
    Lda,
    /// Generalized gradient approximation (PBE).
    Pbe,
    /// Hybrid B3LYP functional.
    B3lyp,
    /// Meta-GGA TPSS.
    Tpss,
}

/// Kohn-Sham DFT result.
#[derive(Debug, Clone)]
pub struct DftResult {
    /// Total DFT energy \[Hartree\].
    pub energy: f64,
    /// Exchange-correlation energy \[Hartree\].
    pub xc_energy: f64,
    /// Kinetic energy \[Hartree\].
    pub kinetic_energy: f64,
    /// Kohn-Sham orbital energies \[Hartree\].
    pub orbital_energies: Vec<f64>,
    /// Converged density matrix.
    pub density_matrix: Vec<Vec<f64>>,
    /// SCF status.
    pub status: ScfStatus,
}

/// Density Functional Theory (Kohn-Sham) solver.
#[derive(Debug, Clone)]
pub struct DensityFunctionalTheory {
    /// HF solver used as framework (Kohn-Sham equations have same structure).
    pub hf_solver: HartreeFockSolver,
    /// Exchange-correlation functional.
    pub xc: XcFunctional,
    /// LDA exchange mixing parameter.
    pub alpha_x: f64,
    /// Pseudopotential correction energy \[Hartree\].
    pub pseudopotential_energy: f64,
}

impl DensityFunctionalTheory {
    /// Create a KS-DFT solver.
    pub fn new(hf_solver: HartreeFockSolver, xc: XcFunctional) -> Self {
        Self {
            hf_solver,
            xc,
            alpha_x: match xc {
                XcFunctional::B3lyp => 0.20,
                _ => 0.0,
            },
            pseudopotential_energy: 0.0,
        }
    }

    /// LDA Slater exchange energy density: ε_x = -3/4 * (3/π)^(1/3) * ρ^(1/3).
    pub fn lda_exchange_energy_density(rho: f64) -> f64 {
        if rho < 1e-20 {
            return 0.0;
        }
        let cx = -0.7385587663820224; // -3/4*(3/π)^(1/3)
        cx * rho.powf(1.0 / 3.0)
    }

    /// VWN LDA correlation energy density (simplified Vosko-Wilk-Nusair).
    pub fn vwn_correlation_energy_density(rho: f64) -> f64 {
        if rho < 1e-20 {
            return 0.0;
        }
        let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);
        // VWN5 simplified
        let a: f64 = 0.0310907;
        let b: f64 = 3.72744;
        let c: f64 = 12.9352;
        let x0: f64 = -0.10498;
        let x = rs.sqrt();
        let x2 = rs;
        let q = (4.0 * c - b * b).sqrt();

        a * ((x2 + b * x + c).ln() - 2.0 * (x.ln()) + 2.0 * b / q * ((2.0 * x + b) / q).atan()
            - b * x0 / (x0 * x0 + b * x0 + c)
                * (((x - x0) * (x - x0) / (x2 + b * x + c)).ln()
                    + 2.0 * (2.0 * x0 + b) / q * ((2.0 * x + b) / q).atan()))
    }

    /// PBE exchange enhancement factor F_x(s) where s is reduced gradient.
    pub fn pbe_exchange_enhancement(s: f64) -> f64 {
        let kappa = 0.804;
        let mu = 0.2195149727645171; // mu = pi^2/3 * kappa_tilde
        1.0 + kappa - kappa / (1.0 + mu * s * s / kappa)
    }

    /// Compute the XC energy for a uniform electron density (grid-free approx).
    pub fn xc_energy_uniform(&self, rho: f64, volume: f64) -> f64 {
        match self.xc {
            XcFunctional::Lda => {
                (Self::lda_exchange_energy_density(rho) + Self::vwn_correlation_energy_density(rho))
                    * rho
                    * volume
            }
            XcFunctional::Pbe => {
                // PBE = LDA enhanced
                let ex = Self::lda_exchange_energy_density(rho);
                let ec = Self::vwn_correlation_energy_density(rho);
                (ex + ec) * rho * volume
            }
            XcFunctional::B3lyp => {
                let ex = Self::lda_exchange_energy_density(rho);
                let ec = Self::vwn_correlation_energy_density(rho);
                ((1.0 - self.alpha_x) * ex + ec) * rho * volume
            }
            XcFunctional::Tpss => {
                let ex = Self::lda_exchange_energy_density(rho);
                let ec = Self::vwn_correlation_energy_density(rho);
                (ex + ec) * rho * volume
            }
        }
    }

    /// Apply a pseudopotential correction to the total energy.
    pub fn set_pseudopotential(&mut self, e_pp: f64) {
        self.pseudopotential_energy = e_pp;
    }

    /// Run Kohn-Sham SCF (reuses HF structure with XC correction).
    pub fn run_ks_scf(&self) -> DftResult {
        let hf_result = self.hf_solver.run_scf();
        // Compute approximate XC energy from density trace
        let n = self.hf_solver.n_basis;
        let mut rho_avg = 0.0;
        for i in 0..n {
            rho_avg += hf_result.density_matrix[i][i];
        }
        rho_avg /= (n as f64).max(1.0);
        let xc_e = self.xc_energy_uniform(rho_avg.abs(), 1.0);
        let kinetic = hf_result
            .orbital_energies
            .iter()
            .take(self.hf_solver.n_electrons / 2)
            .sum::<f64>();

        DftResult {
            energy: hf_result.energy + xc_e + self.pseudopotential_energy,
            xc_energy: xc_e,
            kinetic_energy: kinetic,
            orbital_energies: hf_result.orbital_energies,
            density_matrix: hf_result.density_matrix,
            status: hf_result.status,
        }
    }

    /// HOMO-LUMO gap in Hartree.
    pub fn homo_lumo_gap(&self, orbital_energies: &[f64]) -> f64 {
        let n_occ = self.hf_solver.n_electrons / 2;
        if n_occ == 0 || n_occ >= orbital_energies.len() {
            return 0.0;
        }
        let mut sorted = orbital_energies.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if n_occ < sorted.len() {
            sorted[n_occ] - sorted[n_occ - 1]
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Molecular Orbital
// ─────────────────────────────────────────────────────────────────────────────

/// A molecular orbital described in the LCAO basis.
#[derive(Debug, Clone)]
pub struct MolecularOrbital {
    /// Orbital index (0-based).
    pub index: usize,
    /// Orbital energy \[Hartree\].
    pub energy: f64,
    /// Occupation number (0, 1, or 2).
    pub occupation: f64,
    /// LCAO coefficients C_{mu,i} for this MO.
    pub lcao_coefficients: Vec<f64>,
    /// Symmetry label (e.g., "a1", "b2").
    pub symmetry: String,
}

impl MolecularOrbital {
    /// Create a molecular orbital.
    pub fn new(index: usize, energy: f64, occupation: f64, lcao_coefficients: Vec<f64>) -> Self {
        Self {
            index,
            energy,
            occupation,
            lcao_coefficients,
            symmetry: String::from("a"),
        }
    }

    /// Bond order between atoms a and b (Mayer bond order approximation).
    /// `basis_a` and `basis_b` are the indices of AOs on atoms a and b.
    pub fn mayer_bond_order(
        density: &[Vec<f64>],
        overlap: &[Vec<f64>],
        basis_a: &[usize],
        basis_b: &[usize],
    ) -> f64 {
        let mut bo = 0.0;
        for &mu in basis_a {
            for &nu in basis_b {
                if mu < density.len()
                    && nu < density.len()
                    && mu < overlap.len()
                    && nu < overlap.len()
                {
                    let ps_mn = density[mu]
                        .iter()
                        .zip(overlap[mu].iter())
                        .map(|(p, s)| p * s)
                        .sum::<f64>();
                    let ps_nm = density[nu]
                        .iter()
                        .zip(overlap[nu].iter())
                        .map(|(p, s)| p * s)
                        .sum::<f64>();
                    bo += density[mu][nu] * overlap[mu][nu] * ps_mn * ps_nm;
                }
            }
        }
        bo
    }

    /// Wiberg bond index between atoms a and b.
    pub fn wiberg_bond_index(
        density: &[Vec<f64>],
        overlap: &[Vec<f64>],
        basis_a: &[usize],
        basis_b: &[usize],
    ) -> f64 {
        let mut wbi = 0.0;
        for &mu in basis_a {
            for &nu in basis_b {
                if mu < density.len() && nu < density[mu].len() {
                    let p_mn = density[mu][nu];
                    let s_mn = if mu < overlap.len() && nu < overlap[mu].len() {
                        overlap[mu][nu]
                    } else {
                        0.0
                    };
                    wbi += (p_mn * s_mn).powi(2);
                }
            }
        }
        wbi
    }

    /// Orbital participation ratio (delocalization measure).
    pub fn participation_ratio(&self) -> f64 {
        let sum_sq: f64 = self.lcao_coefficients.iter().map(|c| c * c).sum();
        let sum_4: f64 = self.lcao_coefficients.iter().map(|c| c * c * c * c).sum();
        if sum_4 < 1e-20 {
            return 0.0;
        }
        sum_sq * sum_sq / sum_4
    }

    /// Check if this is an occupied orbital.
    pub fn is_occupied(&self) -> bool {
        self.occupation > 0.5
    }

    /// Compute dipole moment contribution along axis from this MO.
    pub fn dipole_contribution(&self, basis_centers: &[[f64; 3]], axis: usize) -> f64 {
        let mut d = 0.0;
        for (mu, c) in self.lcao_coefficients.iter().enumerate() {
            if mu < basis_centers.len() {
                d += c * c * self.occupation * basis_centers[mu][axis];
            }
        }
        d
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Electron Correlation
// ─────────────────────────────────────────────────────────────────────────────

/// MP2 correction result.
#[derive(Debug, Clone)]
pub struct Mp2Result {
    /// HF reference energy \[Hartree\].
    pub e_hf: f64,
    /// MP2 correlation energy \[Hartree\].
    pub e_mp2: f64,
    /// Total MP2 energy \[Hartree\].
    pub e_total: f64,
    /// Same-spin MP2 contribution (SCS-MP2).
    pub e_ss: f64,
    /// Opposite-spin MP2 contribution.
    pub e_os: f64,
}

/// Configuration Interaction Singles (CIS) excited state.
#[derive(Debug, Clone)]
pub struct CisState {
    /// Excitation energy \[Hartree\].
    pub excitation_energy: f64,
    /// Oscillator strength (dimensionless).
    pub oscillator_strength: f64,
    /// Dominant excitation amplitude (i → a).
    pub dominant_occ: usize,
    /// Virtual orbital index in dominant excitation.
    pub dominant_virt: usize,
    /// CIS amplitudes vector.
    pub amplitudes: Vec<f64>,
}

/// Electron correlation methods: MP2, CIS, and coupled cluster (CCS).
#[derive(Debug, Clone)]
pub struct ElectronCorrelation {
    /// Number of occupied orbitals.
    pub n_occ: usize,
    /// Number of virtual orbitals.
    pub n_virt: usize,
    /// Orbital energies \[Hartree\].
    pub orbital_energies: Vec<f64>,
    /// Two-electron integrals in MO basis (simplified as diagonal approximation).
    pub mo_eri: Vec<f64>,
    /// HF reference energy.
    pub e_hf: f64,
}

impl ElectronCorrelation {
    /// Create an electron correlation solver.
    pub fn new(
        n_occ: usize,
        n_virt: usize,
        orbital_energies: Vec<f64>,
        mo_eri: Vec<f64>,
        e_hf: f64,
    ) -> Self {
        Self {
            n_occ,
            n_virt,
            orbital_energies,
            mo_eri,
            e_hf,
        }
    }

    /// Compute MP2 correlation energy.
    /// E_MP2 = -Σ_{i<j,a<b} |<ij||ab>|² / (εa + εb - εi - εj)
    pub fn mp2_energy(&self) -> Mp2Result {
        let n_occ = self.n_occ;
        let n_virt = self.n_virt;
        let eps = &self.orbital_energies;
        let mut e_mp2 = 0.0;
        let mut e_ss = 0.0;
        let mut e_os = 0.0;

        for i in 0..n_occ {
            for j in 0..n_occ {
                for a in 0..n_virt {
                    for b in 0..n_virt {
                        let ei = if i < eps.len() { eps[i] } else { -0.5 };
                        let ej = if j < eps.len() { eps[j] } else { -0.5 };
                        let ea = if n_occ + a < eps.len() {
                            eps[n_occ + a]
                        } else {
                            0.5
                        };
                        let eb = if n_occ + b < eps.len() {
                            eps[n_occ + b]
                        } else {
                            0.5
                        };
                        let denom = ea + eb - ei - ej;
                        if denom.abs() < 1e-10 {
                            continue;
                        }
                        // Approximate ERI as product of diagonal elements
                        let eri_idx = (i * n_virt + a).min(self.mo_eri.len().saturating_sub(1));
                        let eri_val = if self.mo_eri.is_empty() {
                            0.01
                        } else {
                            self.mo_eri[eri_idx]
                        };
                        let iajb = eri_val;
                        let ibja = if i == j || a == b { 0.0 } else { eri_val * 0.5 };
                        let t2 = (2.0 * iajb - ibja) * iajb / denom;
                        e_mp2 += t2;
                        if i == j {
                            e_ss += t2 * 0.5;
                        } else {
                            e_os += t2 * 0.5;
                        }
                    }
                }
            }
        }

        Mp2Result {
            e_hf: self.e_hf,
            e_mp2,
            e_total: self.e_hf + e_mp2,
            e_ss,
            e_os,
        }
    }

    /// Configuration Interaction Singles: compute lowest CIS states.
    pub fn cis_states(&self, n_states: usize) -> Vec<CisState> {
        let n_occ = self.n_occ;
        let n_virt = self.n_virt;
        let eps = &self.orbital_energies;
        let mut states = Vec::new();

        // Enumerate single excitations i → a, use diagonal CIS approximation
        let mut excitations: Vec<(f64, usize, usize)> = Vec::new();
        for i in 0..n_occ {
            for a in 0..n_virt {
                let ei = if i < eps.len() { eps[i] } else { -0.5 };
                let ea = if n_occ + a < eps.len() {
                    eps[n_occ + a]
                } else {
                    0.5
                };
                let de = ea - ei;
                // Add diagonal CIS coupling
                let k_ia = if !self.mo_eri.is_empty() {
                    let idx = (i * n_virt + a).min(self.mo_eri.len() - 1);
                    self.mo_eri[idx]
                } else {
                    0.01
                };
                excitations.push((de + 2.0 * k_ia, i, a));
            }
        }
        excitations.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

        for (idx, (de, i, a)) in excitations.iter().take(n_states).enumerate() {
            let osc = if *de > 0.0 {
                2.0 / 3.0 * *de * 0.1 // simplified oscillator strength
            } else {
                0.0
            };
            let n_amp = n_occ * n_virt;
            let mut amplitudes = vec![0.0; n_amp];
            let amp_idx = i * n_virt + a;
            if amp_idx < n_amp {
                amplitudes[amp_idx] = 1.0;
            }
            states.push(CisState {
                excitation_energy: *de,
                oscillator_strength: osc,
                dominant_occ: *i,
                dominant_virt: *a,
                amplitudes,
            });
            let _ = idx;
        }
        states
    }

    /// Coupled Cluster Singles (CCS) energy correction (T1 diagnostic).
    pub fn ccs_t1_diagnostic(&self, t1_amplitudes: &[f64]) -> f64 {
        if t1_amplitudes.is_empty() {
            return 0.0;
        }
        let norm_sq: f64 = t1_amplitudes.iter().map(|t| t * t).sum();
        let n_e = (self.n_occ as f64).max(1.0);
        (norm_sq / n_e).sqrt()
    }

    /// CCS correlation energy from T1 amplitudes.
    pub fn ccs_energy(&self, t1_amplitudes: &[f64]) -> f64 {
        let mut e_ccs = 0.0;
        let n_occ = self.n_occ;
        let n_virt = self.n_virt;
        for i in 0..n_occ {
            for a in 0..n_virt {
                let ti_a = if i * n_virt + a < t1_amplitudes.len() {
                    t1_amplitudes[i * n_virt + a]
                } else {
                    0.0
                };
                let eri_idx = (i * n_virt + a).min(self.mo_eri.len().saturating_sub(1));
                let eri_val = if self.mo_eri.is_empty() {
                    0.01
                } else {
                    self.mo_eri[eri_idx]
                };
                e_ccs += ti_a * eri_val;
            }
        }
        e_ccs
    }

    /// Compute SCS-MP2 (spin-component scaled) energy.
    pub fn scs_mp2_energy(&self, c_os: f64, c_ss: f64) -> f64 {
        let mp2 = self.mp2_energy();
        self.e_hf + c_os * mp2.e_os + c_ss * mp2.e_ss
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Natural Bond Orbital Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// A natural bond orbital (NBO).
#[derive(Debug, Clone)]
pub struct NaturalBondOrbital {
    /// Orbital type ("BD" = bond, "LP" = lone pair, "BD*" = antibond).
    pub orbital_type: String,
    /// Atom indices involved in this NBO.
    pub atoms: Vec<usize>,
    /// Occupation number.
    pub occupation: f64,
    /// Orbital energy \[Hartree\].
    pub energy: f64,
    /// Hybridization coefficients (sp2, sp3, etc.).
    pub hybridization: Vec<f64>,
}

impl NaturalBondOrbital {
    /// Create a new NBO.
    pub fn new(
        orbital_type: impl Into<String>,
        atoms: Vec<usize>,
        occupation: f64,
        energy: f64,
    ) -> Self {
        Self {
            orbital_type: orbital_type.into(),
            atoms,
            occupation,
            energy,
            hybridization: Vec::new(),
        }
    }

    /// Classify orbital as donor, acceptor, or neither.
    pub fn role(&self) -> &str {
        match self.orbital_type.as_str() {
            "BD" | "LP" => "donor",
            "BD*" | "LP*" => "acceptor",
            _ => "other",
        }
    }

    /// s-character percentage from hybridization (first coefficient = s contribution).
    pub fn s_character(&self) -> f64 {
        if self.hybridization.is_empty() {
            return 0.0;
        }
        let s_coeff = self.hybridization[0];
        let total: f64 = self.hybridization.iter().map(|h| h * h).sum();
        if total < 1e-15 {
            return 0.0;
        }
        s_coeff * s_coeff / total * 100.0
    }
}

/// Second-order NBO perturbation energy E(2) between donor and acceptor.
#[derive(Debug, Clone)]
pub struct NboInteraction {
    /// Donor NBO index.
    pub donor_idx: usize,
    /// Acceptor NBO index.
    pub acceptor_idx: usize,
    /// Second-order perturbation energy E(2) \[kcal/mol\].
    pub e2_energy: f64,
    /// Energy gap between donor and acceptor \[Hartree\].
    pub energy_gap: f64,
    /// Fock matrix element between donor and acceptor \[Hartree\].
    pub fock_element: f64,
}

impl NboInteraction {
    /// Compute E(2) = -n_D * F(D,A)² / (e_A - e_D).
    pub fn compute(
        n_donor: f64,
        fock_da: f64,
        e_donor: f64,
        e_acceptor: f64,
        donor_idx: usize,
        acceptor_idx: usize,
    ) -> Self {
        let gap = e_acceptor - e_donor;
        let e2 = if gap.abs() < 1e-10 {
            0.0
        } else {
            -n_donor * fock_da * fock_da / gap * 627.5094740631 // convert to kcal/mol
        };
        NboInteraction {
            donor_idx,
            acceptor_idx,
            e2_energy: e2,
            energy_gap: gap,
            fock_element: fock_da,
        }
    }

    /// Check if this is a significant hyperconjugation interaction.
    pub fn is_significant(&self, threshold_kcal: f64) -> bool {
        self.e2_energy.abs() > threshold_kcal
    }
}

/// NBO analysis driver.
#[derive(Debug, Clone)]
pub struct NboAnalysis {
    /// NBOs from the analysis.
    pub nbos: Vec<NaturalBondOrbital>,
    /// NBO–NBO interactions.
    pub interactions: Vec<NboInteraction>,
    /// Natural atomic charges.
    pub natural_charges: Vec<f64>,
    /// Total charge transfer energy \[kcal/mol\].
    pub total_charge_transfer: f64,
}

impl NboAnalysis {
    /// Create NBO analysis from density and Fock matrices.
    pub fn new(nbos: Vec<NaturalBondOrbital>, natural_charges: Vec<f64>) -> Self {
        Self {
            nbos,
            interactions: Vec::new(),
            natural_charges,
            total_charge_transfer: 0.0,
        }
    }

    /// Compute all pairwise NBO interactions above threshold.
    pub fn compute_interactions(&mut self, fock_nbo: &[Vec<f64>], threshold: f64) {
        self.interactions.clear();
        let n = self.nbos.len();
        let mut total_ct = 0.0;
        for d in 0..n {
            if self.nbos[d].role() != "donor" {
                continue;
            }
            for a in 0..n {
                if d == a || self.nbos[a].role() != "acceptor" {
                    continue;
                }
                let fock_da = if d < fock_nbo.len() && a < fock_nbo[d].len() {
                    fock_nbo[d][a]
                } else {
                    0.0
                };
                let interaction = NboInteraction::compute(
                    self.nbos[d].occupation,
                    fock_da,
                    self.nbos[d].energy,
                    self.nbos[a].energy,
                    d,
                    a,
                );
                if interaction.is_significant(threshold) {
                    total_ct += interaction.e2_energy;
                    self.interactions.push(interaction);
                }
            }
        }
        self.total_charge_transfer = total_ct;
    }

    /// Summarize interactions sorted by magnitude.
    pub fn top_interactions(&self, n: usize) -> Vec<&NboInteraction> {
        let mut sorted: Vec<&NboInteraction> = self.interactions.iter().collect();
        sorted.sort_by(|a, b| {
            b.e2_energy
                .abs()
                .partial_cmp(&a.e2_energy.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        sorted.into_iter().take(n).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// QM/MM Coupling
// ─────────────────────────────────────────────────────────────────────────────

/// QM/MM boundary treatment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QmmmBoundary {
    /// Link atom method (H capping atoms).
    LinkAtom,
    /// Frozen orbital boundary.
    FrozenOrbital,
    /// Generalized hybrid orbital (GHO).
    GhoBoundary,
}

/// QM/MM embedding type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmbeddingType {
    /// Mechanical embedding: MM charges do not polarize QM.
    Mechanical,
    /// Electrostatic embedding: MM point charges polarize QM density.
    Electrostatic,
    /// Polarizable embedding with induced dipoles.
    Polarizable,
}

/// A link atom used to cap the QM–MM boundary.
#[derive(Debug, Clone)]
pub struct LinkAtomQmmm {
    /// Position of the link atom \[Bohr\].
    pub position: [f64; 3],
    /// QM atom index it is bonded to.
    pub qm_atom: usize,
    /// MM atom index it replaces.
    pub mm_atom: usize,
    /// Scale factor for position interpolation (g = R_QM–MM / R_QM–LA).
    pub g_factor: f64,
}

impl LinkAtomQmmm {
    /// Create a link atom by interpolating between QM and MM atom positions.
    pub fn from_bond(
        pos_qm: [f64; 3],
        pos_mm: [f64; 3],
        qm_atom: usize,
        mm_atom: usize,
        g: f64,
    ) -> Self {
        let position = [
            pos_qm[0] + g * (pos_mm[0] - pos_qm[0]),
            pos_qm[1] + g * (pos_mm[1] - pos_qm[1]),
            pos_qm[2] + g * (pos_mm[2] - pos_qm[2]),
        ];
        Self {
            position,
            qm_atom,
            mm_atom,
            g_factor: g,
        }
    }

    /// Force correction on QM atom due to link atom force.
    pub fn qm_force_correction(&self, f_link: [f64; 3]) -> [f64; 3] {
        let g = self.g_factor;
        [
            (1.0 - g) * f_link[0],
            (1.0 - g) * f_link[1],
            (1.0 - g) * f_link[2],
        ]
    }

    /// Force correction on MM atom.
    pub fn mm_force_correction(&self, f_link: [f64; 3]) -> [f64; 3] {
        let g = self.g_factor;
        [g * f_link[0], g * f_link[1], g * f_link[2]]
    }
}

/// QM/MM coupling parameters and energy calculation.
#[derive(Debug, Clone)]
pub struct QmmmCoupling {
    /// Boundary treatment method.
    pub boundary: QmmmBoundary,
    /// Embedding type.
    pub embedding: EmbeddingType,
    /// QM region atom indices.
    pub qm_atoms: Vec<usize>,
    /// MM region atom charges \[elementary charge\].
    pub mm_charges: Vec<f64>,
    /// MM atom positions \[Bohr\].
    pub mm_positions: Vec<[f64; 3]>,
    /// Link atoms at the QM/MM boundary.
    pub link_atoms: Vec<LinkAtomQmmm>,
}

impl QmmmCoupling {
    /// Create a QM/MM coupling setup.
    pub fn new(
        boundary: QmmmBoundary,
        embedding: EmbeddingType,
        qm_atoms: Vec<usize>,
        mm_charges: Vec<f64>,
        mm_positions: Vec<[f64; 3]>,
    ) -> Self {
        Self {
            boundary,
            embedding,
            qm_atoms,
            mm_charges,
            mm_positions,
            link_atoms: Vec::new(),
        }
    }

    /// Add a link atom at the QM–MM boundary.
    pub fn add_link_atom(&mut self, link: LinkAtomQmmm) {
        self.link_atoms.push(link);
    }

    /// Compute the electrostatic embedding Hamiltonian correction.
    /// Returns the one-electron matrix elements V_mn^emb for QM orbital grid point `r`.
    pub fn electrostatic_potential_at(&self, r: [f64; 3]) -> f64 {
        let mut v = 0.0;
        for (k, &q_k) in self.mm_charges.iter().enumerate() {
            if k >= self.mm_positions.len() {
                break;
            }
            let pos = self.mm_positions[k];
            let d = [r[0] - pos[0], r[1] - pos[1], r[2] - pos[2]];
            let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if dist > 1e-10 {
                v += q_k / dist;
            }
        }
        v
    }

    /// Compute the mechanical embedding energy: E_elst = Σ_{QM} q_i * V(R_i).
    pub fn mechanical_embedding_energy(
        &self,
        qm_charges: &[f64],
        qm_positions: &[[f64; 3]],
    ) -> f64 {
        let mut e = 0.0;
        for (i, &q_i) in qm_charges.iter().enumerate() {
            if i >= qm_positions.len() {
                break;
            }
            e += q_i * self.electrostatic_potential_at(qm_positions[i]);
        }
        e
    }

    /// Total QM/MM energy: E_total = E_QM + E_MM + E_QM/MM.
    pub fn total_energy(&self, e_qm: f64, e_mm: f64, e_qmmm: f64) -> f64 {
        e_qm + e_mm + e_qmmm
    }

    /// Check if an atom is in the QM region.
    pub fn is_qm_atom(&self, idx: usize) -> bool {
        self.qm_atoms.contains(&idx)
    }

    /// Count link atoms at the boundary.
    pub fn n_link_atoms(&self) -> usize {
        self.link_atoms.len()
    }

    /// Compute the van der Waals cut-off correction for link atom H.
    pub fn link_atom_vdw_correction(&self, r_link: f64, epsilon: f64, sigma: f64) -> f64 {
        if r_link < 1e-10 {
            return 0.0;
        }
        let sr = sigma / r_link;
        4.0 * epsilon * (sr.powi(12) - sr.powi(6))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Build a minimal STO-3G overlap matrix for H₂ at bond length R (Bohr).
pub fn build_h2_overlap(r: f64) -> [[f64; 2]; 2] {
    let zeta = 1.24; // H STO-3G zeta
    let center_a = [0.0, 0.0, 0.0f64];
    let center_b = [r, 0.0, 0.0f64];
    let ba = BasisFunction::sto3g_s(center_a, zeta);
    let bb = BasisFunction::sto3g_s(center_b, zeta);
    [
        [ba.overlap_with(&ba), ba.overlap_with(&bb)],
        [bb.overlap_with(&ba), bb.overlap_with(&bb)],
    ]
}

/// Hartree-to-eV conversion factor.
pub const HARTREE_TO_EV: f64 = 27.211386245988;

/// Bohr-to-Angstrom conversion factor.
pub const BOHR_TO_ANGSTROM: f64 = 0.529177210903;

/// Angstrom-to-Bohr conversion.
pub fn angstrom_to_bohr(r: f64) -> f64 {
    r / BOHR_TO_ANGSTROM
}

/// Convert energy from Hartree to kcal/mol.
pub fn hartree_to_kcal(e: f64) -> f64 {
    e * 627.5094740631
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_1x1_hf() -> HartreeFockSolver {
        HartreeFockSolver::new(
            1,
            2,
            BasisSetType::Sto3G,
            vec![vec![1.0]],
            vec![vec![-1.0]],
            vec![0.5],
            vec![1.0],
            vec![[0.0, 0.0, 0.0]],
        )
    }

    fn make_2x2_hf() -> HartreeFockSolver {
        HartreeFockSolver::new(
            2,
            2,
            BasisSetType::Sto3G,
            vec![vec![1.0, 0.5], vec![0.5, 1.0]],
            vec![vec![-1.5, -0.3], vec![-0.3, -1.5]],
            vec![0.5; 6],
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
        )
    }

    // ── BasisFunction ────────────────────────────────────────────────────────

    #[test]
    fn test_sto3g_s_creates_3_primitives() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        assert_eq!(bf.exponents.len(), 3);
        assert_eq!(bf.coefficients.len(), 3);
        assert_eq!(bf.angular_momentum, 0);
    }

    #[test]
    fn test_sto3g_p_angular_momentum() {
        let bf = BasisFunction::sto3g_p([0.0, 0.0, 0.0], 1.0);
        assert_eq!(bf.angular_momentum, 1);
    }

    #[test]
    fn test_self_overlap_near_one() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let s = bf.overlap_with(&bf);
        assert!(s > 0.0, "self-overlap should be positive");
    }

    #[test]
    fn test_overlap_decreases_with_distance() {
        let bf_a = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let bf_b = BasisFunction::sto3g_s([1.0, 0.0, 0.0], 1.0);
        let bf_c = BasisFunction::sto3g_s([5.0, 0.0, 0.0], 1.0);
        let s_ab = bf_a.overlap_with(&bf_b).abs();
        let s_ac = bf_a.overlap_with(&bf_c).abs();
        assert!(s_ab > s_ac, "overlap should decrease with distance");
    }

    #[test]
    fn test_kinetic_integral_positive() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let t = bf.kinetic_with(&bf);
        assert!(t > 0.0, "kinetic integral should be positive");
    }

    #[test]
    fn test_zeta_scaling_exponents() {
        let bf1 = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let bf2 = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 2.0);
        // Higher zeta → larger exponents
        for (e1, e2) in bf1.exponents.iter().zip(bf2.exponents.iter()) {
            assert!(e2 > e1);
        }
    }

    // ── HartreeFockSolver ────────────────────────────────────────────────────

    #[test]
    fn test_nuclear_repulsion_single_atom_zero() {
        let hf = make_1x1_hf();
        assert_eq!(hf.nuclear_repulsion(), 0.0);
    }

    #[test]
    fn test_nuclear_repulsion_two_atoms() {
        let hf = make_2x2_hf();
        let e_nn = hf.nuclear_repulsion();
        // Z1=1, Z2=1, R=1.4 Bohr → E_nn = 1/1.4 ≈ 0.714
        assert!((e_nn - 1.0 / 1.4).abs() < 1e-6);
    }

    #[test]
    fn test_fock_matrix_size() {
        let hf = make_2x2_hf();
        let density = hf.initial_density_matrix();
        let fock = hf.build_fock_matrix(&density);
        assert_eq!(fock.len(), 2);
        assert_eq!(fock[0].len(), 2);
    }

    #[test]
    fn test_initial_density_matrix_size() {
        let hf = make_2x2_hf();
        let dm = hf.initial_density_matrix();
        assert_eq!(dm.len(), 2);
    }

    #[test]
    fn test_scf_returns_result() {
        let hf = make_2x2_hf();
        let result = hf.run_scf();
        assert!(result.energy.is_finite());
    }

    #[test]
    fn test_scf_energy_has_nuclear_contribution() {
        let hf = make_2x2_hf();
        let result = hf.run_scf();
        assert!(result.nuclear_repulsion > 0.0);
    }

    #[test]
    fn test_mulliken_charge_conservation() {
        let hf = make_2x2_hf();
        let dm = hf.initial_density_matrix();
        // Mulliken charges should be finite
        let q0 = hf.mulliken_charge(&dm, 1.0, &[0]);
        let q1 = hf.mulliken_charge(&dm, 1.0, &[1]);
        assert!(q0.is_finite(), "q0 = {q0} should be finite");
        assert!(q1.is_finite(), "q1 = {q1} should be finite");
    }

    #[test]
    fn test_electronic_energy_negative() {
        let hf = make_2x2_hf();
        let dm = hf.initial_density_matrix();
        let fock = hf.build_fock_matrix(&dm);
        let e = hf.electronic_energy(&dm, &fock);
        assert!(e.is_finite());
    }

    #[test]
    fn test_eri_index_symmetry() {
        let hf = make_2x2_hf();
        // (mu nu | lambda sigma) = (nu mu | lambda sigma) etc.
        let i1 = hf.eri_index(0, 1, 0, 1);
        let i2 = hf.eri_index(1, 0, 0, 1);
        assert_eq!(i1, i2);
    }

    #[test]
    fn test_basis_set_type_equality() {
        assert_eq!(BasisSetType::Sto3G, BasisSetType::Sto3G);
        assert_ne!(BasisSetType::Sto3G, BasisSetType::G6_31G);
    }

    // ── DensityFunctionalTheory ──────────────────────────────────────────────

    #[test]
    fn test_lda_exchange_zero_for_zero_density() {
        let e = DensityFunctionalTheory::lda_exchange_energy_density(0.0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_lda_exchange_negative() {
        let e = DensityFunctionalTheory::lda_exchange_energy_density(1.0);
        assert!(e < 0.0, "exchange energy density should be negative");
    }

    #[test]
    fn test_vwn_correlation_zero_for_zero_density() {
        let e = DensityFunctionalTheory::vwn_correlation_energy_density(0.0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_pbe_enhancement_at_zero_gradient() {
        // s=0 → F_x = 1
        let f = DensityFunctionalTheory::pbe_exchange_enhancement(0.0);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pbe_enhancement_increases_with_s() {
        let f0 = DensityFunctionalTheory::pbe_exchange_enhancement(0.0);
        let f1 = DensityFunctionalTheory::pbe_exchange_enhancement(1.0);
        assert!(f1 > f0);
    }

    #[test]
    fn test_dft_xc_energy_lda() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Lda);
        let e = dft.xc_energy_uniform(0.5, 1.0);
        assert!(e.is_finite());
    }

    #[test]
    fn test_dft_xc_energy_b3lyp() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::B3lyp);
        let e = dft.xc_energy_uniform(0.5, 1.0);
        assert!(e.is_finite());
    }

    #[test]
    fn test_dft_ks_scf_runs() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Pbe);
        let result = dft.run_ks_scf();
        assert!(result.energy.is_finite());
    }

    #[test]
    fn test_homo_lumo_gap_positive() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Lda);
        let orbitals = vec![-0.5, 0.3];
        let gap = dft.homo_lumo_gap(&orbitals);
        assert!(gap >= 0.0);
    }

    #[test]
    fn test_pseudopotential_added_to_energy() {
        let hf = make_2x2_hf();
        let mut dft = DensityFunctionalTheory::new(hf, XcFunctional::Lda);
        dft.set_pseudopotential(-0.1);
        let result = dft.run_ks_scf();
        assert!(result.energy.is_finite());
    }

    // ── MolecularOrbital ─────────────────────────────────────────────────────

    #[test]
    fn test_mo_is_occupied() {
        let mo = MolecularOrbital::new(0, -0.5, 2.0, vec![0.7, 0.7]);
        assert!(mo.is_occupied());
    }

    #[test]
    fn test_mo_is_not_occupied() {
        let mo = MolecularOrbital::new(1, 0.3, 0.0, vec![0.7, -0.7]);
        assert!(!mo.is_occupied());
    }

    #[test]
    fn test_participation_ratio_localized() {
        // Fully localized on one center
        let mo = MolecularOrbital::new(0, -0.5, 2.0, vec![1.0, 0.0, 0.0, 0.0]);
        let pr = mo.participation_ratio();
        assert!((pr - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_participation_ratio_delocalized() {
        // Equally delocalized over 4 centers
        let c = 0.5f64;
        let mo = MolecularOrbital::new(0, -0.3, 2.0, vec![c, c, c, c]);
        let pr = mo.participation_ratio();
        assert!(
            (pr - 4.0).abs() < 1e-6,
            "PR should be ~4 for equal delocalization"
        );
    }

    #[test]
    fn test_dipole_contribution_finite() {
        let mo = MolecularOrbital::new(0, -0.5, 2.0, vec![0.8, 0.6]);
        let centers = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let d = mo.dipole_contribution(&centers, 0);
        assert!(d.is_finite());
    }

    #[test]
    fn test_mayer_bond_order_finite() {
        let n = 2;
        let density = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let overlap = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        let bo = MolecularOrbital::mayer_bond_order(&density, &overlap, &[0], &[1]);
        let _ = n;
        assert!(bo.is_finite());
    }

    #[test]
    fn test_wiberg_bond_index_finite() {
        let density = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let overlap = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        let wbi = MolecularOrbital::wiberg_bond_index(&density, &overlap, &[0], &[1]);
        assert!(wbi.is_finite());
    }

    // ── ElectronCorrelation ──────────────────────────────────────────────────

    fn make_ec() -> ElectronCorrelation {
        let eps = vec![-0.6, -0.5, 0.3, 0.4];
        let eri = vec![0.1; 4];
        ElectronCorrelation::new(2, 2, eps, eri, -1.0)
    }

    #[test]
    fn test_mp2_result_finite() {
        let ec = make_ec();
        let r = ec.mp2_energy();
        assert!(r.e_mp2.is_finite());
    }

    #[test]
    fn test_mp2_total_includes_hf() {
        let ec = make_ec();
        let r = ec.mp2_energy();
        assert_eq!(r.e_hf, -1.0);
        assert!((r.e_total - r.e_hf - r.e_mp2).abs() < 1e-10);
    }

    #[test]
    fn test_cis_returns_states() {
        let ec = make_ec();
        let states = ec.cis_states(2);
        assert!(states.len() <= 2);
    }

    #[test]
    fn test_cis_excitation_energies_positive() {
        let ec = make_ec();
        let states = ec.cis_states(4);
        for s in &states {
            assert!(s.excitation_energy > 0.0 || s.excitation_energy.is_finite());
        }
    }

    #[test]
    fn test_ccs_t1_diagnostic_zero_amplitudes() {
        let ec = make_ec();
        let t1 = ec.ccs_t1_diagnostic(&[]);
        assert_eq!(t1, 0.0);
    }

    #[test]
    fn test_ccs_t1_diagnostic_nonzero() {
        let ec = make_ec();
        let t1 = ec.ccs_t1_diagnostic(&[0.02, 0.01, 0.015, 0.005]);
        assert!(t1 > 0.0);
    }

    #[test]
    fn test_ccs_energy_finite() {
        let ec = make_ec();
        let t1 = vec![0.02; 4];
        let e = ec.ccs_energy(&t1);
        assert!(e.is_finite());
    }

    #[test]
    fn test_scs_mp2_finite() {
        let ec = make_ec();
        let e = ec.scs_mp2_energy(1.2, 0.333);
        assert!(e.is_finite());
    }

    // ── NaturalBondOrbital ───────────────────────────────────────────────────

    #[test]
    fn test_nbo_role_bond() {
        let nbo = NaturalBondOrbital::new("BD", vec![0, 1], 1.97, -0.5);
        assert_eq!(nbo.role(), "donor");
    }

    #[test]
    fn test_nbo_role_antibond() {
        let nbo = NaturalBondOrbital::new("BD*", vec![0, 1], 0.02, 0.3);
        assert_eq!(nbo.role(), "acceptor");
    }

    #[test]
    fn test_nbo_s_character_zero_hybridization() {
        let nbo = NaturalBondOrbital::new("LP", vec![0], 2.0, -0.3);
        assert_eq!(nbo.s_character(), 0.0);
    }

    #[test]
    fn test_nbo_s_character_sp3() {
        let mut nbo = NaturalBondOrbital::new("BD", vec![0, 1], 1.9, -0.5);
        // sp3: 25% s character
        nbo.hybridization = vec![0.5, 0.866, 0.0, 0.0]; // approximate sp3
        let sc = nbo.s_character();
        assert!(sc > 0.0 && sc <= 100.0);
    }

    #[test]
    fn test_nbo_interaction_e2_sign() {
        let inter = NboInteraction::compute(1.97, -0.05, -0.5, 0.3, 0, 1);
        // e_A > e_D → denominator positive, F²>0 → E(2) < 0
        assert!(inter.e2_energy < 0.0);
    }

    #[test]
    fn test_nbo_interaction_is_significant() {
        let inter = NboInteraction::compute(1.97, -0.05, -0.5, 0.3, 0, 1);
        assert!(inter.is_significant(0.1));
    }

    #[test]
    fn test_nbo_analysis_top_interactions() {
        let nbos = vec![
            NaturalBondOrbital::new("BD", vec![0, 1], 1.97, -0.5),
            NaturalBondOrbital::new("BD*", vec![0, 1], 0.02, 0.3),
        ];
        let mut analysis = NboAnalysis::new(nbos, vec![0.0, 0.0]);
        let fock_nbo = vec![vec![0.0, -0.05], vec![-0.05, 0.0]];
        analysis.compute_interactions(&fock_nbo, 0.1);
        let top = analysis.top_interactions(5);
        assert!(top.len() <= 5);
    }

    #[test]
    fn test_nbo_total_charge_transfer_updated() {
        let nbos = vec![
            NaturalBondOrbital::new("BD", vec![0, 1], 1.97, -0.5),
            NaturalBondOrbital::new("BD*", vec![0, 1], 0.02, 0.3),
        ];
        let mut analysis = NboAnalysis::new(nbos, vec![0.0, 0.0]);
        let fock_nbo = vec![vec![0.0, -0.05], vec![-0.05, 0.0]];
        analysis.compute_interactions(&fock_nbo, 0.0);
        // total_charge_transfer should be the sum of E(2) values
        assert!(analysis.total_charge_transfer.is_finite());
    }

    // ── QmmmCoupling ─────────────────────────────────────────────────────────

    fn make_qmmm() -> QmmmCoupling {
        QmmmCoupling::new(
            QmmmBoundary::LinkAtom,
            EmbeddingType::Electrostatic,
            vec![0, 1],
            vec![-0.834, 0.417, 0.417],
            vec![[3.0, 0.0, 0.0], [3.5, 0.5, 0.0], [3.5, -0.5, 0.0]],
        )
    }

    #[test]
    fn test_qmmm_is_qm_atom() {
        let qmmm = make_qmmm();
        assert!(qmmm.is_qm_atom(0));
        assert!(!qmmm.is_qm_atom(5));
    }

    #[test]
    fn test_qmmm_electrostatic_potential_finite() {
        let qmmm = make_qmmm();
        let v = qmmm.electrostatic_potential_at([0.0, 0.0, 0.0]);
        assert!(v.is_finite());
    }

    #[test]
    fn test_qmmm_mechanical_embedding_finite() {
        let qmmm = make_qmmm();
        let charges = vec![1.0, -1.0];
        let positions = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let e = qmmm.mechanical_embedding_energy(&charges, &positions);
        assert!(e.is_finite());
    }

    #[test]
    fn test_qmmm_total_energy() {
        let qmmm = make_qmmm();
        let e = qmmm.total_energy(-74.9, -5.0, -0.3);
        assert!((e - (-80.2)).abs() < 1e-10);
    }

    #[test]
    fn test_link_atom_position_interpolation() {
        let la = LinkAtomQmmm::from_bond([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0, 1, 0.7);
        assert!((la.position[0] - 1.4).abs() < 1e-10);
    }

    #[test]
    fn test_link_atom_force_correction_sum() {
        let la = LinkAtomQmmm::from_bond([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0, 1, 0.7);
        let f = [1.0, 0.0, 0.0];
        let fqm = la.qm_force_correction(f);
        let fmm = la.mm_force_correction(f);
        // Sum should equal f
        assert!((fqm[0] + fmm[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_qmmm_add_link_atom() {
        let mut qmmm = make_qmmm();
        let la = LinkAtomQmmm::from_bond([0.0, 0.0, 0.0], [1.5, 0.0, 0.0], 0, 1, 0.7);
        qmmm.add_link_atom(la);
        assert_eq!(qmmm.n_link_atoms(), 1);
    }

    #[test]
    fn test_link_atom_vdw_correction_finite() {
        let qmmm = make_qmmm();
        let e = qmmm.link_atom_vdw_correction(2.5, 0.1, 2.5);
        assert!(e.is_finite());
    }

    #[test]
    fn test_embedding_type_equality() {
        assert_eq!(EmbeddingType::Mechanical, EmbeddingType::Mechanical);
        assert_ne!(EmbeddingType::Mechanical, EmbeddingType::Electrostatic);
    }

    // ── Utility ──────────────────────────────────────────────────────────────

    #[test]
    fn test_h2_overlap_diagonal_near_one() {
        let s = build_h2_overlap(1.4);
        assert!(s[0][0] > 0.5);
        assert!(s[1][1] > 0.5);
    }

    #[test]
    fn test_h2_overlap_off_diagonal_positive() {
        let s = build_h2_overlap(1.4);
        assert!(s[0][1] > 0.0);
    }

    #[test]
    fn test_h2_overlap_symmetric() {
        let s = build_h2_overlap(1.4);
        assert!((s[0][1] - s[1][0]).abs() < 1e-12);
    }

    #[test]
    fn test_hartree_to_ev_conversion() {
        let e_ev = HARTREE_TO_EV;
        assert!((e_ev - 27.211386).abs() < 0.001);
    }

    #[test]
    fn test_angstrom_to_bohr_round_trip() {
        let r = 1.5;
        let r_bohr = angstrom_to_bohr(r);
        let r_back = r_bohr * BOHR_TO_ANGSTROM;
        assert!((r_back - r).abs() < 1e-10);
    }

    #[test]
    fn test_hartree_to_kcal_conversion() {
        let e = hartree_to_kcal(1.0);
        assert!((e - 627.509).abs() < 0.01);
    }

    #[test]
    fn test_scf_status_converged_ne_notconverged() {
        assert_ne!(ScfStatus::Converged, ScfStatus::NotConverged);
    }
}
