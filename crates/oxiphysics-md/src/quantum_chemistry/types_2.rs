//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::primitive_eri_contracted;
use super::types::BasisFunction;
use super::types_3::{DftResult, HartreeFockResult, ScfStatus, XcFunctional};

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
        let cx = -0.7385587663820224;
        cx * rho.powf(1.0 / 3.0)
    }
    /// VWN LDA correlation energy density (simplified Vosko-Wilk-Nusair).
    pub fn vwn_correlation_energy_density(rho: f64) -> f64 {
        if rho < 1e-20 {
            return 0.0;
        }
        let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);
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
        let mu = 0.2195149727645171;
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
        let mut fock = self.core_hamiltonian.clone();
        for (mu, fock_row) in fock.iter_mut().enumerate() {
            for (nu, fock_val) in fock_row.iter_mut().enumerate() {
                let mut g_mn = 0.0;
                for (lambda, density_row) in density.iter().enumerate() {
                    for (sigma, &pls) in density_row.iter().enumerate() {
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
                *fock_val += g_mn;
            }
        }
        fock
    }
    /// ERI index for (mu nu | lambda sigma) with 8-fold symmetry.
    pub(crate) fn eri_index(&self, mu: usize, nu: usize, lambda: usize, sigma: usize) -> usize {
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
        let _ = n;
        pair
    }
    /// Compute the initial (guess) density matrix using core Hamiltonian diagonalization.
    pub fn initial_density_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_basis;
        let n_occ = self.n_electrons / 2;
        let coeffs = self.diagonalize_symmetric(&self.core_hamiltonian, &self.overlap);
        Self::build_density_from_coeffs(&coeffs, n_occ, n)
    }
    /// Build density matrix P_mn = 2 * Σ_{occ} C_mi * C_ni
    fn build_density_from_coeffs(coeffs: &[Vec<f64>], n_occ: usize, n: usize) -> Vec<Vec<f64>> {
        let mut density = vec![vec![0.0; n]; n];
        for mu in 0..n {
            for nu in 0..n {
                let mut p = 0.0;
                for coeff in coeffs.iter().take(n_occ.min(coeffs.len())) {
                    if coeff.len() > mu && coeff.len() > nu {
                        p += 2.0 * coeff[mu] * coeff[nu];
                    }
                }
                density[mu][nu] = p;
            }
        }
        density
    }
    /// Generalized symmetric eigenproblem solver: returns eigenvectors as rows sorted by
    /// ascending eigenvalue.  Solves the Roothaan-Hall equation FC = SCε correctly by
    /// using `generalized_symmetric_eigen_n` (Löwdin orthogonalisation via Cholesky) when
    /// the overlap matrix S is non-trivial, falling back to the standard solver otherwise.
    fn diagonalize_symmetric(&self, matrix: &[Vec<f64>], overlap: &[Vec<f64>]) -> Vec<Vec<f64>> {
        use oxiphysics_core::numerical_methods::{
            generalized_symmetric_eigen_n, symmetric_eigen_n,
        };
        let result = generalized_symmetric_eigen_n(matrix, overlap);
        let (_evals, evecs) = if let Some(r) = result {
            r
        } else {
            symmetric_eigen_n(matrix)
        };
        evecs
    }
    /// Diagonalize the generalised eigenvalue problem and return both eigenvalues and eigenvectors.
    ///
    /// Returns `(eigenvalues, eigenvectors)` where `eigenvectors[i]` is the i-th eigenvector
    /// (row) corresponding to `eigenvalues[i]`, sorted in ascending order.
    fn diagonalize_and_get_energies(
        &self,
        matrix: &[Vec<f64>],
        overlap: &[Vec<f64>],
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        use oxiphysics_core::numerical_methods::{
            generalized_symmetric_eigen_n, symmetric_eigen_n,
        };
        if let Some((evals, evecs)) = generalized_symmetric_eigen_n(matrix, overlap) {
            (evals, evecs)
        } else {
            symmetric_eigen_n(matrix)
        }
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
            let (evals, evecs) = self.diagonalize_and_get_energies(&fock, &self.overlap);
            if iter > 0 && (energy - energy_old).abs() < self.energy_threshold {
                status = ScfStatus::Converged;
                let len = n.min(evals.len());
                orbital_energies[..len].copy_from_slice(&evals[..len]);
                mo_coefficients = evecs;
                break;
            }
            energy_old = energy;
            density = Self::build_density_from_coeffs(&evecs, n_occ, n);
            if iter == self.max_iter - 1 {
                let len = n.min(evals.len());
                orbital_energies[..len].copy_from_slice(&evals[..len]);
                mo_coefficients = evecs;
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
            for (nu, &d_val) in density[mu].iter().enumerate().take(n) {
                ps_mu += d_val * self.overlap[mu][nu];
            }
            q_a -= ps_mu;
        }
        q_a
    }
    /// Compute two-electron repulsion integrals (μν|λσ) for a given basis set.
    ///
    /// Uses the exact Boys-function formula for s-type primitive pairs (ss|ss).
    /// For p-type contributions the zeroth-order (monopole) term is computed via
    /// the product-Gaussian center — a valid approximation for STO-3G.
    /// Returns a flat `Vec<f64>` indexed via `Self::eri_index`.
    pub fn compute_eris(&self, basis: &[BasisFunction]) -> Vec<f64> {
        let n = basis.len();
        let n_pair = n * (n + 1) / 2;
        let n_eri = n_pair * (n_pair + 1) / 2;
        let mut eri = vec![0.0_f64; n_eri];
        for mu in 0..n {
            for nu in 0..=mu {
                for lambda in 0..n {
                    for sigma in 0..=lambda {
                        let idx = self.eri_index(mu, nu, lambda, sigma);
                        if idx >= eri.len() {
                            continue;
                        }
                        let val = primitive_eri_contracted(
                            &basis[mu],
                            &basis[nu],
                            &basis[lambda],
                            &basis[sigma],
                        );
                        eri[idx] = val;
                    }
                }
            }
        }
        eri
    }
    /// STO-3G Slater exponent parameters per element (Hehre, Ditchfield, Pople 1969).
    ///
    /// Returns `(zeta_s, Some(zeta_p))` for second-row atoms and `(zeta_s, None)` for hydrogen.
    fn sto3g_zeta(z: u32) -> Option<(f64, Option<f64>)> {
        match z {
            1 => Some((1.24, None)),
            6 => Some((5.67, Some(1.72))),
            7 => Some((6.67, Some(1.95))),
            8 => Some((7.66, Some(2.25))),
            9 => Some((8.65, Some(2.55))),
            16 => Some((6.92, Some(1.89))),
            _ => None,
        }
    }
    /// Number of valence electrons contributed by element `z` in a minimal basis treatment.
    fn valence_electrons(z: u32) -> usize {
        match z {
            1 => 1,
            6 => 4,
            7 => 5,
            8 => 6,
            9 => 7,
            16 => 6,
            _ => 0,
        }
    }
    /// Construct a `HartreeFockSolver` from atomic positions and atomic numbers.
    ///
    /// Builds the STO-3G minimal basis, computes the overlap matrix S, the core Hamiltonian
    /// H_core = T + V_ne (using the Boys function for nuclear attraction integrals), and the
    /// two-electron repulsion integrals (ERIs) from first principles.
    ///
    /// # Arguments
    /// * `atomic_numbers` — element identifiers (1=H, 6=C, 7=N, 8=O, 9=F, 16=S).
    /// * `positions_angstrom` — Cartesian positions in **Ångström**; converted to Bohr internally.
    /// * `total_charge` — net molecular charge (0 for neutral).
    ///
    /// Returns `None` if an unsupported element is encountered.
    pub fn from_atoms(
        atomic_numbers: &[u32],
        positions_angstrom: &[[f64; 3]],
        total_charge: i32,
    ) -> Option<Self> {
        const BOHR: f64 = 1.889_725_988_6;
        let nat = atomic_numbers.len();
        if nat == 0 || nat != positions_angstrom.len() {
            return None;
        }
        let pos_bohr: Vec<[f64; 3]> = positions_angstrom
            .iter()
            .map(|p| [p[0] * BOHR, p[1] * BOHR, p[2] * BOHR])
            .collect();
        let mut basis: Vec<BasisFunction> = Vec::new();
        let mut n_electrons: usize = 0;
        for (a, &z) in atomic_numbers.iter().enumerate() {
            let (zs, zp_opt) = Self::sto3g_zeta(z)?;
            let center = pos_bohr[a];
            basis.push(BasisFunction::sto3g_s(center, zs));
            if let Some(zp) = zp_opt {
                for _ in 0..3 {
                    basis.push(BasisFunction::sto3g_p(center, zp));
                }
            }
            n_electrons += Self::valence_electrons(z);
        }
        let n_electrons = (n_electrons as i64 - total_charge as i64).max(0) as usize;
        let n_basis = basis.len();
        let overlap: Vec<Vec<f64>> = (0..n_basis)
            .map(|i| {
                (0..n_basis)
                    .map(|j| basis[i].overlap_with(&basis[j]))
                    .collect()
            })
            .collect();
        let mut h_core = vec![vec![0.0_f64; n_basis]; n_basis];
        for i in 0..n_basis {
            for j in 0..n_basis {
                let t_ij = basis[i].kinetic_with(&basis[j]);
                let mut v_ij = 0.0;
                for (a, &z) in atomic_numbers.iter().enumerate() {
                    v_ij += basis[i].nuclear_attraction_with(&basis[j], z as f64, pos_bohr[a]);
                }
                h_core[i][j] = t_ij + v_ij;
            }
        }
        let mut solver = Self::new(
            n_basis,
            n_electrons,
            BasisSetType::Sto3G,
            overlap,
            h_core,
            vec![0.0],
            atomic_numbers.iter().map(|&z| z as f64).collect(),
            pos_bohr,
        );
        let eri = solver.compute_eris(&basis);
        solver.eri = eri;
        Some(solver)
    }
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
            -n_donor * fock_da * fock_da / gap * 627.5094740631
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
