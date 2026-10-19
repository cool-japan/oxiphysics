//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::error::QuantumChemistryError;
use super::integrals::electron_repulsion_integral;
use super::kohn_sham::{run_kohn_sham, symmetrize_from_lower};
use super::molecular_grid::MolecularGridSettings;
use super::sto3g_basis::sto3g_atom_basis;
use super::types::BasisFunction;
use super::types_3::{DftResult, HartreeFockResult, ScfStatus, XcFunctional};
use super::xc_functional::{
    B3LYP_EXACT_EXCHANGE, VWN5, XcDual, evaluate_xc, uniform_gas_tau, vwn_epsilon,
};

/// Restricted Kohn–Sham density functional theory solver.
///
/// The Kohn–Sham equations are solved self-consistently on a Becke molecular
/// grid ([`MolecularGrid`]): the exchange-correlation potential matrix V_xc is
/// rebuilt from the current density in every iteration and enters the
/// Kohn–Sham matrix `F = H_core + J + V_xc − (a₀/2) K` (see
/// [`run_ks_scf`](Self::run_ks_scf)).  The integrals (S, H_core, ERIs) and the
/// basis functions are taken from `hf_solver`, which must have been built with
/// [`HartreeFockSolver::from_atoms`] or [`HartreeFockSolver::from_basis`].
#[derive(Debug, Clone)]
pub struct DensityFunctionalTheory {
    /// Integral and basis-set provider (also supplies `max_iter` and the
    /// convergence thresholds of the SCF).
    pub hf_solver: HartreeFockSolver,
    /// Exchange-correlation functional.
    pub xc: XcFunctional,
    /// Fraction a₀ of exact (Hartree–Fock) exchange: 0.20 for B3LYP, 0 for the
    /// pure functionals.  A non-zero value scales the semilocal exchange by
    /// (1 − a₀) and adds −(a₀/2) K to the Kohn–Sham matrix.
    pub alpha_x: f64,
    /// Pseudopotential correction energy \[Hartree\].
    pub pseudopotential_energy: f64,
    /// Resolution of the molecular integration grid.
    pub grid: MolecularGridSettings,
}
impl DensityFunctionalTheory {
    /// Create a KS-DFT solver with the default grid.
    pub fn new(hf_solver: HartreeFockSolver, xc: XcFunctional) -> Self {
        Self {
            hf_solver,
            xc,
            alpha_x: match xc {
                XcFunctional::B3lyp => B3LYP_EXACT_EXCHANGE,
                _ => 0.0,
            },
            pseudopotential_energy: 0.0,
            grid: MolecularGridSettings::default(),
        }
    }
    /// Use a different integration grid.
    pub fn with_grid(mut self, grid: MolecularGridSettings) -> Self {
        self.grid = grid;
        self
    }
    /// LDA Slater exchange energy per electron: ε_x = -3/4 * (3/π)^(1/3) * ρ^(1/3).
    pub fn lda_exchange_energy_density(rho: f64) -> f64 {
        if rho < 1e-20 {
            return 0.0;
        }
        let cx = -0.7385587663820224;
        cx * rho.powf(1.0 / 3.0)
    }
    /// VWN5 (Vosko–Wilk–Nusair functional V, paramagnetic) correlation energy
    /// per electron ε_c(ρ) \[Hartree\].
    pub fn vwn_correlation_energy_density(rho: f64) -> f64 {
        if rho < 1e-20 {
            return 0.0;
        }
        let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);
        vwn_epsilon(XcDual::constant(rs), &VWN5).v
    }
    /// PBE exchange enhancement factor F_x(s) where s is reduced gradient.
    pub fn pbe_exchange_enhancement(s: f64) -> f64 {
        let kappa = 0.804;
        let mu = 0.2195149727645171;
        1.0 + kappa - kappa / (1.0 + mu * s * s / kappa)
    }
    /// Exchange-correlation energy of a *uniform* electron gas of density `rho`
    /// (Bohr⁻³) filling `volume` (Bohr³).
    ///
    /// Evaluates the selected functional at ∇ρ = 0 and τ = τ_unif(ρ); for a
    /// uniform gas the exact exchange of a hybrid equals Slater exchange, so the
    /// a₀ fraction is added back as Slater exchange.  (GGAs and TPSS reduce to
    /// their LSDA limit; B3LYP keeps its LYP and VWN(RPA) correlation.)
    pub fn xc_energy_uniform(&self, rho: f64, volume: f64) -> f64 {
        let point = evaluate_xc(self.xc, self.alpha_x, rho, 0.0, uniform_gas_tau(rho));
        let exact_exchange = self.alpha_x * rho * Self::lda_exchange_energy_density(rho);
        (point.energy_density + exact_exchange) * volume
    }
    /// Apply a pseudopotential correction to the total energy.
    pub fn set_pseudopotential(&mut self, e_pp: f64) {
        self.pseudopotential_energy = e_pp;
    }
    /// Run the restricted Kohn–Sham SCF.
    ///
    /// Starting from the core-Hamiltonian guess, each iteration builds
    /// J, (a₀ K) and V_xc on the molecular grid from the current density,
    /// forms F = H_core + J − (a₀/2) K + V_xc, extrapolates it with DIIS and
    /// diagonalises FC = SCε.  The total energy is
    /// E = Tr(PH) + ½Tr(PJ) − (a₀/4)Tr(PK) + E_xc + E_nn (+ pseudopotential).
    /// Convergence requires |ΔE| < `hf_solver.energy_threshold` and
    /// max|FPS − SPF| < `hf_solver.density_threshold`; otherwise the result is
    /// returned with [`ScfStatus::NotConverged`].
    ///
    /// # Errors
    /// [`QuantumChemistryError::MissingData`] if the solver carries no basis
    /// functions, [`QuantumChemistryError::OpenShellReference`] for an odd
    /// electron count, [`QuantumChemistryError::DimensionMismatch`] for
    /// inconsistent integrals, [`QuantumChemistryError::SingularOverlap`] for a
    /// linearly dependent basis and [`QuantumChemistryError::CoincidentNuclei`]
    /// if no grid can be built.
    pub fn run_ks_scf(&self) -> Result<DftResult, QuantumChemistryError> {
        run_kohn_sham(self)
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
    /// A user-supplied set of contracted Gaussians
    /// (see [`HartreeFockSolver::from_basis`]).
    Custom,
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
    /// Contracted basis functions (filled by [`from_atoms`](Self::from_atoms);
    /// empty when the integrals were supplied directly through [`new`](Self::new)).
    pub basis_functions: Vec<BasisFunction>,
    /// Index of the atom on which each basis function is centred (parallel to `basis_functions`).
    pub basis_atoms: Vec<usize>,
    /// External point charges `(position [Bohr], charge [e])` embedded in the core
    /// Hamiltonian by [`add_external_point_charges`](Self::add_external_point_charges).
    pub external_charges: Vec<([f64; 3], f64)>,
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
            basis_functions: Vec::new(),
            basis_atoms: Vec::new(),
            external_charges: Vec::new(),
        }
    }
    /// Embed external point charges (electrostatic embedding).
    ///
    /// Each charge Q_J at R_J (Bohr) acts on the electrons like a nucleus of charge
    /// Q_J: the core Hamiltonian gains ⟨μ| −Q_J/|r − R_J| |ν⟩ (computed with
    /// [`BasisFunction::nuclear_attraction_with`]) and
    /// [`nuclear_repulsion`](Self::nuclear_repulsion) gains Σ_A Z_A Q_J / R_AJ.
    ///
    /// # Errors
    /// [`QuantumChemistryError::MissingData`] if the solver has no basis functions
    /// (it was not built with [`from_atoms`](Self::from_atoms)),
    /// [`QuantumChemistryError::InvalidExternalCharge`] for a non-finite charge
    /// or one that coincides with a nucleus.
    pub fn add_external_point_charges(
        &mut self,
        charges: &[([f64; 3], f64)],
    ) -> Result<(), QuantumChemistryError> {
        if charges.is_empty() {
            return Ok(());
        }
        if self.basis_functions.len() != self.n_basis {
            return Err(QuantumChemistryError::MissingData("basis functions"));
        }
        for (index, (pos, q)) in charges.iter().enumerate() {
            let coincides = self.atomic_positions.iter().any(|r| {
                (r[0] - pos[0]).powi(2) + (r[1] - pos[1]).powi(2) + (r[2] - pos[2]).powi(2) < 1e-12
            });
            if !q.is_finite() || pos.iter().any(|x| !x.is_finite()) || coincides {
                return Err(QuantumChemistryError::InvalidExternalCharge { index });
            }
        }
        let v_ext = self.external_potential_matrix(charges);
        for (h_row, v_row) in self.core_hamiltonian.iter_mut().zip(v_ext.iter()) {
            for (h, v) in h_row.iter_mut().zip(v_row.iter()) {
                *h += v;
            }
        }
        self.external_charges.extend_from_slice(charges);
        Ok(())
    }
    /// One-electron matrix V_μν = Σ_J ⟨μ| −Q_J/|r − R_J| |ν⟩ of the given point
    /// charges in the solver's basis (zero-sized if the basis is unknown).
    pub fn external_potential_matrix(&self, charges: &[([f64; 3], f64)]) -> Vec<Vec<f64>> {
        let n = self.basis_functions.len();
        let mut v = vec![vec![0.0_f64; n]; n];
        for (i, bf_i) in self.basis_functions.iter().enumerate() {
            for (j, bf_j) in self.basis_functions.iter().enumerate().skip(i) {
                let v_ij: f64 = charges
                    .iter()
                    .map(|(pos, q)| bf_i.nuclear_attraction_with(bf_j, *q, *pos))
                    .sum();
                v[i][j] = v_ij;
                v[j][i] = v_ij;
            }
        }
        v
    }
    /// Compute the nuclear repulsion energy \[Hartree\]: nucleus–nucleus repulsion
    /// plus the interaction of the nuclei with any embedded external point charges.
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
        for (z_a, r_a) in self.atomic_charges.iter().zip(self.atomic_positions.iter()) {
            for (pos, q) in &self.external_charges {
                let r = ((r_a[0] - pos[0]).powi(2)
                    + (r_a[1] - pos[1]).powi(2)
                    + (r_a[2] - pos[2]).powi(2))
                .sqrt();
                if r > 1e-10 {
                    e_nn += z_a * q / r;
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
    pub(super) fn build_density_from_coeffs(
        coeffs: &[Vec<f64>],
        n_occ: usize,
        n: usize,
    ) -> Vec<Vec<f64>> {
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
    /// Every unique integral (8-fold permutational symmetry) is evaluated exactly
    /// with the McMurchie–Davidson scheme
    /// ([`electron_repulsion_integral`](super::integrals::electron_repulsion_integral)),
    /// for s, p and higher Cartesian Gaussians alike.
    /// Returns a flat `Vec<f64>` indexed via `Self::eri_index`.
    pub fn compute_eris(&self, basis: &[BasisFunction]) -> Vec<f64> {
        let n = basis.len();
        let n_pair = n * (n + 1) / 2;
        let n_eri = n_pair * (n_pair + 1) / 2;
        let mut eri = vec![0.0_f64; n_eri];
        for mu in 0..n {
            for nu in 0..=mu {
                let mn = mu * (mu + 1) / 2 + nu;
                for lambda in 0..n {
                    for sigma in 0..=lambda {
                        let ls = lambda * (lambda + 1) / 2 + sigma;
                        if ls > mn {
                            continue;
                        }
                        let idx = self.eri_index(mu, nu, lambda, sigma);
                        if let Some(slot) = eri.get_mut(idx) {
                            *slot = electron_repulsion_integral(
                                &basis[mu],
                                &basis[nu],
                                &basis[lambda],
                                &basis[sigma],
                            );
                        }
                    }
                }
            }
        }
        eri
    }
    /// Construct a `HartreeFockSolver` from atomic positions and atomic numbers.
    ///
    /// Builds the all-electron STO-3G basis ([`sto3g_atom_basis`]: 1s for H/He,
    /// 1s 2s 2p for Li–Ne, 1s 2s 2p 3s 3p for Na–Ar, with genuine Cartesian
    /// p_x/p_y/p_z components), and evaluates S, H_core = T + V_ne and the ERIs
    /// exactly (see [`from_basis`](Self::from_basis)).  The electron count is
    /// Σ Z − `total_charge`, consistent with the full nuclear charges used in
    /// V_ne and the nuclear repulsion.
    ///
    /// # Arguments
    /// * `atomic_numbers` — element identifiers (1 = H … 18 = Ar).
    /// * `positions_angstrom` — Cartesian positions in **Ångström**; converted to Bohr internally.
    /// * `total_charge` — net molecular charge (0 for neutral).
    ///
    /// Returns `None` if an unsupported element is encountered, the input is
    /// empty or inconsistent, or the charge exceeds the total nuclear charge.
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
        let mut basis_atoms: Vec<usize> = Vec::new();
        for (a, (&z, &center)) in atomic_numbers.iter().zip(pos_bohr.iter()).enumerate() {
            let shells = sto3g_atom_basis(z, center)?;
            basis_atoms.extend(std::iter::repeat_n(a, shells.len()));
            basis.extend(shells);
        }
        let nuclear_total: i64 = atomic_numbers.iter().map(|&z| i64::from(z)).sum();
        let n_electrons = usize::try_from(nuclear_total - i64::from(total_charge)).ok()?;
        let nuclear_charges: Vec<f64> = atomic_numbers.iter().map(|&z| f64::from(z)).collect();
        Self::from_basis(
            basis,
            basis_atoms,
            nuclear_charges,
            pos_bohr,
            n_electrons,
            BasisSetType::Sto3G,
        )
        .ok()
    }
    /// Construct a solver for an arbitrary contracted Cartesian Gaussian basis.
    ///
    /// Computes the overlap matrix S, the core Hamiltonian H_core = T + V_ne and
    /// all two-electron integrals exactly (McMurchie–Davidson).
    ///
    /// * `basis_atoms[μ]` — index of the atom carrying basis function μ.
    /// * `nuclear_charges`, `positions_bohr` — nuclei \[e, Bohr\].
    /// * `n_electrons` — total number of electrons.
    ///
    /// # Errors
    /// [`QuantumChemistryError::DimensionMismatch`] if `basis_atoms` does not
    /// match `basis`, the nuclear arrays differ in length, or an entry of
    /// `basis_atoms` is not a valid atom index; [`QuantumChemistryError::MissingData`]
    /// for an empty basis.
    pub fn from_basis(
        basis: Vec<BasisFunction>,
        basis_atoms: Vec<usize>,
        nuclear_charges: Vec<f64>,
        positions_bohr: Vec<[f64; 3]>,
        n_electrons: usize,
        basis_type: BasisSetType,
    ) -> Result<Self, QuantumChemistryError> {
        if basis.is_empty() {
            return Err(QuantumChemistryError::MissingData("basis functions"));
        }
        if basis_atoms.len() != basis.len() {
            return Err(QuantumChemistryError::DimensionMismatch {
                what: "basis-function atom indices",
                expected: basis.len(),
                found: basis_atoms.len(),
            });
        }
        if positions_bohr.len() != nuclear_charges.len() {
            return Err(QuantumChemistryError::DimensionMismatch {
                what: "nuclear positions",
                expected: nuclear_charges.len(),
                found: positions_bohr.len(),
            });
        }
        if let Some(&bad) = basis_atoms.iter().find(|&&a| a >= nuclear_charges.len()) {
            return Err(QuantumChemistryError::DimensionMismatch {
                what: "atom index of a basis function",
                expected: nuclear_charges.len(),
                found: bad,
            });
        }
        let n_basis = basis.len();
        let mut overlap = vec![vec![0.0_f64; n_basis]; n_basis];
        let mut h_core = vec![vec![0.0_f64; n_basis]; n_basis];
        for i in 0..n_basis {
            for j in i..n_basis {
                let s_ij = basis[i].overlap_with(&basis[j]);
                let t_ij = basis[i].kinetic_with(&basis[j]);
                let v_ij: f64 = nuclear_charges
                    .iter()
                    .zip(positions_bohr.iter())
                    .map(|(&z, &c)| basis[i].nuclear_attraction_with(&basis[j], z, c))
                    .sum();
                overlap[i][j] = s_ij;
                overlap[j][i] = s_ij;
                h_core[i][j] = t_ij + v_ij;
                h_core[j][i] = t_ij + v_ij;
            }
        }
        let mut solver = Self::new(
            n_basis,
            n_electrons,
            basis_type,
            overlap,
            h_core,
            Vec::new(),
            nuclear_charges,
            positions_bohr,
        );
        solver.eri = solver.compute_eris(&basis);
        solver.basis_functions = basis;
        solver.basis_atoms = basis_atoms;
        Ok(solver)
    }
    /// Kinetic-energy matrix T_μν = ⟨μ| −½∇² |ν⟩ of the solver's basis
    /// (zero-sized if the basis functions are unknown).
    pub fn kinetic_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.basis_functions.len();
        let mut t: Vec<Vec<f64>> = self
            .basis_functions
            .iter()
            .enumerate()
            .map(|(i, bf_i)| {
                let mut row = vec![0.0_f64; n];
                for (slot, bf_j) in row.iter_mut().zip(self.basis_functions.iter()).take(i + 1) {
                    *slot = bf_i.kinetic_with(bf_j);
                }
                row
            })
            .collect();
        symmetrize_from_lower(&mut t);
        t
    }
    /// Coulomb matrix J_μν = Σ_λσ P_λσ (μν|λσ).
    pub fn coulomb_matrix(&self, density: &[Vec<f64>]) -> Vec<Vec<f64>> {
        self.two_electron_matrix(density, false)
    }
    /// Exchange matrix K_μν = Σ_λσ P_λσ (μλ|νσ).
    pub fn exchange_matrix(&self, density: &[Vec<f64>]) -> Vec<Vec<f64>> {
        self.two_electron_matrix(density, true)
    }
    fn two_electron_matrix(&self, density: &[Vec<f64>], exchange: bool) -> Vec<Vec<f64>> {
        let n = self.n_basis;
        let mut out = vec![vec![0.0_f64; n]; n];
        for (mu, out_row) in out.iter_mut().enumerate() {
            for (nu, out_val) in out_row.iter_mut().enumerate().take(mu + 1) {
                let mut acc = 0.0;
                for (lambda, row) in density.iter().enumerate().take(n) {
                    for (sigma, &p) in row.iter().enumerate().take(n) {
                        let idx = if exchange {
                            self.eri_index(mu, lambda, nu, sigma)
                        } else {
                            self.eri_index(mu, nu, lambda, sigma)
                        };
                        acc += p * self.eri.get(idx).copied().unwrap_or(0.0);
                    }
                }
                *out_val = acc;
            }
        }
        symmetrize_from_lower(&mut out);
        out
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
