// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Restricted Kohn–Sham self-consistent field on a Becke molecular grid.
//!
//! The Kohn–Sham matrix is
//!
//! ```text
//! F = H_core + J[P] − (a₀/2) K[P] + V_xc[P]
//! V_xc,μν = ∫ [ v_ρ φ_μ φ_ν + 2 v_σ ∇ρ·∇(φ_μ φ_ν) + ½ v_τ ∇φ_μ·∇φ_ν ] dr
//! ```
//!
//! with `ρ = Σ P_μν φ_μ φ_ν`, `∇ρ = 2 Σ P_μν φ_μ ∇φ_ν`,
//! `τ = ½ Σ P_μν ∇φ_μ·∇φ_ν`, and the total energy
//!
//! ```text
//! E = Tr(P H) + ½ Tr(P J) − (a₀/4) Tr(P K) + E_xc + E_nn.
//! ```
//!
//! The iterations are accelerated with Pulay's DIIS on the commutator
//! `FPS − SPF`.

use oxiphysics_core::numerical_methods::generalized_symmetric_eigen_n;

use super::error::QuantumChemistryError;
use super::molecular_grid::MolecularGrid;
use super::types::BasisFunction;
use super::types_2::{DensityFunctionalTheory, HartreeFockSolver};
use super::types_3::{DftResult, ScfStatus, XcFunctional};
use super::xc_functional::{evaluate_xc, xc_needs_gradient, xc_needs_tau};

/// Maximum number of stored DIIS vectors.
const DIIS_DEPTH: usize = 8;

type Matrix = Vec<Vec<f64>>;

/// Exchange-correlation energy, potential matrix and integrated electron number.
pub(crate) struct XcContribution {
    /// Semilocal XC energy \[Hartree\].
    pub(crate) energy: f64,
    /// V_xc matrix in the AO basis.
    pub(crate) matrix: Matrix,
    /// ∫ρ over the grid.
    pub(crate) electrons: f64,
}

/// Evaluate the semilocal XC energy and potential matrix of density `p`.
pub(crate) fn xc_contribution(
    grid: &MolecularGrid,
    basis: &[BasisFunction],
    p: &[Vec<f64>],
    functional: XcFunctional,
    exact_exchange: f64,
) -> XcContribution {
    let n = basis.len();
    let gradient = xc_needs_gradient(functional);
    let meta = xc_needs_tau(functional);
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut energy = 0.0;
    let mut electrons = 0.0;
    let mut phi = vec![0.0_f64; n];
    let mut dphi = vec![[0.0_f64; 3]; n];
    let mut p_phi = vec![0.0_f64; n];
    let mut a = vec![0.0_f64; n];
    for (&point, &w) in grid.points.iter().zip(grid.weights.iter()) {
        for (mu, bf) in basis.iter().enumerate() {
            let (value, grad) = bf.value_and_gradient(point);
            phi[mu] = value;
            dphi[mu] = grad;
        }
        // (Pφ)_μ = Σ_ν P_μν φ_ν
        for (mu, row) in p.iter().enumerate().take(n) {
            p_phi[mu] = row.iter().zip(phi.iter()).map(|(x, y)| x * y).sum();
        }
        let rho: f64 = phi.iter().zip(p_phi.iter()).map(|(x, y)| x * y).sum();
        if rho <= 1e-14 {
            continue;
        }
        let mut grad_rho = [0.0_f64; 3];
        if gradient {
            for (pp, dp) in p_phi.iter().zip(dphi.iter()) {
                for k in 0..3 {
                    grad_rho[k] += 2.0 * pp * dp[k];
                }
            }
        }
        let sigma =
            grad_rho[0] * grad_rho[0] + grad_rho[1] * grad_rho[1] + grad_rho[2] * grad_rho[2];
        let mut tau = 0.0;
        if meta {
            for (mu, row) in p.iter().enumerate().take(n) {
                for (nu, &pmn) in row.iter().enumerate().take(n) {
                    let dot = dphi[mu][0] * dphi[nu][0]
                        + dphi[mu][1] * dphi[nu][1]
                        + dphi[mu][2] * dphi[nu][2];
                    tau += 0.5 * pmn * dot;
                }
            }
        }
        let v = evaluate_xc(functional, exact_exchange, rho, sigma, tau);
        energy += w * v.energy_density;
        electrons += w * rho;
        // a_μ = w(½ v_ρ φ_μ + 2 v_σ ∇ρ·∇φ_μ); V_μν += a_μ φ_ν + φ_μ a_ν.
        for mu in 0..n {
            let grad_dot =
                grad_rho[0] * dphi[mu][0] + grad_rho[1] * dphi[mu][1] + grad_rho[2] * dphi[mu][2];
            a[mu] = w * (0.5 * v.v_rho * phi[mu] + 2.0 * v.v_sigma * grad_dot);
        }
        let half_tau = 0.5 * w * v.v_tau;
        for mu in 0..n {
            for nu in 0..=mu {
                let mut val = a[mu] * phi[nu] + phi[mu] * a[nu];
                if meta {
                    val += half_tau
                        * (dphi[mu][0] * dphi[nu][0]
                            + dphi[mu][1] * dphi[nu][1]
                            + dphi[mu][2] * dphi[nu][2]);
                }
                matrix[mu][nu] += val;
            }
        }
    }
    symmetrize_from_lower(&mut matrix);
    XcContribution {
        energy,
        matrix,
        electrons,
    }
}

/// Copy the lower triangle of a square matrix onto its upper triangle.
pub(super) fn symmetrize_from_lower(m: &mut [Vec<f64>]) {
    for i in 1..m.len() {
        let (upper_rows, lower_rows) = m.split_at_mut(i);
        if let Some(row_i) = lower_rows.first() {
            for (upper_row, &value) in upper_rows.iter_mut().zip(row_i.iter()) {
                if let Some(slot) = upper_row.get_mut(i) {
                    *slot = value;
                }
            }
        }
    }
}

/// `Tr(A B)` for symmetric matrices.
fn trace_product(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(ra, rb)| ra.iter().zip(rb.iter()).map(|(x, y)| x * y).sum::<f64>())
        .sum()
}

fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Matrix {
    let n = a.len();
    let m = b.first().map_or(0, Vec::len);
    let mut out = vec![vec![0.0_f64; m]; n];
    for (i, row) in a.iter().enumerate() {
        for (k, &aik) in row.iter().enumerate() {
            if aik == 0.0 {
                continue;
            }
            for (j, &bkj) in b[k].iter().enumerate() {
                out[i][j] += aik * bkj;
            }
        }
    }
    out
}

/// Solve the small dense system `A x = b` by Gaussian elimination with partial pivoting.
fn solve_linear(mut a: Matrix, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| {
            a[i][col]
                .abs()
                .partial_cmp(&a[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            if factor == 0.0 {
                continue;
            }
            let (top, bottom) = a.split_at_mut(row);
            if let (Some(pivot_row), Some(target)) = (top.get(col), bottom.first_mut()) {
                for (x, p) in target[col..].iter_mut().zip(pivot_row[col..].iter()) {
                    *x -= factor * p;
                }
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0_f64; n];
    for row in (0..n).rev() {
        let mut acc = b[row];
        for k in (row + 1)..n {
            acc -= a[row][k] * x[k];
        }
        x[row] = acc / a[row][row];
    }
    Some(x)
}

/// Pulay DIIS extrapolation of Fock matrices.
struct Diis {
    focks: Vec<Matrix>,
    errors: Vec<Vec<f64>>,
}

impl Diis {
    fn new() -> Self {
        Self {
            focks: Vec::new(),
            errors: Vec::new(),
        }
    }

    /// Store `(F, e)` and return the extrapolated Fock matrix.
    fn extrapolate(&mut self, fock: &Matrix, error: Vec<f64>) -> Matrix {
        if self.focks.len() == DIIS_DEPTH {
            self.focks.remove(0);
            self.errors.remove(0);
        }
        self.focks.push(fock.clone());
        self.errors.push(error);
        let m = self.focks.len();
        if m < 2 {
            return fock.clone();
        }
        let mut b = vec![vec![0.0_f64; m + 1]; m + 1];
        for (i, ei) in self.errors.iter().enumerate() {
            for (j, ej) in self.errors.iter().enumerate().take(i + 1) {
                let dot: f64 = ei.iter().zip(ej.iter()).map(|(x, y)| x * y).sum();
                b[i][j] = dot;
                b[j][i] = dot;
            }
            b[i][m] = -1.0;
            b[m][i] = -1.0;
        }
        // Scale for conditioning.
        let scale = (0..m).map(|i| b[i][i]).fold(0.0_f64, f64::max);
        if scale > 0.0 {
            for row in b.iter_mut().take(m) {
                for x in row.iter_mut().take(m) {
                    *x /= scale;
                }
            }
        }
        let mut rhs = vec![0.0_f64; m + 1];
        rhs[m] = -1.0;
        match solve_linear(b, rhs) {
            Some(c) => {
                let n = fock.len();
                let mut out = vec![vec![0.0_f64; n]; n];
                for (ci, f) in c.iter().take(m).zip(self.focks.iter()) {
                    for (orow, frow) in out.iter_mut().zip(f.iter()) {
                        for (o, x) in orow.iter_mut().zip(frow.iter()) {
                            *o += ci * x;
                        }
                    }
                }
                out
            }
            None => {
                // Singular subspace: restart from the newest vector.
                let newest_f = self.focks.pop();
                let newest_e = self.errors.pop();
                self.focks.clear();
                self.errors.clear();
                if let (Some(f), Some(e)) = (newest_f, newest_e) {
                    self.focks.push(f);
                    self.errors.push(e);
                }
                fock.clone()
            }
        }
    }
}

/// Run the restricted Kohn–Sham SCF for `dft`.
pub(crate) fn run_kohn_sham(
    dft: &DensityFunctionalTheory,
) -> Result<DftResult, QuantumChemistryError> {
    let hf: &HartreeFockSolver = &dft.hf_solver;
    let n = hf.n_basis;
    if hf.basis_functions.len() != n || n == 0 {
        return Err(QuantumChemistryError::MissingData("basis functions"));
    }
    if hf.overlap.len() != n || hf.core_hamiltonian.len() != n {
        return Err(QuantumChemistryError::DimensionMismatch {
            what: "one-electron matrices",
            expected: n,
            found: hf.overlap.len().min(hf.core_hamiltonian.len()),
        });
    }
    let n_pair = n * (n + 1) / 2;
    if hf.eri.len() != n_pair * (n_pair + 1) / 2 {
        return Err(QuantumChemistryError::DimensionMismatch {
            what: "packed AO ERIs",
            expected: n_pair * (n_pair + 1) / 2,
            found: hf.eri.len(),
        });
    }
    if !hf.n_electrons.is_multiple_of(2) {
        return Err(QuantumChemistryError::OpenShellReference(hf.n_electrons));
    }
    let n_occ = hf.n_electrons / 2;
    if n_occ > n {
        return Err(QuantumChemistryError::DimensionMismatch {
            what: "occupied orbitals vs basis functions",
            expected: n,
            found: n_occ,
        });
    }
    let grid = MolecularGrid::new(&hf.atomic_charges, &hf.atomic_positions, dft.grid)?;
    let a0 = dft.alpha_x;
    let h = &hf.core_hamiltonian;
    let s = &hf.overlap;
    let e_nn = hf.nuclear_repulsion();
    let diagonalise = |f: &[Vec<f64>]| -> Result<(Vec<f64>, Matrix), QuantumChemistryError> {
        generalized_symmetric_eigen_n(f, s).ok_or(QuantumChemistryError::SingularOverlap)
    };
    // Core-Hamiltonian guess.
    let (_, guess) = diagonalise(h)?;
    let mut density = HartreeFockSolver::build_density_from_coeffs(&guess, n_occ, n);
    let mut diis = Diis::new();
    let mut energy_old = f64::NAN;
    let mut status = ScfStatus::NotConverged;
    let mut iterations = 0;
    let mut final_parts = None;
    for iter in 0..hf.max_iter.max(1) {
        iterations = iter + 1;
        let j = hf.coulomb_matrix(&density);
        let k = if a0 != 0.0 {
            Some(hf.exchange_matrix(&density))
        } else {
            None
        };
        let xc = xc_contribution(&grid, &hf.basis_functions, &density, dft.xc, a0);
        let mut fock = h.clone();
        for mu in 0..n {
            for nu in 0..n {
                let mut val = j[mu][nu] + xc.matrix[mu][nu];
                if let Some(k) = &k {
                    val -= 0.5 * a0 * k[mu][nu];
                }
                fock[mu][nu] += val;
            }
        }
        let exact_exchange = k
            .as_ref()
            .map_or(0.0, |k| -0.25 * a0 * trace_product(&density, k));
        let energy = trace_product(&density, h)
            + 0.5 * trace_product(&density, &j)
            + exact_exchange
            + xc.energy
            + e_nn
            + dft.pseudopotential_energy;
        // Commutator FPS − SPF.
        let fps = mat_mul(&mat_mul(&fock, &density), s);
        let mut error = Vec::with_capacity(n * n);
        let mut error_max = 0.0_f64;
        for (mu, row) in fps.iter().enumerate() {
            for (nu, &x) in row.iter().enumerate() {
                let e = x - fps[nu][mu];
                error_max = error_max.max(e.abs());
                error.push(e);
            }
        }
        let converged = energy.is_finite()
            && (energy - energy_old).abs() < hf.energy_threshold
            && error_max < hf.density_threshold;
        if converged {
            status = ScfStatus::Converged;
            let (eps, coeffs) = diagonalise(&fock)?;
            final_parts = Some((energy, xc, exact_exchange, eps, coeffs));
            break;
        }
        energy_old = energy;
        let extrapolated = diis.extrapolate(&fock, error);
        let (eps, coeffs) = diagonalise(&extrapolated)?;
        if iter + 1 == hf.max_iter.max(1) {
            final_parts = Some((energy, xc, exact_exchange, eps, coeffs));
            break;
        }
        density = HartreeFockSolver::build_density_from_coeffs(&coeffs, n_occ, n);
    }
    let (energy, xc, exact_exchange, orbital_energies, mo_coefficients) =
        final_parts.ok_or(QuantumChemistryError::ScfNotConverged)?;
    let kinetic_energy = trace_product(&density, &hf.kinetic_matrix());
    Ok(DftResult {
        energy,
        xc_energy: xc.energy + exact_exchange,
        kinetic_energy,
        orbital_energies,
        density_matrix: density,
        status,
        mo_coefficients,
        iterations,
        integrated_electrons: xc.electrons,
    })
}
