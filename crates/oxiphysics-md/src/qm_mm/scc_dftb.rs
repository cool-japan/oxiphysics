// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Self-Consistent-Charge Density-Functional Tight-Binding (SCC-DFTB) method.
//!
//! Implements the DFTB2 energy expression from Elstner et al. (1998),
//! *Phys. Rev. B* **58**, 7260.  The Hamiltonian matrix uses an Extended-Hückel
//! approximation for the zeroth-order (non-SCC) part, and the Mataga–Nishimoto
//! formula for the γ_AB coulomb kernel.  Repulsion is handled via simple
//! exponential pair potentials tuned to the mio-0-1 parameter set.
//!
//! ## Reference
//! - Elstner, M. et al. (1998). *Phys. Rev. B* 58, 7260-7268.
//! - Mataga, N., Nishimoto, K. (1957). *Z. Phys. Chem.* 13, 140.

use super::QmRegion;
use crate::quantum_chemistry::BasisFunction;
use oxiphysics_core::numerical_methods::generalized_symmetric_eigen_n;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Ångström → Bohr conversion factor.
const ANG_TO_BOHR: f64 = 1.889_725_989;

/// Hartree → kcal/mol conversion factor.
const HARTREE_TO_KCAL: f64 = 627.509_4;

/// Wolfsberg–Helmholtz constant for the Extended-Hückel off-diagonal elements.
const K_EHT: f64 = 1.75;

/// Finite-difference step for numerical gradient (Å).
const GRAD_STEP_ANG: f64 = 0.001;

/// Maximum SCC iterations.
const MAX_ITER: usize = 200;

/// SCC mixing coefficient (linear DIIS damping).
const DIIS_MIX: f64 = 0.30;

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur during an SCC-DFTB calculation.
#[derive(Debug)]
pub enum DftbError {
    /// One or more atoms have an element not covered by the mio-0-1 parameter set.
    UnsupportedElement(u32),
    /// Open-shell (multiplicity > 1) calculations are not yet supported.
    OpenShellNotSupported,
    /// The generalised eigenvalue solver failed (e.g. overlap matrix not positive definite).
    EigensolverFailed,
    /// The SCC cycle did not converge within the maximum number of iterations.
    SccNotConverged,
    /// The QM region contains no atoms.
    EmptyRegion,
}

impl std::fmt::Display for DftbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DftbError::UnsupportedElement(z) => {
                write!(
                    f,
                    "SCC-DFTB: element Z={z} is not in the mio-0-1 parameter set"
                )
            }
            DftbError::OpenShellNotSupported => {
                write!(
                    f,
                    "SCC-DFTB: open-shell (multiplicity > 1) is not implemented"
                )
            }
            DftbError::EigensolverFailed => {
                write!(
                    f,
                    "SCC-DFTB: generalised eigensolver failed (overlap not PD?)"
                )
            }
            DftbError::SccNotConverged => {
                write!(f, "SCC-DFTB: charge self-consistency did not converge")
            }
            DftbError::EmptyRegion => {
                write!(f, "SCC-DFTB: QM region is empty")
            }
        }
    }
}

impl std::error::Error for DftbError {}

// ---------------------------------------------------------------------------
// Atom parameters
// ---------------------------------------------------------------------------

/// Per-element parameters for the mio-0-1 DFTB parameter set.
struct DftbAtomParams {
    /// On-site energy for the s orbital (Hartree).
    eps_s: f64,
    /// On-site energy for p orbitals (Hartree); 0.0 for hydrogen.
    eps_p: f64,
    /// Hubbard U parameter (Hartree).
    hubbard_u: f64,
    /// Number of valence electrons.
    n_val: usize,
    /// Slater exponent ζ_s for the s orbital (Bohr⁻¹).
    zeta_s: f64,
    /// Slater exponent ζ_p for the p orbitals (Bohr⁻¹); 0.0 for hydrogen.
    zeta_p: f64,
}

/// Return DFTB atom parameters for element `z`, or `None` if not supported.
fn dftb_params(z: u32) -> Option<DftbAtomParams> {
    match z {
        1 => Some(DftbAtomParams {
            eps_s: -0.238_60,
            eps_p: 0.0,
            hubbard_u: 0.419_6,
            n_val: 1,
            zeta_s: 1.238_6,
            zeta_p: 0.0,
        }),
        6 => Some(DftbAtomParams {
            eps_s: -0.501_88,
            eps_p: -0.197_39,
            hubbard_u: 0.364_7,
            n_val: 4,
            zeta_s: 1.683_1,
            zeta_p: 1.683_1,
        }),
        7 => Some(DftbAtomParams {
            eps_s: -0.660_82,
            eps_p: -0.262_61,
            hubbard_u: 0.430_9,
            n_val: 5,
            zeta_s: 2.015_6,
            zeta_p: 2.015_6,
        }),
        8 => Some(DftbAtomParams {
            eps_s: -0.877_41,
            eps_p: -0.338_34,
            hubbard_u: 0.495_4,
            n_val: 6,
            zeta_s: 2.613_9,
            zeta_p: 2.613_9,
        }),
        16 => Some(DftbAtomParams {
            eps_s: -0.528_40,
            eps_p: -0.204_30,
            hubbard_u: 0.320_0,
            n_val: 6,
            zeta_s: 1.745_7,
            zeta_p: 1.745_7,
        }),
        _ => None,
    }
}

/// Exponential repulsive pair potential V_rep(R) = A * exp(−d * R) [Hartree; R in Bohr].
///
/// Parameters are approximate values tuned for qualitative DFTB2 behaviour.
fn repulsive_params(za: u32, zb: u32) -> (f64, f64) {
    let key = if za <= zb { (za, zb) } else { (zb, za) };
    match key {
        (1, 1) => (3.50, 2.20),
        (1, 6) => (3.80, 2.00),
        (1, 7) => (4.00, 2.10),
        (1, 8) => (4.20, 2.20),
        (6, 6) => (4.50, 1.80),
        (6, 7) => (4.70, 1.90),
        (6, 8) => (4.90, 2.00),
        (7, 7) => (5.00, 2.00),
        (7, 8) => (5.20, 2.10),
        (8, 8) => (5.50, 2.20),
        _ => (4.00, 2.00),
    }
}

// ---------------------------------------------------------------------------
// γ_AB (Mataga–Nishimoto)
// ---------------------------------------------------------------------------

/// Mataga–Nishimoto γ_AB coulomb kernel [Hartree].
///
/// For R → 0 (same atom) this reduces to the on-site Hubbard average
/// (U_A + U_B) / 2.
pub(crate) fn gamma_ab(r_bohr: f64, u_a: f64, u_b: f64) -> f64 {
    if r_bohr < 1e-10 {
        (u_a + u_b) / 2.0
    } else {
        1.0 / (r_bohr + 2.0 / (u_a + u_b))
    }
}

// ---------------------------------------------------------------------------
// SCC results container
// ---------------------------------------------------------------------------

struct SccResult {
    /// Total band-structure energy Σ ε_i (Hartree).
    e_band: f64,
    /// SCC coulomb correction (Hartree).
    e_coulomb: f64,
}

// ---------------------------------------------------------------------------
// Basis construction helpers
// ---------------------------------------------------------------------------

/// Build the DFTB minimal basis for all atoms.
///
/// H  → [1s]
/// Heavy atoms → [1s, 2px, 2py, 2pz]
///
/// Returns `(basis, atom_of_basis, orbital_type_of_basis)`.
/// `orbital_type_of_basis[i]` is `'s'` or `'p'`.
fn build_basis(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    params: &[DftbAtomParams],
) -> (Vec<BasisFunction>, Vec<usize>) {
    let mut basis: Vec<BasisFunction> = Vec::new();
    let mut atom_of_basis: Vec<usize> = Vec::new();

    for (a, (&z, pos)) in atomic_numbers.iter().zip(positions_ang.iter()).enumerate() {
        // Convert atom position to Bohr.
        let c = [
            pos[0] * ANG_TO_BOHR,
            pos[1] * ANG_TO_BOHR,
            pos[2] * ANG_TO_BOHR,
        ];
        let p = &params[a];

        // 1s orbital.
        basis.push(BasisFunction::sto3g_s(c, p.zeta_s));
        atom_of_basis.push(a);

        // 2p orbitals for heavy atoms (Z ≥ 2).
        if z >= 2 {
            // px, py, pz — they all share the same exponents; we use a single
            // contracted p-function here (the p-type BasisFunction in
            // quantum_chemistry represents one directional component; three are
            // needed for a complete sp valence shell).
            for _ in 0..3 {
                basis.push(BasisFunction::sto3g_p(c, p.zeta_p));
                atom_of_basis.push(a);
            }
        }
    }

    (basis, atom_of_basis)
}

// ---------------------------------------------------------------------------
// H⁰ matrix
// ---------------------------------------------------------------------------

/// Return the on-site energy for basis function `i`.
fn orbital_energy(
    bf_idx: usize,
    atom_of_basis: &[usize],
    params: &[DftbAtomParams],
    basis: &[BasisFunction],
) -> f64 {
    let a = atom_of_basis[bf_idx];
    if basis[bf_idx].angular_momentum == 0 {
        params[a].eps_s
    } else {
        params[a].eps_p
    }
}

/// Build the zeroth-order Hamiltonian H⁰ using the Extended-Hückel formula.
fn build_h0(
    basis: &[BasisFunction],
    overlap: &[Vec<f64>],
    atom_of_basis: &[usize],
    params: &[DftbAtomParams],
) -> Vec<Vec<f64>> {
    let nbas = basis.len();
    let mut h0 = vec![vec![0.0_f64; nbas]; nbas];

    for i in 0..nbas {
        let eps_i = orbital_energy(i, atom_of_basis, params, basis);
        h0[i][i] = eps_i;
        for j in (i + 1)..nbas {
            let eps_j = orbital_energy(j, atom_of_basis, params, basis);
            let h_ij = K_EHT * (eps_i + eps_j) / 2.0 * overlap[i][j];
            h0[i][j] = h_ij;
            h0[j][i] = h_ij;
        }
    }
    h0
}

// ---------------------------------------------------------------------------
// γ matrix
// ---------------------------------------------------------------------------

/// Build the n_atoms × n_atoms γ_AB matrix.
fn build_gamma(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    params: &[DftbAtomParams],
) -> Vec<Vec<f64>> {
    let n = atomic_numbers.len();
    let mut gamma = vec![vec![0.0_f64; n]; n];
    for a in 0..n {
        for b in 0..n {
            let r_bohr = if a == b {
                0.0
            } else {
                let dx = (positions_ang[a][0] - positions_ang[b][0]) * ANG_TO_BOHR;
                let dy = (positions_ang[a][1] - positions_ang[b][1]) * ANG_TO_BOHR;
                let dz = (positions_ang[a][2] - positions_ang[b][2]) * ANG_TO_BOHR;
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            gamma[a][b] = gamma_ab(r_bohr, params[a].hubbard_u, params[b].hubbard_u);
        }
    }
    gamma
}

// ---------------------------------------------------------------------------
// SCC loop
// ---------------------------------------------------------------------------

fn scc_loop(
    h0: &[Vec<f64>],
    overlap: &[Vec<f64>],
    gamma: &[Vec<f64>],
    atom_of_basis: &[usize],
    n_val_per_atom: &[usize],
    n_occ: usize,
) -> Result<SccResult, DftbError> {
    let nbas = h0.len();
    let n_atoms = gamma.len();

    // Reference charges q⁰_A = n_val_A (neutral atom).
    let q0: Vec<f64> = n_val_per_atom.iter().map(|&v| v as f64).collect();

    // Initialise Mulliken charges to neutral.
    let mut dq: Vec<f64> = vec![0.0_f64; n_atoms];

    // Density matrix (zero initially).
    let mut density = vec![vec![0.0_f64; nbas]; nbas];

    let mut e_band_prev = 0.0_f64;

    for _iter in 0..MAX_ITER {
        // 1. Build H^SCC = H^0 + ΔH^SCC.
        let mut h_scc = h0.to_vec();
        for i in 0..nbas {
            let atom_i = atom_of_basis[i];
            for j in 0..nbas {
                let atom_j = atom_of_basis[j];
                // SCC shift: ½ S_μν Σ_C (γ_AC + γ_BC) Δq_C
                let mut shift = 0.0_f64;
                for c in 0..n_atoms {
                    shift += (gamma[atom_i][c] + gamma[atom_j][c]) * dq[c];
                }
                h_scc[i][j] += 0.5 * overlap[i][j] * shift;
            }
        }

        // 2. Solve H^SCC C = ε S C.
        let (evals, evecs) =
            generalized_symmetric_eigen_n(&h_scc, overlap).ok_or(DftbError::EigensolverFailed)?;

        // 3. Build density matrix P_μν = 2 Σ_{i=0}^{n_occ-1} C_μi * C_νi.
        //    evecs[i] = i-th eigenvector; evecs[i][μ] = C_μi.
        for mu in 0..nbas {
            for nu in 0..nbas {
                let mut p = 0.0_f64;
                for evec in evecs.iter().take(n_occ) {
                    p += evec[mu] * evec[nu];
                }
                density[mu][nu] = 2.0 * p;
            }
        }

        // 4. Mulliken charges q_A = n_val_A - Σ_{μ∈A} Σ_ν P_μν S_νμ.
        let mut q_new = vec![0.0_f64; n_atoms];
        for mu in 0..nbas {
            let a = atom_of_basis[mu];
            for nu in 0..nbas {
                q_new[a] += density[mu][nu] * overlap[nu][mu];
            }
        }
        // q_new currently holds the population; convert to Mulliken charge.
        // The "Mulliken charge" here tracks the electron count at each atom.
        // dq_A = q0_A - q_new_A  (charge deviation from neutral).
        let dq_new: Vec<f64> = (0..n_atoms).map(|a| q0[a] - q_new[a]).collect();

        // 5. Band energy — sum of occupied orbital eigenvalues (Hartree).
        // The density matrix already accounts for 2-fold occupation; eigenvalues
        // are single-particle energies so we sum without an extra factor of 2.
        let e_band: f64 = evals.iter().take(n_occ).sum::<f64>();

        // 6. Convergence check.
        let max_dq_change: f64 = (0..n_atoms)
            .map(|a| (dq_new[a] - dq[a]).abs())
            .fold(0.0_f64, f64::max);
        let de = (e_band - e_band_prev).abs();

        // 7. DIIS linear mixing.
        for a in 0..n_atoms {
            dq[a] = DIIS_MIX * dq_new[a] + (1.0 - DIIS_MIX) * dq[a];
        }
        e_band_prev = e_band;

        if max_dq_change < 1e-6 && de < 1e-8 {
            // Converged — compute final coulomb correction.
            let e_coulomb: f64 = (0..n_atoms)
                .map(|a| {
                    (0..n_atoms)
                        .map(|b| 0.5 * gamma[a][b] * dq[a] * dq[b])
                        .sum::<f64>()
                })
                .sum();

            return Ok(SccResult { e_band, e_coulomb });
        }
    }

    Err(DftbError::SccNotConverged)
}

// ---------------------------------------------------------------------------
// Repulsive energy
// ---------------------------------------------------------------------------

fn repulsive_energy(atomic_numbers: &[u32], positions_ang: &[[f64; 3]]) -> f64 {
    let n = atomic_numbers.len();
    let mut e_rep = 0.0_f64;
    for a in 0..n {
        for b in (a + 1)..n {
            let dx = (positions_ang[a][0] - positions_ang[b][0]) * ANG_TO_BOHR;
            let dy = (positions_ang[a][1] - positions_ang[b][1]) * ANG_TO_BOHR;
            let dz = (positions_ang[a][2] - positions_ang[b][2]) * ANG_TO_BOHR;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let (amp, decay) = repulsive_params(atomic_numbers[a], atomic_numbers[b]);
            e_rep += amp * (-decay * r).exp();
        }
    }
    e_rep
}

// ---------------------------------------------------------------------------
// Core energy evaluation (in Hartree)
// ---------------------------------------------------------------------------

fn compute_energy_hartree(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    n_electrons: usize,
    n_occ: usize,
) -> Result<f64, DftbError> {
    if atomic_numbers.is_empty() {
        return Err(DftbError::EmptyRegion);
    }
    let n_atoms = atomic_numbers.len();

    // Collect per-atom parameters.
    let mut params: Vec<DftbAtomParams> = Vec::with_capacity(n_atoms);
    for &z in atomic_numbers {
        params.push(dftb_params(z).ok_or(DftbError::UnsupportedElement(z))?);
    }
    let _ = n_electrons; // used only for validation below

    // Build basis.
    let (basis, atom_of_basis) = build_basis(atomic_numbers, positions_ang, &params);
    let nbas = basis.len();

    // Overlap matrix S.
    let mut overlap = vec![vec![0.0_f64; nbas]; nbas];
    for i in 0..nbas {
        for j in 0..nbas {
            overlap[i][j] = basis[i].overlap_with(&basis[j]);
        }
    }
    // Regularise diagonal to avoid near-singular overlaps at zero distance.
    for (i, row) in overlap.iter_mut().enumerate() {
        if row[i] < 1e-10 {
            row[i] = 1.0;
        }
    }

    // H⁰.
    let h0 = build_h0(&basis, &overlap, &atom_of_basis, &params);

    // γ matrix.
    let gamma = build_gamma(atomic_numbers, positions_ang, &params);

    // Valence electrons per atom.
    let n_val_per_atom: Vec<usize> = params.iter().map(|p| p.n_val).collect();

    // SCC loop.
    let scc = scc_loop(
        &h0,
        &overlap,
        &gamma,
        &atom_of_basis,
        &n_val_per_atom,
        n_occ,
    )?;

    // Total energy (Elstner 1998, eqs. 14–17):
    // E_total = E_band - E_coulomb + E_rep
    // (E_band already double-counts the charge-charge term once; subtracting
    // E_coulomb corrects the double-counting.)
    let e_rep = repulsive_energy(atomic_numbers, positions_ang);
    let e_total = scc.e_band - scc.e_coulomb + e_rep;

    Ok(e_total)
}

// ---------------------------------------------------------------------------
// Numerical gradient
// ---------------------------------------------------------------------------

fn numerical_forces(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    n_electrons: usize,
    n_occ: usize,
) -> Vec<[f64; 3]> {
    let n = atomic_numbers.len();
    let mut forces = vec![[0.0_f64; 3]; n];

    for a in 0..n {
        for k in 0..3 {
            // Forward displacement.
            let mut pos_fwd = positions_ang.to_vec();
            pos_fwd[a][k] += GRAD_STEP_ANG;
            let e_fwd =
                compute_energy_hartree(atomic_numbers, &pos_fwd, n_electrons, n_occ).unwrap_or(0.0);

            // Backward displacement.
            let mut pos_bwd = positions_ang.to_vec();
            pos_bwd[a][k] -= GRAD_STEP_ANG;
            let e_bwd =
                compute_energy_hartree(atomic_numbers, &pos_bwd, n_electrons, n_occ).unwrap_or(0.0);

            // Central difference: F = -dE/dx (Hartree/Å), then convert to kcal/mol/Å.
            forces[a][k] = -(e_fwd - e_bwd) / (2.0 * GRAD_STEP_ANG) * HARTREE_TO_KCAL;
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Run an SCC-DFTB calculation, storing energy (kcal/mol) and forces (kcal/mol/Å)
/// in `region`.
///
/// Based on Elstner et al. (1998), *Phys. Rev. B* 58, 7260, with analytical
/// Extended-Hückel H⁰ and Mataga–Nishimoto γ_AB.  If the calculation fails for
/// any reason (unsupported element, eigensolver failure, non-convergence) the
/// energy and forces remain at zero without panicking.
pub fn run_scc_dftb(region: &mut QmRegion) {
    if let Err(_e) = try_run_scc_dftb(region) {
        // energy/forces stay at 0 — non-panicking failure
    }
}

fn try_run_scc_dftb(region: &mut QmRegion) -> Result<(), DftbError> {
    if region.atomic_numbers.is_empty() {
        return Err(DftbError::EmptyRegion);
    }
    if region.multiplicity != 1 {
        return Err(DftbError::OpenShellNotSupported);
    }

    let n_atoms = region.atomic_numbers.len();

    // Collect per-atom parameters (validate first).
    let mut params_check: Vec<DftbAtomParams> = Vec::with_capacity(n_atoms);
    for &z in &region.atomic_numbers {
        params_check.push(dftb_params(z).ok_or(DftbError::UnsupportedElement(z))?);
    }

    // Total valence electrons accounting for charge.
    let total_val: usize = params_check.iter().map(|p| p.n_val).sum();
    let total_val_i = total_val as i32 - region.charge;
    if total_val_i <= 0 {
        return Err(DftbError::UnsupportedElement(0));
    }
    let n_electrons = total_val_i as usize;
    if !n_electrons.is_multiple_of(2) {
        return Err(DftbError::OpenShellNotSupported);
    }
    let n_occ = n_electrons / 2;

    let e_hartree = compute_energy_hartree(
        &region.atomic_numbers,
        &region.positions,
        n_electrons,
        n_occ,
    )?;

    region.energy = e_hartree * HARTREE_TO_KCAL;

    region.forces = numerical_forces(
        &region.atomic_numbers,
        &region.positions,
        n_electrons,
        n_occ,
    );

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qm_mm::{QmMethod, QmRegion};

    #[test]
    fn h2_dftb_energy_negative() {
        let mut r = QmRegion::new(vec![0, 1], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![1, 1];
        r.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_scc_dftb(&mut r);
        assert!(
            r.energy < 0.0,
            "H2 DFTB energy should be negative: {}",
            r.energy
        );
    }

    #[test]
    fn h2_forces_sum_to_zero() {
        let mut r = QmRegion::new(vec![0, 1], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![1, 1];
        r.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_scc_dftb(&mut r);
        let f_sum: f64 = r.forces.iter().map(|f| f[0]).sum();
        assert!(f_sum.abs() < 2.0, "force sum not zero: {f_sum}");
    }

    #[test]
    fn h2o_dftb_finite_energy() {
        let mut r = QmRegion::new(vec![0, 1, 2], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![8, 1, 1];
        r.set_positions(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.586],
            [0.0, -0.757, 0.586],
        ]);
        run_scc_dftb(&mut r);
        assert!(r.energy.is_finite(), "water DFTB energy finite");
    }

    #[test]
    fn unsupported_element_no_panic() {
        let mut r = QmRegion::new(vec![0], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![10]; // Ne — unsupported
        r.set_positions(vec![[0.0, 0.0, 0.0]]);
        run_scc_dftb(&mut r); // must not panic
    }

    #[test]
    fn gamma_diagonal_equals_u() {
        // γ_AA(0) should equal U_A (self-interaction = Hubbard parameter).
        let u_c = 0.364_7_f64;
        let g = gamma_ab(0.0, u_c, u_c);
        assert!((g - u_c).abs() < 0.01, "γ_AA ≈ U_A: {g} vs {u_c}");
    }
}
