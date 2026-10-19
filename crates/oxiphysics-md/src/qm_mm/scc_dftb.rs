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
//! ## Electrostatic embedding
//!
//! External point charges Q_J (from [`QmRegion::external_charges`]) create the
//! potential φ_A = Σ_J Q_J / R_AJ at every nucleus.  Following the standard
//! QM/MM-SCC-DFTB coupling (Cui et al. 2001) the energy gains
//! E_ext = Σ_A (q⁰_A − q_A) φ_A and the Hamiltonian the Mulliken-form shift
//! −½ S_μν (φ_A + φ_B).
//!
//! ## Reference
//! - Elstner, M. et al. (1998). *Phys. Rev. B* 58, 7260-7268.
//! - Cui, Q. et al. (2001). *J. Phys. Chem. B* 105, 569-585.
//! - Mataga, N., Nishimoto, K. (1957). *Z. Phys. Chem.* 13, 140.

use super::QmRegion;
use super::embedding::{
    ANG_TO_BOHR, HARTREE_TO_KCAL, begin_run, external_potential, finite_difference_forces,
    forces_from_potential_derivatives,
};
use super::error::QmError;
use super::sp_overlap::{OrbType, build_overlap};
use crate::quantum_chemistry::{BasisFunction, CartesianAxis};
use oxiphysics_core::numerical_methods::generalized_symmetric_eigen_n;

// ---------------------------------------------------------------------------
// Numerical constants
// ---------------------------------------------------------------------------

/// Wolfsberg–Helmholtz constant for the Extended-Hückel off-diagonal elements.
const K_EHT: f64 = 1.75;

/// Finite-difference step for numerical gradient (Å).
const GRAD_STEP_ANG: f64 = 0.001;

/// Maximum SCC iterations.
const MAX_ITER: usize = 500;

/// SCC linear charge-mixing coefficient.
const CHARGE_MIX: f64 = 0.30;

/// Convergence threshold on the Mulliken charges (e).
const CHARGE_TOL: f64 = 1e-9;

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
// Basis construction helpers
// ---------------------------------------------------------------------------

/// Build the DFTB minimal basis for all atoms.
///
/// H  → [1s]
/// Heavy atoms → [s, px, py, pz]
///
/// Returns `(basis, orbital types, atom_of_basis)`.
fn build_basis(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    params: &[DftbAtomParams],
) -> (Vec<BasisFunction>, Vec<OrbType>, Vec<usize>) {
    let mut basis: Vec<BasisFunction> = Vec::new();
    let mut orb_types: Vec<OrbType> = Vec::new();
    let mut atom_of_basis: Vec<usize> = Vec::new();

    for (a, (&z, pos)) in atomic_numbers.iter().zip(positions_ang.iter()).enumerate() {
        // Convert atom position to Bohr.
        let c = [
            pos[0] * ANG_TO_BOHR,
            pos[1] * ANG_TO_BOHR,
            pos[2] * ANG_TO_BOHR,
        ];
        let p = &params[a];

        basis.push(BasisFunction::sto3g_s(c, p.zeta_s));
        orb_types.push(OrbType::S);
        atom_of_basis.push(a);

        // Directional p orbitals for heavy atoms (Z ≥ 2).
        if z >= 2 {
            for (ot, axis) in OrbType::P_SHELL.into_iter().zip(CartesianAxis::ALL) {
                basis.push(BasisFunction::sto3g_p(c, p.zeta_p, axis));
                orb_types.push(ot);
                atom_of_basis.push(a);
            }
        }
    }

    (basis, orb_types, atom_of_basis)
}

// ---------------------------------------------------------------------------
// H⁰ matrix
// ---------------------------------------------------------------------------

/// Return the on-site energy (Hartree) of an orbital of type `ot` on an atom
/// with parameters `p`.
fn orbital_energy(ot: OrbType, p: &DftbAtomParams) -> f64 {
    match ot {
        OrbType::S => p.eps_s,
        OrbType::Px | OrbType::Py | OrbType::Pz => p.eps_p,
    }
}

/// Build the zeroth-order Hamiltonian H⁰ using the Extended-Hückel formula.
fn build_h0(
    orb_types: &[OrbType],
    overlap: &[Vec<f64>],
    atom_of_basis: &[usize],
    params: &[DftbAtomParams],
) -> Vec<Vec<f64>> {
    let nbas = orb_types.len();
    let mut h0 = vec![vec![0.0_f64; nbas]; nbas];

    for i in 0..nbas {
        let eps_i = orbital_energy(orb_types[i], &params[atom_of_basis[i]]);
        h0[i][i] = eps_i;
        for j in (i + 1)..nbas {
            let eps_j = orbital_energy(orb_types[j], &params[atom_of_basis[j]]);
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

/// Converged SCC state.
struct SccState {
    /// Converged density matrix P_μν (including the factor 2 of double occupation).
    density: Vec<Vec<f64>>,
    /// Mulliken charge fluctuations Δq_A = q_A − q⁰_A (electron population minus
    /// neutral valence count).
    dq: Vec<f64>,
    /// Number of SCC iterations used (inspected by the convergence tests).
    #[cfg(test)]
    iterations: usize,
}

/// Mulliken electron populations q_A = Σ_{μ∈A} Σ_ν P_μν S_νμ.
fn mulliken_populations(
    density: &[Vec<f64>],
    overlap: &[Vec<f64>],
    atom_of_basis: &[usize],
    n_atoms: usize,
) -> Vec<f64> {
    let mut pop = vec![0.0_f64; n_atoms];
    for (mu, &a) in atom_of_basis.iter().enumerate() {
        for (nu, p_row) in density[mu].iter().enumerate() {
            pop[a] += p_row * overlap[nu][mu];
        }
    }
    pop
}

/// Self-consistent-charge cycle.
///
/// H_μν = H⁰_μν + ½ S_μν (V_A + V_B) − ½ S_μν (φ_A + φ_B),
/// V_A = Σ_C γ_AC Δq_C,  Δq_C = q_C − q⁰_C,
///
/// so that an electron-rich atom (Δq > 0) raises its own orbital energies
/// (Elstner 1998, Eq. 14).
fn scc_loop(
    h0: &[Vec<f64>],
    overlap: &[Vec<f64>],
    gamma: &[Vec<f64>],
    atom_of_basis: &[usize],
    q0: &[f64],
    n_occ: usize,
    phi: &[f64],
) -> Result<SccState, QmError> {
    let nbas = h0.len();
    let n_atoms = gamma.len();

    let mut dq: Vec<f64> = vec![0.0_f64; n_atoms];

    let mut iterations = 0_usize;
    while iterations < MAX_ITER {
        iterations += 1;
        // 1. Build H^SCC.
        let v: Vec<f64> = (0..n_atoms)
            .map(|a| (0..n_atoms).map(|c| gamma[a][c] * dq[c]).sum::<f64>() - phi[a])
            .collect();
        let mut h_scc = h0.to_vec();
        for i in 0..nbas {
            let a = atom_of_basis[i];
            for j in 0..nbas {
                let b = atom_of_basis[j];
                h_scc[i][j] += 0.5 * overlap[i][j] * (v[a] + v[b]);
            }
        }

        // 2. Solve H^SCC C = ε S C.
        let (_evals, evecs) =
            generalized_symmetric_eigen_n(&h_scc, overlap).ok_or(QmError::EigensolverFailed)?;

        // 3. Density matrix P_μν = 2 Σ_{i<n_occ} C_μi C_νi (evecs[i][μ] = C_μi).
        let mut density = vec![vec![0.0_f64; nbas]; nbas];
        for evec in evecs.iter().take(n_occ) {
            for mu in 0..nbas {
                for nu in 0..nbas {
                    density[mu][nu] += 2.0 * evec[mu] * evec[nu];
                }
            }
        }

        // 4. New charge fluctuations.
        let pop = mulliken_populations(&density, overlap, atom_of_basis, n_atoms);
        let dq_out: Vec<f64> = pop.iter().zip(q0.iter()).map(|(p, q)| p - q).collect();

        let max_change = dq_out
            .iter()
            .zip(dq.iter())
            .map(|(o, i)| (o - i).abs())
            .fold(0.0_f64, f64::max);
        if max_change < CHARGE_TOL {
            return Ok(SccState {
                density,
                dq: dq_out,
                #[cfg(test)]
                iterations,
            });
        }

        // 5. Linear charge mixing.
        for (d, o) in dq.iter_mut().zip(dq_out.iter()) {
            *d = CHARGE_MIX * o + (1.0 - CHARGE_MIX) * *d;
        }
    }

    Err(QmError::ScfNotConverged)
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

/// Result of a single SCC-DFTB calculation.
struct DftbPoint {
    /// Total energy including the embedding interaction (Hartree).
    energy: f64,
    /// Electrostatic interaction with the external charges (Hartree).
    embedding: f64,
    /// Mulliken net atomic charges q⁰_A − q_A (e).
    net_charges: Vec<f64>,
    /// SCC iterations used (inspected by the convergence tests).
    #[cfg(test)]
    iterations: usize,
}

/// SCC-DFTB energy for the given geometry, total charge, multiplicity and
/// external potential φ_A (Hartree/e).
///
/// E = Tr(P H⁰) + ½ Σ_AB γ_AB Δq_A Δq_B − Σ_A Δq_A φ_A + E_rep
fn compute_point(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    phi: &[f64],
) -> Result<DftbPoint, QmError> {
    if atomic_numbers.is_empty() {
        return Err(QmError::EmptyRegion);
    }
    let n_atoms = atomic_numbers.len();
    if positions_ang.len() != n_atoms || phi.len() != n_atoms {
        return Err(QmError::DimensionMismatch {
            what: "SCC-DFTB positions / external potential vs atomic numbers",
            expected: n_atoms,
            found: positions_ang.len().min(phi.len()),
        });
    }

    let params: Vec<DftbAtomParams> = atomic_numbers
        .iter()
        .map(|&z| dftb_params(z).ok_or(QmError::UnsupportedElement(z)))
        .collect::<Result<_, _>>()?;
    if multiplicity != 1 {
        return Err(QmError::OpenShellNotSupported);
    }

    let (basis, orb_types, atom_of_basis) = build_basis(atomic_numbers, positions_ang, &params);
    let nbas = basis.len();

    // Electron count (signed so that a large positive charge cannot wrap).
    let total_val: i64 = params.iter().map(|p| p.n_val as i64).sum();
    let n_electrons = total_val - i64::from(charge);
    if n_electrons < 0 || n_electrons > 2 * nbas as i64 {
        return Err(QmError::InvalidElectronCount {
            electrons: n_electrons,
            n_basis: nbas,
        });
    }
    if n_electrons % 2 != 0 {
        return Err(QmError::OpenShellNotSupported);
    }
    let n_occ = (n_electrons / 2) as usize;

    let overlap = build_overlap(&basis, &orb_types);
    let h0 = build_h0(&orb_types, &overlap, &atom_of_basis, &params);
    let gamma = build_gamma(atomic_numbers, positions_ang, &params);
    let q0: Vec<f64> = params.iter().map(|p| p.n_val as f64).collect();

    let scc = scc_loop(&h0, &overlap, &gamma, &atom_of_basis, &q0, n_occ, phi)?;

    let e_h0: f64 = scc
        .density
        .iter()
        .zip(h0.iter())
        .map(|(p_row, h_row)| {
            p_row
                .iter()
                .zip(h_row.iter())
                .map(|(p, h)| p * h)
                .sum::<f64>()
        })
        .sum();
    let e_coulomb: f64 = gamma
        .iter()
        .zip(scc.dq.iter())
        .map(|(g_row, dq_a)| {
            0.5 * dq_a
                * g_row
                    .iter()
                    .zip(scc.dq.iter())
                    .map(|(g, dq_b)| g * dq_b)
                    .sum::<f64>()
        })
        .sum();
    let net_charges: Vec<f64> = scc.dq.iter().map(|d| -d).collect();
    let embedding: f64 = net_charges.iter().zip(phi.iter()).map(|(q, f)| q * f).sum();
    let e_rep = repulsive_energy(atomic_numbers, positions_ang);

    Ok(DftbPoint {
        energy: e_h0 + e_coulomb + embedding + e_rep,
        embedding,
        net_charges,
        #[cfg(test)]
        iterations: scc.iterations,
    })
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Run an SCC-DFTB calculation, storing energy and embedding energy (kcal/mol),
/// forces on the QM atoms and on the external charges (kcal/mol/Å) and Mulliken
/// charges in `region`.
///
/// Based on Elstner et al. (1998), *Phys. Rev. B* 58, 7260, with an
/// Extended-Hückel H⁰ and Mataga–Nishimoto γ_AB.  Nuclear forces are central
/// finite differences (±0.001 Å) of the total energy in the field of the
/// external charges.  The SCC energy is variational in the MO coefficients, so
/// the forces on the external charges follow from the Hellmann–Feynman theorem:
/// ∂E/∂φ_A = q⁰_A − q_A.
///
/// # Errors
/// Any failure (empty region, unsupported element, open shell, invalid electron
/// count, invalid external charge, eigensolver failure or SCC non-convergence,
/// including at a displaced geometry) is returned; the outputs of `region` are
/// then left cleared.
pub fn run_scc_dftb(region: &mut QmRegion) -> Result<(), QmError> {
    let phi = begin_run(region)?;
    let point = compute_point(
        &region.atomic_numbers,
        &region.positions,
        region.charge,
        region.multiplicity,
        &phi,
    )?;

    let forces = finite_difference_forces(&region.positions, GRAD_STEP_ANG, |pos| {
        let phi_disp = external_potential(pos, &region.external_charges)?;
        compute_point(
            &region.atomic_numbers,
            pos,
            region.charge,
            region.multiplicity,
            &phi_disp,
        )
        .map(|p| p.energy)
    })?;

    let external_forces = forces_from_potential_derivatives(
        &region.positions,
        &region.external_charges,
        &point.net_charges,
    );

    region.energy = (point.energy - point.embedding) * HARTREE_TO_KCAL;
    region.embedding_energy = point.embedding * HARTREE_TO_KCAL;
    region.forces = forces;
    region.partial_charges = point.net_charges;
    region.external_charge_forces = external_forces;
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qm_mm::{QmMethod, QmRegion};

    fn h2() -> QmRegion {
        let mut r = QmRegion::new(vec![0, 1], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![1, 1];
        r.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        r
    }

    fn water() -> QmRegion {
        let mut r = QmRegion::new(vec![0, 1, 2], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![8, 1, 1];
        r.set_positions(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.586],
            [0.0, -0.757, 0.586],
        ]);
        r
    }

    #[test]
    fn h2_dftb_energy_negative() {
        let mut r = h2();
        run_scc_dftb(&mut r).expect("DFTB H2");
        assert!(
            r.energy < 0.0,
            "H2 DFTB energy should be negative: {}",
            r.energy
        );
        // Homonuclear: no charge transfer.
        assert!(r.partial_charges.iter().all(|q| q.abs() < 1e-8));
    }

    #[test]
    fn h2_forces_sum_to_zero() {
        let mut r = h2();
        run_scc_dftb(&mut r).expect("DFTB H2");
        let f_sum: f64 = r.forces.iter().map(|f| f[0]).sum();
        assert!(f_sum.abs() < 1e-4, "force sum not zero: {f_sum}");
        assert!((r.forces[0][0] + r.forces[1][0]).abs() < 1e-4);
    }

    #[test]
    fn h2o_dftb_charges_have_correct_sign_and_converge_fast() {
        let r = water();
        let phi = vec![0.0; 3];
        let point =
            compute_point(&r.atomic_numbers, &r.positions, 0, 1, &phi).expect("water SCC-DFTB");
        let q = &point.net_charges;
        assert!(q[0] < 0.0 && q[0] > -1.0, "O charge {q:?}");
        assert!(q[1] > 0.0 && (q[1] - q[2]).abs() < 1e-8, "H charges {q:?}");
        assert!(q.iter().sum::<f64>().abs() < 1e-8);
        assert!(point.iterations < 50, "SCC iterations {}", point.iterations);
        assert!(point.energy.is_finite());
    }

    #[test]
    fn h2o_dftb_run_succeeds() {
        let mut r = water();
        run_scc_dftb(&mut r).expect("water DFTB");
        assert!(r.energy.is_finite() && r.energy != 0.0);
        // Symmetric hydrogens feel mirror-image forces.
        assert!((r.forces[1][1] + r.forces[2][1]).abs() < 1e-3);
        assert!((r.forces[1][2] - r.forces[2][2]).abs() < 1e-3);
    }

    #[test]
    fn unsupported_element_returns_error() {
        let mut r = QmRegion::new(vec![0], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![10]; // Ne — unsupported
        r.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(run_scc_dftb(&mut r), Err(QmError::UnsupportedElement(10)));
        assert_eq!(r.energy, 0.0);
    }

    #[test]
    fn open_shell_and_bad_charges_are_errors() {
        let mut r = QmRegion::new(vec![0], QmMethod::SccDftb, 0, 1);
        r.atomic_numbers = vec![1];
        r.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(run_scc_dftb(&mut r), Err(QmError::OpenShellNotSupported));

        let mut r = h2();
        r.charge = -4;
        assert_eq!(
            run_scc_dftb(&mut r),
            Err(QmError::InvalidElectronCount {
                electrons: 6,
                n_basis: 2
            })
        );
        r.charge = 4;
        assert_eq!(
            run_scc_dftb(&mut r),
            Err(QmError::InvalidElectronCount {
                electrons: -2,
                n_basis: 2
            })
        );
    }

    #[test]
    fn gamma_diagonal_equals_u() {
        // γ_AA(0) should equal U_A (self-interaction = Hubbard parameter).
        let u_c = 0.364_7_f64;
        let g = gamma_ab(0.0, u_c, u_c);
        assert!((g - u_c).abs() < 0.01, "γ_AA ≈ U_A: {g} vs {u_c}");
    }

    #[test]
    fn zero_external_charge_reproduces_vacuum() {
        let mut vac = water();
        run_scc_dftb(&mut vac).expect("vacuum");
        let mut emb = water();
        emb.external_charges = vec![([0.5, 0.0, -3.0], 0.0)];
        run_scc_dftb(&mut emb).expect("embedded");
        assert!((vac.energy - emb.energy).abs() < 1e-9);
        assert_eq!(emb.embedding_energy, 0.0);
    }

    #[test]
    fn h2_is_polarised_by_an_external_charge() {
        let mut r = h2();
        r.external_charges = vec![([3.0, 0.0, 0.0], 1.0)];
        run_scc_dftb(&mut r).expect("embedded H2");
        let q = &r.partial_charges;
        // Electrons are pulled towards the positive charge at +x.
        assert!(q[1] < -1e-3 && q[0] > 1e-3, "charges {q:?}");
        assert!((q[0] + q[1]).abs() < 1e-8);
        assert!(r.dipole_moment()[0] < 0.0);
        assert!(r.embedding_energy < 0.0);
    }

    #[test]
    fn h2o_dipole_grows_in_external_charge_field() {
        let mut vac = water();
        run_scc_dftb(&mut vac).expect("vacuum");
        let mu_vac = vac.dipole_moment()[2];
        assert!(mu_vac > 0.0);
        let mut along = water();
        along.external_charges = vec![([0.0, 0.0, -3.5], 1.0)];
        run_scc_dftb(&mut along).expect("field");
        let mu_along = along.dipole_moment()[2];
        assert!(
            mu_along > mu_vac + 0.01,
            "dipole must grow: {mu_vac} -> {mu_along}"
        );
    }

    #[test]
    fn energy_derivative_wrt_external_charge_is_qm_potential() {
        // Hellmann–Feynman: dE/dQ_J = Σ_A (q⁰_A − q_A) / R_AJ.
        let pos_j = [0.7, -0.4, -3.2];
        let q_j = 0.6;
        let run_with = |q: f64| -> QmRegion {
            let mut r = water();
            r.external_charges = vec![(pos_j, q)];
            run_scc_dftb(&mut r).expect("embedded water");
            r
        };
        let r0 = run_with(q_j);
        let h = 1e-4;
        let fd = (run_with(q_j + h).total_energy_in_field()
            - run_with(q_j - h).total_energy_in_field())
            / (2.0 * h);
        let analytic: f64 = r0
            .partial_charges
            .iter()
            .zip(r0.positions.iter())
            .map(|(q, p)| q / (crate::qm_mm::embedding::distance(p, &pos_j) * ANG_TO_BOHR))
            .sum::<f64>()
            * HARTREE_TO_KCAL;
        assert!(
            (fd - analytic).abs() < 1e-4 * analytic.abs().max(1.0),
            "dE/dQ: FD {fd} vs analytic {analytic}"
        );
    }

    #[test]
    fn external_charge_force_matches_finite_difference() {
        let q_pos = [0.4, 0.3, -3.0];
        let mut r = water();
        r.external_charges = vec![(q_pos, 0.8)];
        run_scc_dftb(&mut r).expect("embedded water");
        let f = r.external_charge_forces[0];
        let h = 1e-4;
        for k in 0..3 {
            let mut e = [0.0; 2];
            for (slot, sign) in e.iter_mut().zip([1.0, -1.0]) {
                let mut rr = water();
                let mut p = q_pos;
                p[k] += sign * h;
                rr.external_charges = vec![(p, 0.8)];
                run_scc_dftb(&mut rr).expect("displaced charge");
                *slot = rr.total_energy_in_field();
            }
            let fd = -(e[0] - e[1]) / (2.0 * h);
            assert!(
                (fd - f[k]).abs() < 1e-4 + 1e-4 * fd.abs(),
                "component {k}: analytic {} vs FD {fd}",
                f[k]
            );
        }
        for k in 0..3 {
            let total: f64 = r.forces.iter().map(|g| g[k]).sum::<f64>() + f[k];
            assert!(total.abs() < 1e-3, "net force component {k}: {total}");
        }
    }
}
