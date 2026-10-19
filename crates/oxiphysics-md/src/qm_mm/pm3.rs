// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PM3 semi-empirical quantum chemistry engine.
//!
//! Implements the PM3 Hamiltonian (Stewart 1989, *J. Comput. Chem.* **10**, 209–220)
//! within the NDDO (Neglect of Diatomic Differential Overlap) approximation.
//!
//! # Key features
//! - STO-3G-contracted minimal basis (s on H; s, px, py, pz on C/N/O/F/S)
//! - Klopman-Ohno two-centre electron-repulsion integrals
//! - SCF loop solved via Löwdin-orthogonalised generalised eigenvalue problem
//!   using [`oxiphysics_core::numerical_methods::generalized_symmetric_eigen_n`]
//! - Core-core repulsion with PM3 Gaussian correction terms
//! - Numerical Hellmann-Feynman forces via ±0.001 Å finite differences
//!
//! # Supported elements
//! H (1), C (6), N (7), O (8), F (9), S (16).
//!
//! # References
//! Stewart, J. J. P. (1989). *J. Comput. Chem.* **10**, 209–220.

use super::QmRegion;
use super::embedding::{
    ANG_TO_BOHR, HARTREE_TO_KCAL, begin_run, external_potential, finite_difference_forces,
    forces_from_potential_derivatives,
};
use super::error::QmError;
use super::sp_overlap::{OrbType, build_overlap};
use crate::quantum_chemistry::{BasisFunction, CartesianAxis};
use oxiphysics_core::numerical_methods::generalized_symmetric_eigen_n;

// Type aliases to reduce type complexity in function signatures.
type BuildBasisResult =
    Result<(Vec<BasisFunction>, Vec<OrbType>, Vec<usize>, Vec<Pm3Params>), QmError>;
type ScfLoopResult = Result<(Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>), QmError>;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// 1 Hartree in eV.
const HARTREE_TO_EV: f64 = 27.211_386_02;
/// Finite-difference step for nuclear forces (Å).
const FD_STEP_ANG: f64 = 0.001;
/// Finite-difference step in the external potential φ_A (Hartree/e) used for the
/// forces on external charges.
const FD_STEP_PHI: f64 = 1e-3;
/// Maximum number of SCF iterations.
const MAX_SCF_ITER: usize = 500;

// ---------------------------------------------------------------------------
// PM3 parameters (Stewart 1989)
// ---------------------------------------------------------------------------

/// PM3 parameters for a single element.
///
/// All energetic parameters are stored in eV; Slater exponents in Bohr⁻¹.
#[derive(Clone, Copy)]
struct Pm3Params {
    /// Core orbital energy for s (eV).
    uss: f64,
    /// Core orbital energy for p (eV); 0 for H.
    upp: f64,
    /// Slater exponent for s orbital (Bohr⁻¹).
    zeta_s: f64,
    /// Slater exponent for p orbital (Bohr⁻¹); 0 for H.
    zeta_p: f64,
    /// Resonance integral β_s (eV).
    beta_s: f64,
    /// Resonance integral β_p (eV); 0 for H.
    beta_p: f64,
    /// Core-core repulsion α (Å⁻¹).
    alpha: f64,
    /// One-centre ss two-electron integral g_ss (eV).
    g_ss: f64,
    /// Number of valence electrons.
    n_val: usize,
    /// Gaussian correction pairs [(a_k, b_k, c_k in Å); 2].
    gaussians: [(f64, f64, f64); 2],
}

/// Return PM3 parameters for element `z`, or `None` if unsupported.
fn pm3_params(z: u32) -> Option<Pm3Params> {
    match z {
        1 => Some(Pm3Params {
            // H — Stewart 1989 Table I
            uss: -13.073_21,
            upp: 0.0,
            zeta_s: 1.268_64,
            zeta_p: 0.0,
            beta_s: -5.626_05,
            beta_p: 0.0,
            alpha: 3.356_10,
            g_ss: 14.794_2,
            n_val: 1,
            gaussians: [
                (1.128_89, 5.096_40, 1.236_60),
                (-1.060_43, 6.001_90, 1.893_20),
            ],
        }),
        6 => Some(Pm3Params {
            // C
            uss: -51.713_17,
            upp: -39.614_38,
            zeta_s: 1.565_38,
            zeta_p: 1.842_00,
            beta_s: -11.910_03,
            beta_p: -9.802_70,
            alpha: 2.707_90,
            g_ss: 11.200_7,
            n_val: 4,
            gaussians: [
                (0.051_20, 1.095_60, 2.026_00),
                (0.033_40, 2.071_50, 3.146_00),
            ],
        }),
        7 => Some(Pm3Params {
            // N
            uss: -71.932_13,
            upp: -57.172_30,
            zeta_s: 2.028_17,
            zeta_p: 2.313_60,
            beta_s: -14.062_51,
            beta_p: -20.043_81,
            alpha: 3.160_90,
            g_ss: 11.904_7,
            n_val: 5,
            gaussians: [
                (2.473_20, 3.246_60, 1.119_50),
                (-0.382_10, 2.003_90, 1.480_70),
            ],
        }),
        8 => Some(Pm3Params {
            // O
            uss: -99.644_36,
            upp: -77.797_50,
            zeta_s: 2.626_05,
            zeta_p: 2.576_23,
            beta_s: -29.272_77,
            beta_p: -29.272_77,
            alpha: 4.583_30,
            g_ss: 13.655_2,
            n_val: 6,
            gaussians: [
                (-1.139_65, 5.026_50, 0.847_63),
                (1.287_70, 7.002_30, 1.459_10),
            ],
        }),
        9 => Some(Pm3Params {
            // F
            uss: -136.105_68,
            upp: -104.889_54,
            zeta_s: 3.776_12,
            zeta_p: 2.494_67,
            beta_s: -48.291_69,
            beta_p: -36.508_38,
            alpha: 5.403_60,
            g_ss: 10.496_7,
            n_val: 7,
            gaussians: [
                (0.720_01, 4.829_50, 0.964_40),
                (-0.029_70, 1.555_80, 3.474_30),
            ],
        }),
        16 => Some(Pm3Params {
            // S
            uss: -55.578_64,
            upp: -48.173_80,
            zeta_s: 1.470_99,
            zeta_p: 1.725_90,
            beta_s: -7.930_50,
            beta_p: -6.816_01,
            alpha: 2.574_79,
            g_ss: 9.557_7,
            n_val: 6,
            gaussians: [
                (-0.358_50, 1.964_20, 0.387_50),
                (-0.013_80, 0.793_80, 3.565_00),
            ],
        }),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Basis construction
// ---------------------------------------------------------------------------

/// β resonance integral (eV) of an orbital of type `ot` on an atom with parameters `p`.
fn orbital_beta(ot: OrbType, p: &Pm3Params) -> f64 {
    match ot {
        OrbType::S => p.beta_s,
        OrbType::Px | OrbType::Py | OrbType::Pz => p.beta_p,
    }
}

/// On-site orbital energy (U_ss or U_pp) in Hartree.
fn orbital_u(ot: OrbType, p: &Pm3Params) -> f64 {
    match ot {
        OrbType::S => p.uss / HARTREE_TO_EV,
        OrbType::Px | OrbType::Py | OrbType::Pz => p.upp / HARTREE_TO_EV,
    }
}

/// Internal: build the STO-3G basis for the QM atoms.
///
/// Returns `(basis_functions, orb_types, atom_of_basis, params_per_basis)`.
/// `atom_of_basis[μ]` is the atom index to which basis function μ belongs.
fn build_basis(atomic_numbers: &[u32], positions: &[[f64; 3]]) -> BuildBasisResult {
    let mut basis: Vec<BasisFunction> = Vec::new();
    let mut orb_types: Vec<OrbType> = Vec::new();
    let mut atom_of_basis: Vec<usize> = Vec::new();
    let mut params_per_basis: Vec<Pm3Params> = Vec::new();

    for (iatom, (&z, pos_ang)) in atomic_numbers.iter().zip(positions.iter()).enumerate() {
        let p = pm3_params(z).ok_or(QmError::UnsupportedElement(z))?;
        // Convert Å position to Bohr for BasisFunction
        let center_bohr = [
            pos_ang[0] * ANG_TO_BOHR,
            pos_ang[1] * ANG_TO_BOHR,
            pos_ang[2] * ANG_TO_BOHR,
        ];

        // s orbital — always present
        basis.push(BasisFunction::sto3g_s(center_bohr, p.zeta_s));
        orb_types.push(OrbType::S);
        atom_of_basis.push(iatom);
        params_per_basis.push(p);

        // p orbitals — only for heavy atoms
        if p.zeta_p > 0.0 {
            for (ot, axis) in OrbType::P_SHELL.into_iter().zip(CartesianAxis::ALL) {
                basis.push(BasisFunction::sto3g_p(center_bohr, p.zeta_p, axis));
                orb_types.push(ot);
                atom_of_basis.push(iatom);
                params_per_basis.push(p);
            }
        }
    }

    Ok((basis, orb_types, atom_of_basis, params_per_basis))
}

// ---------------------------------------------------------------------------
// Klopman-Ohno two-centre electron-repulsion integral (Hartree)
// ---------------------------------------------------------------------------

/// Klopman-Ohno two-centre electron-repulsion integral.
///
/// γ(R, g_ss_A, g_ss_B) = 1 / sqrt(R² + a²)
/// where a = 1/g_ss_A + 1/g_ss_B  (Bohr, since g_ss in Hartree and we work in a.u.)
///
/// Returns value in Hartree.
#[inline]
fn klopman_ohno(r_bohr: f64, g_ss_a_hartree: f64, g_ss_b_hartree: f64) -> f64 {
    let a = 1.0 / g_ss_a_hartree + 1.0 / g_ss_b_hartree;
    1.0 / (r_bohr * r_bohr + a * a).sqrt()
}

/// Build the atom-atom γ matrix (Hartree).
fn build_gamma(positions: &[[f64; 3]], params: &[Pm3Params]) -> Vec<Vec<f64>> {
    let n = params.len().min(positions.len());
    let mut gamma = vec![vec![0.0_f64; n]; n];
    for a in 0..n {
        let ga = params[a].g_ss / HARTREE_TO_EV;
        gamma[a][a] = ga;
        for b in (a + 1)..n {
            let gb = params[b].g_ss / HARTREE_TO_EV;
            let dx = positions[a][0] - positions[b][0];
            let dy = positions[a][1] - positions[b][1];
            let dz = positions[a][2] - positions[b][2];
            let r_bohr = (dx * dx + dy * dy + dz * dz).sqrt() * ANG_TO_BOHR;
            let g = klopman_ohno(r_bohr, ga, gb);
            gamma[a][b] = g;
            gamma[b][a] = g;
        }
    }
    gamma
}

// ---------------------------------------------------------------------------
// Core Hamiltonian (NDDO)
// ---------------------------------------------------------------------------

/// Build the one-electron core Hamiltonian H_core (Hartree).
///
/// Diagonal:   H_μμ = U_μμ + Σ_{B≠A} (-Z_B * γ(R_AB))
/// Off-diag:   H_μν = 0.5 * (β_μ + β_ν) * S_μν   (two-centre, same-centre = 0)
///
/// Electrostatic embedding adds the external potential energy of an electron,
/// −φ_A at nucleus A, in the Mulliken form −½ S_μν (φ_A + φ_B).  This keeps the
/// generalised eigenproblem gauge invariant (a uniform φ shifts every orbital
/// energy by −φ and leaves the density unchanged) and reduces to −φ_A on the
/// diagonal.
fn build_h_core(
    orb_types: &[OrbType],
    atom_of_basis: &[usize],
    params_per_basis: &[Pm3Params],
    atom_params: &[Pm3Params],
    overlap: &[Vec<f64>],
    gamma: &[Vec<f64>],
    phi: &[f64],
) -> Vec<Vec<f64>> {
    let n = orb_types.len();
    let n_atoms = atom_params.len();
    let mut h = vec![vec![0.0_f64; n]; n];

    // Core charges Z_A = n_val
    let core_charge: Vec<f64> = atom_params.iter().map(|p| p.n_val as f64).collect();

    for mu in 0..n {
        let a = atom_of_basis[mu];
        let p_mu = &params_per_basis[mu];
        let ot_mu = orb_types[mu];

        for nu in 0..n {
            let b = atom_of_basis[nu];
            let p_nu = &params_per_basis[nu];
            let ot_nu = orb_types[nu];

            if a == b {
                // Same atom
                if mu == nu {
                    // Diagonal: U_μμ + Σ_{C≠A} (-Z_C * γ_AC)
                    let u_val = orbital_u(ot_mu, p_mu);
                    let mut v_sum = 0.0;
                    for c in 0..n_atoms {
                        if c != a {
                            v_sum -= core_charge[c] * gamma[a][c];
                        }
                    }
                    h[mu][nu] = u_val + v_sum;
                }
                // Off-diagonal same atom: 0 in NDDO (orthogonal AOs on same centre)
            } else {
                // Two-centre
                let beta_mu = orbital_beta(ot_mu, p_mu) / HARTREE_TO_EV;
                let beta_nu = orbital_beta(ot_nu, p_nu) / HARTREE_TO_EV;
                h[mu][nu] = 0.5 * (beta_mu + beta_nu) * overlap[mu][nu];
            }
            // Electrostatic embedding in the external point charges.
            h[mu][nu] -= 0.5 * overlap[mu][nu] * (phi[a] + phi[b]);
        }
    }
    h
}

// ---------------------------------------------------------------------------
// G matrix (NDDO monopole approximation)
// ---------------------------------------------------------------------------

/// Build the two-electron G matrix in the NDDO monopole approximation (Hartree).
///
/// Uses Mulliken atomic populations q_B = Σ_{ν∈B} (P·S)_νν to weight the
/// Coulomb terms, which avoids inconsistencies between the STO-3G overlap
/// and the ZDO approximation in the G-matrix.
///
/// **Diagonal (μ on A, μ==ν):**
/// G_μμ = (q_A - 0.5*(PS)_μμ) * γ_AA + Σ_{B≠A} q_B * γ_AB
///
/// **Off-diagonal same-atom:** 0
///
/// **Off-diagonal cross-atom (A ≠ B):**
/// G_μν = -0.5 * P_μν * γ_AB
fn build_g_matrix(
    density: &[Vec<f64>],
    overlap: &[Vec<f64>],
    atom_of_basis: &[usize],
    gamma: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n = density.len();
    let n_atoms = *atom_of_basis.iter().max().unwrap_or(&0) + 1;

    // Compute Mulliken orbital populations: (P·S)_μμ
    let ps_diag: Vec<f64> = (0..n)
        .map(|mu| {
            (0..n)
                .map(|nu| density[mu][nu] * overlap[nu][mu])
                .sum::<f64>()
        })
        .collect();

    // Mulliken atomic population: q_B = Σ_{ν ∈ B} (PS)_νν
    let mut atom_pop = vec![0.0_f64; n_atoms];
    for lam in 0..n {
        atom_pop[atom_of_basis[lam]] += ps_diag[lam];
    }

    let mut g = vec![vec![0.0_f64; n]; n];
    for mu in 0..n {
        let a = atom_of_basis[mu];
        for nu in 0..=mu {
            let b = atom_of_basis[nu];
            let g_val = if a == b {
                if mu == nu {
                    // Diagonal Coulomb + exchange-self
                    let j: f64 = (0..n_atoms).map(|c| atom_pop[c] * gamma[a][c]).sum();
                    let k_self = 0.5 * ps_diag[mu] * gamma[a][a];
                    j - k_self
                } else {
                    // Same-atom off-diagonal: zero in NDDO monopole
                    0.0_f64
                }
            } else {
                // Two-centre: exchange only
                -0.5 * density[mu][nu] * gamma[a][b]
            };
            g[mu][nu] = g_val;
            g[nu][mu] = g_val;
        }
    }
    g
}

// ---------------------------------------------------------------------------
// Density matrix from occupied MOs
// ---------------------------------------------------------------------------

/// Build density matrix P from occupied MO coefficients.
///
/// P_μν = 2 * Σ_{i=0}^{n_occ-1}  C_{μi} * C_{νi}
fn build_density(evecs: &[Vec<f64>], n_occ: usize, nbasis: usize) -> Vec<Vec<f64>> {
    let mut p = vec![vec![0.0_f64; nbasis]; nbasis];
    for ev in &evecs[..n_occ] {
        for mu in 0..nbasis {
            for nu in 0..nbasis {
                p[mu][nu] += 2.0 * ev[mu] * ev[nu];
            }
        }
    }
    p
}

// ---------------------------------------------------------------------------
// SCF loop
// ---------------------------------------------------------------------------

/// SCF loop: iterate Fock build → diagonalise → update density until convergence.
///
/// Uses linear density-matrix mixing (damping factor 0.4) to stabilise convergence.
/// Convergence is declared when |ΔE| < 1e-8 Hartree AND |ΔP|_max < 1e-6.
///
/// Returns `(density, orbital_energies_hartree, mo_coefficients)`.
fn scf_loop(
    h_core: &[Vec<f64>],
    overlap: &[Vec<f64>],
    n_occ: usize,
    atom_of_basis: &[usize],
    gamma: &[Vec<f64>],
    max_iter: usize,
) -> ScfLoopResult {
    let nbasis = h_core.len();

    // Initial guess: diagonalise core Hamiltonian
    let (_evals0, evecs0) =
        generalized_symmetric_eigen_n(h_core, overlap).ok_or(QmError::EigensolverFailed)?;

    let mut density = if n_occ > 0 {
        build_density(&evecs0, n_occ, nbasis)
    } else {
        vec![vec![0.0_f64; nbasis]; nbasis]
    };

    let mut evals: Vec<f64>;
    let mut evecs: Vec<Vec<f64>>;
    let mut e_prev = f64::MAX;

    // Linear mixing coefficient: new_density = α * P_new + (1-α) * P_old
    // α = 0.5 balances speed (faster convergence) vs stability (no oscillation).
    const ALPHA: f64 = 0.5;

    for iter in 0..max_iter {
        // Build Fock = H_core + G[P]
        let g = build_g_matrix(&density, overlap, atom_of_basis, gamma);
        let mut fock = vec![vec![0.0_f64; nbasis]; nbasis];
        for mu in 0..nbasis {
            for nu in 0..nbasis {
                fock[mu][nu] = h_core[mu][nu] + g[mu][nu];
            }
        }

        // Solve FC = SCε
        let (new_evals, new_evecs) =
            generalized_symmetric_eigen_n(&fock, overlap).ok_or(QmError::EigensolverFailed)?;

        // New density from current Fock eigenvectors
        let p_new = if n_occ > 0 {
            build_density(&new_evecs, n_occ, nbasis)
        } else {
            vec![vec![0.0_f64; nbasis]; nbasis]
        };

        // Linear mixing: reduce oscillation
        let mixed_density: Vec<Vec<f64>> = if iter == 0 {
            p_new.clone()
        } else {
            let mut mixed = vec![vec![0.0_f64; nbasis]; nbasis];
            for mu in 0..nbasis {
                for nu in 0..nbasis {
                    mixed[mu][nu] = ALPHA * p_new[mu][nu] + (1.0 - ALPHA) * density[mu][nu];
                }
            }
            mixed
        };

        // Electronic energy with the mixed density and the current Fock
        let mut e_elec = 0.0_f64;
        for mu in 0..nbasis {
            for nu in 0..nbasis {
                e_elec += mixed_density[mu][nu] * (h_core[mu][nu] + fock[mu][nu]);
            }
        }
        e_elec *= 0.5;

        // Convergence: compare density change and energy change
        let dp_max = {
            let mut max_val = 0.0_f64;
            for m in 0..nbasis {
                for n in 0..nbasis {
                    let d = (density[m][n] - mixed_density[m][n]).abs();
                    if d > max_val {
                        max_val = d;
                    }
                }
            }
            max_val
        };

        let de = (e_elec - e_prev).abs();
        e_prev = e_elec;
        density = mixed_density;
        evals = new_evals;
        evecs = new_evecs;

        if de < 1e-8 && dp_max < 1e-6 {
            return Ok((density, evals, evecs));
        }
    }

    Err(QmError::ScfNotConverged)
}

// ---------------------------------------------------------------------------
// Core-core repulsion
// ---------------------------------------------------------------------------

/// PM3 core-core repulsion energy (Hartree).
///
/// E_cc = Σ_{A<B} Z_A Z_B / R_AB * f(R_AB)
///
/// where f(R_AB) = 1 + exp(-α_A R) + exp(-α_B R)
///               + Σ_k [a_k^A exp(-b_k^A (R-c_k^A)²) + a_k^B exp(-b_k^B (R-c_k^B)²)]
fn core_core_repulsion(positions: &[[f64; 3]], params: &[Pm3Params]) -> f64 {
    let n = params.len().min(positions.len());
    let mut e_cc = 0.0_f64;

    for a in 0..n {
        let pa = &params[a];
        let za = pa.n_val as f64;

        for (b, pb) in params.iter().enumerate().take(n).skip(a + 1) {
            let zb = pb.n_val as f64;

            let dx = positions[a][0] - positions[b][0];
            let dy = positions[a][1] - positions[b][1];
            let dz = positions[a][2] - positions[b][2];
            let r_ang = (dx * dx + dy * dy + dz * dz).sqrt();
            if r_ang < 1e-10 {
                continue;
            }
            let r_bohr = r_ang * ANG_TO_BOHR;

            // Smooth switching function
            let mut f_r = 1.0 + (-pa.alpha * r_ang).exp() + (-pb.alpha * r_ang).exp();

            // Gaussian correction terms (PM3 specific)
            for &(ak, bk, ck) in &pa.gaussians {
                let delta = r_ang - ck;
                f_r += ak * (-bk * delta * delta).exp();
            }
            for &(ak, bk, ck) in &pb.gaussians {
                let delta = r_ang - ck;
                f_r += ak * (-bk * delta * delta).exp();
            }

            // Bare nuclear repulsion (Hartree): Z_A Z_B / R_AB
            e_cc += za * zb / r_bohr * f_r;
        }
    }
    e_cc
}

// ---------------------------------------------------------------------------
// Single-point PM3 energy (Hartree)
// ---------------------------------------------------------------------------

/// Result of a single PM3 SCF calculation.
struct Pm3Point {
    /// Total energy including the embedding interaction (Hartree).
    energy: f64,
    /// Electrostatic interaction with the external charges, Σ_A (Z_A − pop_A) φ_A (Hartree).
    embedding: f64,
    /// Mulliken net atomic charges Z_A − pop_A (e).
    net_charges: Vec<f64>,
}

/// Compute the PM3 total energy for the given geometry and external potential
/// φ_A (Hartree/e) at each nucleus.
fn single_point(
    atomic_numbers: &[u32],
    positions: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    phi: &[f64],
) -> Result<Pm3Point, QmError> {
    if atomic_numbers.is_empty() {
        return Err(QmError::EmptyRegion);
    }
    if positions.len() != atomic_numbers.len() || phi.len() != atomic_numbers.len() {
        return Err(QmError::DimensionMismatch {
            what: "PM3 positions / external potential vs atomic numbers",
            expected: atomic_numbers.len(),
            found: positions.len().min(phi.len()),
        });
    }
    let atom_params: Vec<Pm3Params> = atomic_numbers
        .iter()
        .map(|&z| pm3_params(z).ok_or(QmError::UnsupportedElement(z)))
        .collect::<Result<_, _>>()?;
    if multiplicity != 1 {
        return Err(QmError::OpenShellNotSupported);
    }

    let (basis, orb_types, atom_of_basis, params_per_basis) =
        build_basis(atomic_numbers, positions)?;
    let nbasis = basis.len();

    // Electron count (signed, so that an excessive positive charge cannot wrap).
    let n_val_total: i64 = atom_params.iter().map(|p| p.n_val as i64).sum();
    let n_electrons = n_val_total - i64::from(charge);
    if n_electrons < 0 || n_electrons > 2 * nbasis as i64 {
        return Err(QmError::InvalidElectronCount {
            electrons: n_electrons,
            n_basis: nbasis,
        });
    }
    if n_electrons % 2 != 0 {
        return Err(QmError::OpenShellNotSupported);
    }
    let n_occ = (n_electrons / 2) as usize;

    let overlap = build_overlap(&basis, &orb_types);
    let gamma = build_gamma(positions, &atom_params);
    let h_core = build_h_core(
        &orb_types,
        &atom_of_basis,
        &params_per_basis,
        &atom_params,
        &overlap,
        &gamma,
        phi,
    );

    let (density, _evals, _evecs) = scf_loop(
        &h_core,
        &overlap,
        n_occ,
        &atom_of_basis,
        &gamma,
        MAX_SCF_ITER,
    )?;

    // Electronic energy
    let g = build_g_matrix(&density, &overlap, &atom_of_basis, &gamma);
    let mut e_elec = 0.0_f64;
    for mu in 0..nbasis {
        for nu in 0..nbasis {
            let fock = h_core[mu][nu] + g[mu][nu];
            e_elec += density[mu][nu] * (h_core[mu][nu] + fock);
        }
    }
    e_elec *= 0.5;

    // Mulliken populations and net charges.
    let mut pop = vec![0.0_f64; atom_params.len()];
    for mu in 0..nbasis {
        let a = atom_of_basis[mu];
        for nu in 0..nbasis {
            pop[a] += density[mu][nu] * overlap[nu][mu];
        }
    }
    let net_charges: Vec<f64> = atom_params
        .iter()
        .zip(pop.iter())
        .map(|(p, &q)| p.n_val as f64 - q)
        .collect();

    // Core–external-charge interaction Σ_A Z_A φ_A (the electron part −Σ pop_A φ_A
    // is already contained in e_elec through H_core).
    let e_core_ext: f64 = atom_params
        .iter()
        .zip(phi.iter())
        .map(|(p, &f)| p.n_val as f64 * f)
        .sum();
    let embedding: f64 = net_charges.iter().zip(phi.iter()).map(|(q, f)| q * f).sum();

    let e_cc = core_core_repulsion(positions, &atom_params);
    Ok(Pm3Point {
        energy: e_elec + e_cc + e_core_ext,
        embedding,
        net_charges,
    })
}

/// Total PM3 energy (Hartree) at `positions` in the field of `external_charges`.
fn energy_at(region: &QmRegion, positions: &[[f64; 3]]) -> Result<f64, QmError> {
    let phi = external_potential(positions, &region.external_charges)?;
    single_point(
        &region.atomic_numbers,
        positions,
        region.charge,
        region.multiplicity,
        &phi,
    )
    .map(|p| p.energy)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Run a PM3 semi-empirical SCF calculation, storing energy and embedding
/// energy (kcal/mol), forces on the QM atoms and on the external charges
/// (kcal/mol/Å) and Mulliken charges into `region`.
///
/// Nuclear forces are central finite differences (±0.001 Å) of the total energy
/// in the field of the external charges.  Forces on the external charges use
/// the chain rule through the nuclear potentials, F_J = −Σ_A (∂E/∂φ_A) ∂φ_A/∂R_J,
/// with ∂E/∂φ_A from central differences in φ_A.
///
/// # Errors
/// Any failure (empty region, unsupported element, open shell, invalid electron
/// count, invalid external charge, eigensolver failure or SCF non-convergence,
/// including at a displaced geometry) is returned; the outputs of `region` are
/// then left cleared.
pub fn run_pm3(region: &mut QmRegion) -> Result<(), QmError> {
    let phi = begin_run(region)?;
    let point = single_point(
        &region.atomic_numbers,
        &region.positions,
        region.charge,
        region.multiplicity,
        &phi,
    )?;

    let forces =
        finite_difference_forces(&region.positions, FD_STEP_ANG, |pos| energy_at(region, pos))?;

    let external_forces = if region.external_charges.is_empty() {
        Vec::new()
    } else {
        let mut de_dphi = vec![0.0_f64; phi.len()];
        let mut phi_shift = phi.clone();
        for (a, slot) in de_dphi.iter_mut().enumerate() {
            phi_shift[a] = phi[a] + FD_STEP_PHI;
            let e_plus = single_point(
                &region.atomic_numbers,
                &region.positions,
                region.charge,
                region.multiplicity,
                &phi_shift,
            )?
            .energy;
            phi_shift[a] = phi[a] - FD_STEP_PHI;
            let e_minus = single_point(
                &region.atomic_numbers,
                &region.positions,
                region.charge,
                region.multiplicity,
                &phi_shift,
            )?
            .energy;
            phi_shift[a] = phi[a];
            *slot = (e_plus - e_minus) / (2.0 * FD_STEP_PHI);
        }
        forces_from_potential_derivatives(&region.positions, &region.external_charges, &de_dphi)
    };

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
    use crate::qm_mm::QmMethod;

    fn water(method: QmMethod) -> QmRegion {
        let mut region = QmRegion::new(vec![0, 1, 2], method, 0, 1);
        region.atomic_numbers = vec![8, 1, 1];
        region.set_positions(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.586],
            [0.0, -0.757, 0.586],
        ]);
        region
    }

    #[test]
    fn h2_pm3_energy_negative() {
        // H-H bond at 0.74 Å (experimental)
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_pm3(&mut region).expect("PM3 H2");
        assert!(
            region.energy < 0.0,
            "H2 PM3 energy should be negative: {}",
            region.energy
        );
        assert_eq!(region.embedding_energy, 0.0);
    }

    #[test]
    fn h2_forces_sum_to_zero() {
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_pm3(&mut region).expect("PM3 H2");
        let fx: f64 = region.forces.iter().map(|f| f[0]).sum();
        let fy: f64 = region.forces.iter().map(|f| f[1]).sum();
        let fz: f64 = region.forces.iter().map(|f| f[2]).sum();
        assert!(fx.abs() < 1e-3, "force x sum not zero: {fx}");
        assert!(fy.abs() < 1e-6, "force y sum not zero: {fy}");
        assert!(fz.abs() < 1e-6, "force z sum not zero: {fz}");
        // The bond force is along x and antisymmetric.
        assert!(region.forces[0][0].abs() > 1e-3);
    }

    #[test]
    fn h2o_pm3_converges() {
        let mut region = water(QmMethod::Pm3);
        run_pm3(&mut region).expect("PM3 water");
        assert!(region.energy.is_finite(), "water PM3 energy must be finite");
        assert!(region.energy < 0.0, "water PM3 energy should be negative");
        let q = &region.partial_charges;
        assert!(q.iter().all(|c| c.abs() < 1.0), "charges {q:?}");
        assert!((q[1] - q[2]).abs() < 1e-6, "symmetric H charges: {q:?}");
        assert!(q.iter().sum::<f64>().abs() < 1e-8, "neutral: {q:?}");
        // C2v symmetry: no dipole component perpendicular to the symmetry axis.
        let mu = region.dipole_moment();
        assert!(mu[0].abs() < 1e-8 && mu[1].abs() < 1e-6, "dipole {mu:?}");
    }

    #[test]
    fn unsupported_element_returns_error() {
        let mut region = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![92];
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(run_pm3(&mut region), Err(QmError::UnsupportedElement(92)));
    }

    #[test]
    fn failed_run_clears_stale_outputs() {
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_pm3(&mut region).expect("PM3 H2");
        assert!(region.energy != 0.0);
        // Now an unsupported element: the previous results must not survive.
        region.atomic_numbers = vec![1, 10];
        assert_eq!(run_pm3(&mut region), Err(QmError::UnsupportedElement(10)));
        assert_eq!(region.energy, 0.0);
        assert!(region.forces.iter().all(|f| *f == [0.0; 3]));
        assert!(region.partial_charges.iter().all(|&q| q == 0.0));
    }

    #[test]
    fn excessive_charge_is_an_error_not_a_panic() {
        for charge in [4, -4, 3] {
            let mut region = QmRegion::new(vec![0, 1], QmMethod::Pm3, charge, 1);
            region.atomic_numbers = vec![1, 1];
            region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
            let result = run_pm3(&mut region);
            match charge {
                3 => assert_eq!(
                    result,
                    Err(QmError::InvalidElectronCount {
                        electrons: -1,
                        n_basis: 2
                    })
                ),
                4 => assert_eq!(
                    result,
                    Err(QmError::InvalidElectronCount {
                        electrons: -2,
                        n_basis: 2
                    })
                ),
                _ => assert_eq!(
                    result,
                    Err(QmError::InvalidElectronCount {
                        electrons: 6,
                        n_basis: 2
                    })
                ),
            }
        }
    }

    #[test]
    fn zero_external_charge_reproduces_vacuum() {
        let mut vac = water(QmMethod::Pm3);
        run_pm3(&mut vac).expect("vacuum");
        let mut emb = water(QmMethod::Pm3);
        emb.external_charges = vec![([0.3, 0.2, -3.0], 0.0)];
        run_pm3(&mut emb).expect("embedded");
        assert!((vac.energy - emb.energy).abs() < 1e-9);
        assert_eq!(emb.embedding_energy, 0.0);
        assert_eq!(emb.external_charge_forces, vec![[0.0; 3]]);
    }

    #[test]
    fn h2o_dipole_grows_in_external_charge_field() {
        let mut vac = water(QmMethod::Pm3);
        run_pm3(&mut vac).expect("vacuum");
        let mu_vac = vac.dipole_moment()[2];

        // A positive charge below the oxygen creates a field along +z, which
        // pulls electron density towards −z and increases μ_z; the opposite
        // charge decreases it (positive polarisability along the field).
        let mut along = water(QmMethod::Pm3);
        along.external_charges = vec![([0.0, 0.0, -3.5], 1.0)];
        run_pm3(&mut along).expect("field along dipole");
        let mu_along = along.dipole_moment()[2];

        let mut against = water(QmMethod::Pm3);
        against.external_charges = vec![([0.0, 0.0, -3.5], -1.0)];
        run_pm3(&mut against).expect("field against dipole");
        let mu_against = against.dipole_moment()[2];

        assert!(
            mu_along > mu_vac + 0.01,
            "dipole must grow along the field: {mu_vac} -> {mu_along}"
        );
        assert!(
            mu_against < mu_vac - 0.01,
            "dipole must shrink against the field: {mu_vac} -> {mu_against}"
        );
        // Near-linear response: the two induced dipoles are nearly opposite.
        let induced_along = mu_along - mu_vac;
        let induced_against = mu_against - mu_vac;
        assert!(
            (induced_along + induced_against).abs() < 0.3 * induced_along.abs(),
            "induced dipoles {induced_along} / {induced_against}"
        );
    }

    #[test]
    fn external_charge_force_matches_finite_difference() {
        let q_pos = [0.4, 0.3, -3.0];
        let mut region = water(QmMethod::Pm3);
        region.external_charges = vec![(q_pos, 0.8)];
        run_pm3(&mut region).expect("embedded water");
        let f = region.external_charge_forces[0];

        let h = 1e-3;
        for k in 0..3 {
            let mut plus = water(QmMethod::Pm3);
            let mut p = q_pos;
            p[k] += h;
            plus.external_charges = vec![(p, 0.8)];
            run_pm3(&mut plus).expect("plus");
            let mut minus = water(QmMethod::Pm3);
            let mut m = q_pos;
            m[k] -= h;
            minus.external_charges = vec![(m, 0.8)];
            run_pm3(&mut minus).expect("minus");
            let fd = -(plus.total_energy_in_field() - minus.total_energy_in_field()) / (2.0 * h);
            assert!(
                (fd - f[k]).abs() < 2e-3 + 1e-2 * fd.abs(),
                "component {k}: analytic {} vs FD {fd}",
                f[k]
            );
        }

        // Newton's third law: the total force on QM atoms and charges vanishes.
        for k in 0..3 {
            let total: f64 = region.forces.iter().map(|g| g[k]).sum::<f64>() + f[k];
            assert!(total.abs() < 5e-3, "net force component {k}: {total}");
        }
    }
}
