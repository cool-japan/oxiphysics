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
use crate::quantum_chemistry::BasisFunction;
use oxiphysics_core::numerical_methods::generalized_symmetric_eigen_n;

// Type aliases to reduce type complexity in function signatures.
type BuildBasisResult =
    Result<(Vec<BasisFunction>, Vec<OrbType>, Vec<usize>, Vec<Pm3Params>), Pm3Error>;
type ScfLoopResult = Result<(Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>), Pm3Error>;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Ångström → Bohr (1 Å = 1/0.529177 Bohr).
const ANG_TO_BOHR: f64 = 1.889_726_125_5;
/// 1 Hartree in eV.
const HARTREE_TO_EV: f64 = 27.211_386_02;
/// 1 Hartree in kcal/mol.
const HARTREE_TO_KCAL: f64 = 627.509_474;

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur during a PM3 calculation.
#[derive(Debug)]
pub enum Pm3Error {
    /// One or more atoms have an atomic number with no PM3 parameters.
    UnsupportedElement(u32),
    /// Open-shell (multiplicity ≠ 1) systems are not supported by this implementation.
    OpenShellNotSupported,
    /// The generalised symmetric eigensolver returned `None` (overlap not positive definite).
    EigensolverFailed,
    /// SCF did not converge within the maximum number of iterations.
    ScfNotConverged,
    /// The QM region contains no atoms.
    EmptyRegion,
}

impl std::fmt::Display for Pm3Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedElement(z) => {
                write!(f, "PM3: no parameters for element Z={z}")
            }
            Self::OpenShellNotSupported => {
                write!(
                    f,
                    "PM3: open-shell systems (multiplicity ≠ 1) are not supported"
                )
            }
            Self::EigensolverFailed => {
                write!(
                    f,
                    "PM3: generalised eigensolver failed (overlap matrix not positive definite?)"
                )
            }
            Self::ScfNotConverged => {
                write!(f, "PM3: SCF did not converge")
            }
            Self::EmptyRegion => {
                write!(f, "PM3: QM region contains no atoms")
            }
        }
    }
}

impl std::error::Error for Pm3Error {}

// ---------------------------------------------------------------------------
// PM3 parameters (Stewart 1989)
// ---------------------------------------------------------------------------

/// PM3 parameters for a single element.
///
/// All energetic parameters are stored in eV; Slater exponents in Bohr⁻¹.
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

/// Orbital angular momentum type for a single basis function slot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OrbType {
    /// s-type orbital.
    S,
    /// px-type orbital.
    Px,
    /// py-type orbital.
    Py,
    /// pz-type orbital.
    Pz,
}

impl OrbType {
    /// Return the β resonance integral (in eV) from the atom's PM3 parameters.
    fn beta(self, p: &Pm3Params) -> f64 {
        match self {
            OrbType::S => p.beta_s,
            OrbType::Px | OrbType::Py | OrbType::Pz => p.beta_p,
        }
    }

    /// Return the on-site orbital energy (U_ss or U_pp) in Hartree.
    fn u_energy(self, p: &Pm3Params) -> f64 {
        match self {
            OrbType::S => p.uss / HARTREE_TO_EV,
            OrbType::Px | OrbType::Py | OrbType::Pz => p.upp / HARTREE_TO_EV,
        }
    }
}

/// Internal: build the STO-3G basis for the QM region.
///
/// Returns `(basis_functions, orb_types, atom_of_basis, params_per_basis)`.
/// `atom_of_basis[μ]` is the atom index to which basis function μ belongs.
fn build_basis(region: &QmRegion) -> BuildBasisResult {
    let mut basis: Vec<BasisFunction> = Vec::new();
    let mut orb_types: Vec<OrbType> = Vec::new();
    let mut atom_of_basis: Vec<usize> = Vec::new();
    let mut params_per_basis: Vec<Pm3Params> = Vec::new();

    for (iatom, &z) in region.atomic_numbers.iter().enumerate() {
        let p = pm3_params(z).ok_or(Pm3Error::UnsupportedElement(z))?;
        // Convert Å position to Bohr for BasisFunction
        let pos_ang = region.positions[iatom];
        let center_bohr = [
            pos_ang[0] * ANG_TO_BOHR,
            pos_ang[1] * ANG_TO_BOHR,
            pos_ang[2] * ANG_TO_BOHR,
        ];

        // s orbital — always present
        // Re-fetch params for owned storage (pm3_params does not implement Clone).
        let p_owned = pm3_params(z).ok_or(Pm3Error::UnsupportedElement(z))?;
        basis.push(BasisFunction::sto3g_s(center_bohr, p.zeta_s));
        orb_types.push(OrbType::S);
        atom_of_basis.push(iatom);
        params_per_basis.push(p_owned);

        // p orbitals — only for heavy atoms
        if p.zeta_p > 0.0 {
            for ot in &[OrbType::Px, OrbType::Py, OrbType::Pz] {
                let pp = pm3_params(z).ok_or(Pm3Error::UnsupportedElement(z))?;
                basis.push(BasisFunction::sto3g_p(center_bohr, p.zeta_p));
                orb_types.push(*ot);
                atom_of_basis.push(iatom);
                params_per_basis.push(pp);
            }
        }
    }

    Ok((basis, orb_types, atom_of_basis, params_per_basis))
}

// ---------------------------------------------------------------------------
// Overlap matrix — normalized STO-3G contractions
// ---------------------------------------------------------------------------

/// Primitive Gaussian normalization constant for an s-type orbital.
///
/// N(α) = (2α/π)^(3/4)
#[inline]
fn norm_s(alpha: f64) -> f64 {
    use std::f64::consts::PI;
    (2.0 * alpha / PI).powf(0.75)
}

/// Primitive Gaussian normalization constant for a p-type orbital.
///
/// N(α) = (128 α^5 / π^3)^(1/4)
/// Equivalently: (2α/π)^(3/4) * (4α)^(1/2) / π^(1/4) · (π/4α)^(1/4)
/// Standard form: N_p(α) = (2α/π)^(3/4) * (4α)^(1/2) * π^(-1/4)
#[inline]
fn norm_p(alpha: f64) -> f64 {
    use std::f64::consts::PI;
    // N_p = (128 α^5 / π^3)^(1/4)
    (128.0 * alpha.powi(5) / PI.powi(3)).powf(0.25)
}

/// Contracted Gaussian primitive overlap for a specific (l_i, l_j, axis) combination.
///
/// `axis` is 0=x, 1=y, 2=z for the p-orbital direction; ignored for s-s.
/// `diff` is the vector from center_j to center_i (Bohr).
/// Returns the contracted GTO overlap for the specified pair.
fn contracted_overlap_ij(
    alphas_i: &[f64],
    coeffs_i: &[f64],
    ot_i: OrbType,
    alphas_j: &[f64],
    coeffs_j: &[f64],
    ot_j: OrbType,
    diff: [f64; 3],
    r2: f64,
) -> f64 {
    use std::f64::consts::PI;
    let mut s = 0.0_f64;

    for (&ai, &ci) in alphas_i.iter().zip(coeffs_i.iter()) {
        let ni = if ot_i == OrbType::S {
            norm_s(ai)
        } else {
            norm_p(ai)
        };
        for (&aj, &cj) in alphas_j.iter().zip(coeffs_j.iter()) {
            let nj = if ot_j == OrbType::S {
                norm_s(aj)
            } else {
                norm_p(aj)
            };
            let gamma = ai + aj;
            let k_ss = (PI / gamma).powf(1.5) * (-(ai * aj / gamma) * r2).exp();
            let pref = ci * ni * cj * nj;

            let val = match (ot_i, ot_j) {
                (OrbType::S, OrbType::S) => {
                    // s-s: Gaussian product formula
                    k_ss
                }
                (OrbType::S, op) | (op, OrbType::S) => {
                    // s-p or p-s: <s_A | p_k^B> = K_AB * (P_k - X_B_k) * (π/γ)^(3/2)
                    // P_k - X_B_k = α_A * (X_A_k - X_B_k) / γ  for (s on A, p on B)
                    // For (p on A, s on B): <p_k^A | s_B> = K_AB * (P_k - X_A_k) = α_B*(X_B_k-X_A_k)/γ
                    let axis = match op {
                        OrbType::Px => 0_usize,
                        OrbType::Py => 1,
                        OrbType::Pz => 2,
                        OrbType::S => unreachable!(),
                    };
                    // diff = center_i - center_j
                    let (alpha_s, sign) = if ot_i == OrbType::S {
                        // s on i, p on j: P_k - X_j_k = α_i * diff[axis] / γ
                        (ai, 1.0_f64)
                    } else {
                        // p on i, s on j: P_k - X_i_k = α_j * (-diff[axis]) / γ
                        (aj, -1.0_f64)
                    };
                    let proj = alpha_s * diff[axis] * sign / gamma;
                    k_ss * proj
                }
                (ot_a, ot_b) => {
                    // p-p: three cases — same axis or different axis
                    let axis_a = match ot_a {
                        OrbType::Px => 0,
                        OrbType::Py => 1,
                        OrbType::Pz => 2,
                        OrbType::S => unreachable!(),
                    };
                    let axis_b = match ot_b {
                        OrbType::Px => 0,
                        OrbType::Py => 1,
                        OrbType::Pz => 2,
                        OrbType::S => unreachable!(),
                    };

                    if axis_a == axis_b {
                        // Same axis: <pk_A | pk_B> = K * ((Pk-XA)*(Pk-XB) + 1/(2γ)) * (π/γ)^1.5
                        // Pk - XA = α_B*(XB_k - XA_k)/γ = -α_B*diff_k/γ
                        // Pk - XB = α_A*(XA_k - XB_k)/γ = α_A*diff_k/γ
                        let d_k = diff[axis_a];
                        let cross = (-aj * d_k / gamma) * (ai * d_k / gamma);
                        k_ss * (cross + 1.0 / (2.0 * gamma))
                    } else {
                        // Different axes: <pk_A | pl_B> = K * (Pk-XA)*(Pl-XB) * (π/γ)^1.5
                        // Pk - XA_k = -α_B*diff_k/γ
                        // Pl - XB_l = α_A*diff_l/γ
                        let d_k = diff[axis_a];
                        let d_l = diff[axis_b];
                        let cross = (-aj * d_k / gamma) * (ai * d_l / gamma);
                        k_ss * cross
                    }
                }
            };
            s += pref * val;
        }
    }
    s
}

/// Build the N×N normalized overlap matrix S.
///
/// Uses full angular-momentum-aware overlap formulas (s-s, s-p, p-s, p-p)
/// with explicit px/py/pz distinction from `orb_types`.
/// Off-diagonal same-atom p-p (e.g. px-py on same atom): exactly 0 by symmetry.
/// Diagonal elements are renormalized to 1 to account for STO-3G contraction conventions.
fn build_overlap(basis: &[BasisFunction], orb_types: &[OrbType]) -> Vec<Vec<f64>> {
    let n = basis.len();
    let mut s = vec![vec![0.0_f64; n]; n];

    // Compute raw overlaps
    for i in 0..n {
        let ot_i = orb_types[i];
        s[i][i] = contracted_overlap_ij(
            &basis[i].exponents,
            &basis[i].coefficients,
            ot_i,
            &basis[i].exponents,
            &basis[i].coefficients,
            ot_i,
            [0.0; 3],
            0.0,
        );

        for j in (i + 1)..n {
            let ot_j = orb_types[j];
            let diff = [
                basis[i].center[0] - basis[j].center[0],
                basis[i].center[1] - basis[j].center[1],
                basis[i].center[2] - basis[j].center[2],
            ];
            let r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
            let sij = contracted_overlap_ij(
                &basis[i].exponents,
                &basis[i].coefficients,
                ot_i,
                &basis[j].exponents,
                &basis[j].coefficients,
                ot_j,
                diff,
                r2,
            );
            s[i][j] = sij;
            s[j][i] = sij;
        }
    }

    // Renormalize: S'_ij = S_ij / sqrt(S_ii * S_jj) so that S'_ii = 1.
    // This handles the STO-3G contraction normalization convention.
    let norms: Vec<f64> = (0..n)
        .map(|i| {
            let v = s[i][i];
            if v > 0.0 { v.sqrt() } else { 1.0 }
        })
        .collect();
    for i in 0..n {
        for j in 0..n {
            s[i][j] /= norms[i] * norms[j];
        }
    }

    s
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
fn build_gamma(region: &QmRegion, params: &[Option<Pm3Params>]) -> Vec<Vec<f64>> {
    let n = region.n_atoms();
    let mut gamma = vec![vec![0.0_f64; n]; n];
    for a in 0..n {
        let ga = params[a]
            .as_ref()
            .map_or(14.4 / HARTREE_TO_EV, |p| p.g_ss / HARTREE_TO_EV);
        for b in a..n {
            let gb = params[b]
                .as_ref()
                .map_or(14.4 / HARTREE_TO_EV, |p| p.g_ss / HARTREE_TO_EV);
            let r_bohr = if a == b {
                0.0
            } else {
                let dx = region.positions[a][0] - region.positions[b][0];
                let dy = region.positions[a][1] - region.positions[b][1];
                let dz = region.positions[a][2] - region.positions[b][2];
                ((dx * dx + dy * dy + dz * dz).sqrt()) * ANG_TO_BOHR
            };
            let g = if a == b {
                ga
            } else {
                klopman_ohno(r_bohr, ga, gb)
            };
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
fn build_h_core(
    region: &QmRegion,
    basis: &[BasisFunction],
    orb_types: &[OrbType],
    atom_of_basis: &[usize],
    params_per_basis: &[Pm3Params],
    overlap: &[Vec<f64>],
    gamma: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n = basis.len();
    let n_atoms = region.n_atoms();
    let mut h = vec![vec![0.0_f64; n]; n];

    // Precompute core charges Z_A = n_val
    let core_charge: Vec<f64> = region
        .atomic_numbers
        .iter()
        .filter_map(|&z| pm3_params(z).map(|p| p.n_val as f64))
        .collect();

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
                    let u_val = ot_mu.u_energy(p_mu);
                    let mut v_sum = 0.0;
                    for c in 0..n_atoms {
                        if c != a && c < core_charge.len() {
                            v_sum -= core_charge[c] * gamma[a][c];
                        }
                    }
                    h[mu][nu] = u_val + v_sum;
                }
                // Off-diagonal same atom: 0 in NDDO (orthogonal AOs on same centre)
            } else {
                // Two-centre
                let beta_mu = ot_mu.beta(p_mu) / HARTREE_TO_EV;
                let beta_nu = ot_nu.beta(p_nu) / HARTREE_TO_EV;
                h[mu][nu] = 0.5 * (beta_mu + beta_nu) * overlap[mu][nu];
            }
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
        generalized_symmetric_eigen_n(h_core, overlap).ok_or(Pm3Error::EigensolverFailed)?;

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
            generalized_symmetric_eigen_n(&fock, overlap).ok_or(Pm3Error::EigensolverFailed)?;

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

    Err(Pm3Error::ScfNotConverged)
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
fn core_core_repulsion(region: &QmRegion) -> f64 {
    let n = region.n_atoms();
    let params: Vec<Option<Pm3Params>> = region
        .atomic_numbers
        .iter()
        .map(|&z| pm3_params(z))
        .collect();

    let mut e_cc = 0.0_f64;

    for a in 0..n {
        let pa = match &params[a] {
            Some(p) => p,
            None => continue,
        };
        let za = pa.n_val as f64;

        for (b, pb_opt) in params.iter().enumerate().skip(a + 1) {
            let pb = match pb_opt {
                Some(p) => p,
                None => continue,
            };
            let zb = pb.n_val as f64;

            let dx = region.positions[a][0] - region.positions[b][0];
            let dy = region.positions[a][1] - region.positions[b][1];
            let dz = region.positions[a][2] - region.positions[b][2];
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

/// Compute the PM3 total electronic energy for the given positions.
///
/// Returns energy in Hartree on success.
fn single_point_energy(region: &QmRegion) -> Result<f64, Pm3Error> {
    if region.n_atoms() == 0 {
        return Err(Pm3Error::EmptyRegion);
    }

    // Check all elements are supported
    for &z in &region.atomic_numbers {
        pm3_params(z).ok_or(Pm3Error::UnsupportedElement(z))?;
    }

    // Total valence electrons
    let n_val_total: usize = region
        .atomic_numbers
        .iter()
        .filter_map(|&z| pm3_params(z).map(|p| p.n_val))
        .sum();

    let n_electrons = (n_val_total as i64 - region.charge as i64) as usize;
    if !n_electrons.is_multiple_of(2) || region.multiplicity != 1 {
        return Err(Pm3Error::OpenShellNotSupported);
    }
    let n_occ = n_electrons / 2;

    let (basis, orb_types, atom_of_basis, params_per_basis) = build_basis(region)?;

    // Precompute per-atom params for γ
    let atom_params: Vec<Option<Pm3Params>> = region
        .atomic_numbers
        .iter()
        .map(|&z| pm3_params(z))
        .collect();

    let overlap = build_overlap(&basis, &orb_types);
    let gamma = build_gamma(region, &atom_params);
    let h_core = build_h_core(
        region,
        &basis,
        &orb_types,
        &atom_of_basis,
        &params_per_basis,
        &overlap,
        &gamma,
    );

    let (density, _evals, _evecs) =
        scf_loop(&h_core, &overlap, n_occ, &atom_of_basis, &gamma, 500)?;

    // Electronic energy
    let nbasis = basis.len();
    let g = build_g_matrix(&density, &overlap, &atom_of_basis, &gamma);
    let mut fock = vec![vec![0.0_f64; nbasis]; nbasis];
    for mu in 0..nbasis {
        for nu in 0..nbasis {
            fock[mu][nu] = h_core[mu][nu] + g[mu][nu];
        }
    }

    let mut e_elec = 0.0_f64;
    for mu in 0..nbasis {
        for nu in 0..nbasis {
            e_elec += density[mu][nu] * (h_core[mu][nu] + fock[mu][nu]);
        }
    }
    e_elec *= 0.5;

    let e_cc = core_core_repulsion(region);
    Ok(e_elec + e_cc)
}

// ---------------------------------------------------------------------------
// Numerical gradient
// ---------------------------------------------------------------------------

/// Compute forces on all atoms via central-difference finite differentiation.
///
/// F_ix = -(E(+δ) - E(-δ)) / (2δ)  for δ = 0.001 Å
///
/// This requires 6·N SCF calls, which is acceptable for small QM regions.
fn numerical_gradient(region: &QmRegion) -> Vec<[f64; 3]> {
    let n = region.n_atoms();
    let delta = 0.001_f64; // Å
    let mut forces = vec![[0.0_f64; 3]; n];

    for (i, f_i) in forces.iter_mut().enumerate() {
        for (k, f_ik) in f_i.iter_mut().enumerate() {
            // Forward displacement
            let mut reg_plus = clone_region(region);
            reg_plus.positions[i][k] += delta;

            // Backward displacement
            let mut reg_minus = clone_region(region);
            reg_minus.positions[i][k] -= delta;

            let e_plus = single_point_energy(&reg_plus).unwrap_or(0.0);
            let e_minus = single_point_energy(&reg_minus).unwrap_or(0.0);

            // Force = -dE/dx, in Hartree/Å → convert to kcal/mol/Å
            let grad_hartree_per_ang = (e_plus - e_minus) / (2.0 * delta);
            *f_ik = -grad_hartree_per_ang * HARTREE_TO_KCAL;
        }
    }
    forces
}

/// Minimal clone of a `QmRegion` for finite-difference displacements.
fn clone_region(region: &QmRegion) -> QmRegion {
    let mut r = QmRegion::new(
        region.atom_indices.clone(),
        region.method,
        region.charge,
        region.multiplicity,
    );
    r.positions = region.positions.clone();
    r.atomic_numbers = region.atomic_numbers.clone();
    r
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Internal PM3 calculation returning `Err` on any failure.
fn try_run_pm3(region: &mut QmRegion) -> Result<(), Pm3Error> {
    if region.n_atoms() == 0 {
        return Err(Pm3Error::EmptyRegion);
    }

    // Check element support and multiplicity first
    for &z in &region.atomic_numbers {
        pm3_params(z).ok_or(Pm3Error::UnsupportedElement(z))?;
    }
    if region.multiplicity != 1 {
        return Err(Pm3Error::OpenShellNotSupported);
    }

    let e_hartree = single_point_energy(region)?;
    region.energy = e_hartree * HARTREE_TO_KCAL;

    region.forces = numerical_gradient(region);
    Ok(())
}

/// Run a PM3 semi-empirical SCF calculation, storing energy (kcal/mol) and
/// forces (kcal/mol/Å) into `region`.  Silently falls back to zero if the
/// element is unsupported or SCF fails.
pub fn run_pm3(region: &mut QmRegion) {
    if let Err(_e) = try_run_pm3(region) {
        // Energy/forces stay at their default (0.0) — caller can detect this via
        // checking region.energy == 0.0 or via ScfStatus from a higher-level wrapper.
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qm_mm::QmMethod;

    #[test]
    fn h2_pm3_energy_negative() {
        // H-H bond at 0.74 Å (experimental)
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_pm3(&mut region);
        assert!(
            region.energy < 0.0,
            "H2 PM3 energy should be negative: {}",
            region.energy
        );
    }

    #[test]
    fn h2_forces_sum_to_zero() {
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_pm3(&mut region);
        let fx: f64 = region.forces.iter().map(|f| f[0]).sum();
        let fy: f64 = region.forces.iter().map(|f| f[1]).sum();
        let fz: f64 = region.forces.iter().map(|f| f[2]).sum();
        assert!(fx.abs() < 1.0, "force x sum not zero: {fx}");
        assert!(fy.abs() < 0.1, "force y sum not zero: {fy}");
        assert!(fz.abs() < 0.1, "force z sum not zero: {fz}");
    }

    #[test]
    fn h2o_pm3_converges() {
        // Water molecule: O at origin, two H at ~0.96 Å, ~104.5° angle
        let mut region = QmRegion::new(vec![0, 1, 2], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![8, 1, 1];
        region.set_positions(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.586],
            [0.0, -0.757, 0.586],
        ]);
        run_pm3(&mut region);
        assert!(region.energy.is_finite(), "water PM3 energy must be finite");
        assert!(region.energy < 0.0, "water PM3 energy should be negative");
    }

    #[test]
    fn unsupported_element_does_not_panic() {
        // Ne (Z=10) is not supported; run_pm3 should not panic, just leave energy=0
        let mut region = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        region.atomic_numbers = vec![10];
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        run_pm3(&mut region); // must not panic
    }
}
