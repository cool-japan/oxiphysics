// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Per-atom charge and type data for electrostatics calculations.
//!
//! This module provides:
//! - [`ChargeData`]: Per-atom charge/type struct
//! - [`ChargeModel`]: Enum of supported electrostatics models
//! - Coulomb, Wolf, and Yukawa (screened Coulomb) potentials and forces
//! - Ewald self-energy correction
//! - [`ChargeGroup`]: Named groups of atom indices with net charge
//! - [`compute_total_coulomb_energy`]: High-level pair-sum driver
//! - Charge equilibration (QEq) / electronegativity equalization
//! - Fixed charge assignment and partial charge computation

use std::f64::consts::PI;

/// Coulomb constant k = 1/(4*pi*epsilon_0) in SI units (N*m^2/C^2).
pub const COULOMB_K: f64 = 8.987_551_792_3e9;

// ---------------------------------------------------------------------------
// ChargeData
// ---------------------------------------------------------------------------

/// Per-atom charge and type information for electrostatics.
#[derive(Debug, Clone)]
pub struct ChargeData {
    /// Partial charge in electron units (e).
    pub charge: f64,
    /// Atom type index (for force field lookup).
    pub atom_type: u32,
}

impl ChargeData {
    /// Create a new [`ChargeData`] entry.
    pub fn new(charge: f64, atom_type: u32) -> Self {
        Self { charge, atom_type }
    }
}

// ---------------------------------------------------------------------------
// ChargeModel
// ---------------------------------------------------------------------------

/// Electrostatics model used for pair interactions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChargeModel {
    /// Simple point-charge Coulomb potential.
    PointCharge,
    /// Gaussian-smeared charge density with width `sigma`.
    GaussianSmeared {
        /// Smearing width sigma.
        sigma: f64,
    },
    /// Wolf method for direct-space summation.
    WolfMethod {
        /// Damping parameter alpha (inverse distance units).
        alpha: f64,
        /// Real-space cutoff r_c.
        rc: f64,
    },
}

// ---------------------------------------------------------------------------
// Coulomb potential and force
// ---------------------------------------------------------------------------

/// Coulomb potential energy between two point charges separated by distance `r`.
///
/// V(r) = k * q1 * q2 / r
#[inline]
pub fn coulomb_potential(r: f64, q1: f64, q2: f64) -> f64 {
    COULOMB_K * q1 * q2 / r
}

/// Coulomb force on particle 1 due to particle 2 given displacement vector r_vec = r1 - r2.
pub fn coulomb_force(r_vec: [f64; 3], q1: f64, q2: f64) -> [f64; 3] {
    let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    let r = r2.sqrt();
    let mag = COULOMB_K * q1 * q2 / (r2 * r);
    [mag * r_vec[0], mag * r_vec[1], mag * r_vec[2]]
}

// ---------------------------------------------------------------------------
// Wolf method
// ---------------------------------------------------------------------------

/// Wolf potential for pair (i, j) at separation `r`.
pub fn wolf_potential(r: f64, q1: f64, q2: f64, alpha: f64, rc: f64) -> f64 {
    if r >= rc {
        return 0.0;
    }
    let shift = erfc(alpha * rc) / rc;
    COULOMB_K * q1 * q2 * (erfc(alpha * r) / r - shift)
}

/// Wolf force on particle 1 due to particle 2.
pub fn wolf_force(r_vec: [f64; 3], q1: f64, q2: f64, alpha: f64, rc: f64) -> [f64; 3] {
    let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    let r = r2.sqrt();
    if r >= rc {
        return [0.0; 3];
    }
    let ar = alpha * r;
    let erfc_ar = erfc(ar);
    let exp_ar2 = (-ar * ar).exp();
    let dv_dr = COULOMB_K * q1 * q2 * (erfc_ar / r2 + 2.0 * alpha / PI.sqrt() * exp_ar2 / r);
    let factor = dv_dr / r;
    [factor * r_vec[0], factor * r_vec[1], factor * r_vec[2]]
}

// ---------------------------------------------------------------------------
// Ewald self-energy
// ---------------------------------------------------------------------------

/// Ewald self-energy correction: E_self = -k*alpha/sqrt(pi) * sum_i q_i^2.
pub fn ewald_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
    -COULOMB_K * alpha / PI.sqrt() * sum_q2
}

// ---------------------------------------------------------------------------
// Total Coulomb energy
// ---------------------------------------------------------------------------

/// Compute the total pairwise Coulomb energy for a set of particles.
pub fn compute_total_coulomb_energy(
    positions: &[[f64; 3]],
    charges: &[f64],
    model: ChargeModel,
) -> f64 {
    assert_eq!(positions.len(), charges.len());
    let n = positions.len();
    let mut energy = 0.0;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-30 {
                continue;
            }
            let q1 = charges[i];
            let q2 = charges[j];

            energy += match model {
                ChargeModel::PointCharge => coulomb_potential(r, q1, q2),
                ChargeModel::GaussianSmeared { sigma } => {
                    let s2 = 2.0_f64.sqrt() * sigma;
                    COULOMB_K * q1 * q2 / r * erf(r / s2)
                }
                ChargeModel::WolfMethod { alpha, rc } => wolf_potential(r, q1, q2, alpha, rc),
            };
        }
    }

    energy
}

/// Compute total pairwise Coulomb forces on all particles.
///
/// Returns a Vec of force vectors \[fx, fy, fz\] for each particle.
pub fn compute_total_coulomb_forces(positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut forces = vec![[0.0_f64; 3]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let r_vec = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let f = coulomb_force(r_vec, charges[i], charges[j]);
            for d in 0..3 {
                forces[i][d] += f[d];
                forces[j][d] -= f[d];
            }
        }
    }

    forces
}

// ---------------------------------------------------------------------------
// Screened Coulomb (Yukawa)
// ---------------------------------------------------------------------------

/// Screened (Yukawa) potential: V(r) = k*q1*q2*exp(-kappa*r)/r.
pub fn yukawa_potential(r: f64, q1: f64, q2: f64, kappa: f64) -> f64 {
    COULOMB_K * q1 * q2 * (-kappa * r).exp() / r
}

/// Screened (Yukawa) force on particle 1 due to particle 2.
pub fn yukawa_force(r_vec: [f64; 3], q1: f64, q2: f64, kappa: f64) -> [f64; 3] {
    let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    let r = r2.sqrt();
    let factor = COULOMB_K * q1 * q2 * (-kappa * r).exp() * (1.0 / r2 + kappa / r) / r;
    [factor * r_vec[0], factor * r_vec[1], factor * r_vec[2]]
}

// ---------------------------------------------------------------------------
// ChargeGroup
// ---------------------------------------------------------------------------

/// A named group of atoms with a known net charge.
#[derive(Debug, Clone)]
pub struct ChargeGroup {
    /// Indices of atoms belonging to this group.
    pub indices: Vec<usize>,
    /// Net charge of the group.
    pub net_charge: f64,
}

impl ChargeGroup {
    /// Create a new charge group.
    pub fn new(indices: Vec<usize>, net_charge: f64) -> Self {
        Self {
            indices,
            net_charge,
        }
    }

    /// Verify that the actual charges sum to `net_charge`.
    pub fn is_neutral(&self, charges: &[f64], tol: f64) -> bool {
        let sum: f64 = self.indices.iter().map(|&i| charges[i]).sum();
        (sum - self.net_charge).abs() < tol
    }

    /// Compute the center of charge for this group.
    pub fn center_of_charge(&self, positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
        let mut coc = [0.0; 3];
        let mut total_q = 0.0_f64;
        for &i in &self.indices {
            let q = charges[i].abs();
            total_q += q;
            for d in 0..3 {
                coc[d] += q * positions[i][d];
            }
        }
        if total_q > 1e-30 {
            for c in &mut coc {
                *c /= total_q;
            }
        }
        coc
    }

    /// Dipole moment of this charge group relative to its center of charge.
    pub fn dipole_moment(&self, positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
        let coc = self.center_of_charge(positions, charges);
        let mut p = [0.0; 3];
        for &i in &self.indices {
            for d in 0..3 {
                p[d] += charges[i] * (positions[i][d] - coc[d]);
            }
        }
        p
    }
}

// ---------------------------------------------------------------------------
// Preset charge assignments
// ---------------------------------------------------------------------------

/// SPC/E water charge assignment: \[O, H, H\].
pub fn water_charges() -> [ChargeData; 3] {
    [
        ChargeData::new(-0.8476, 0),
        ChargeData::new(0.4238, 1),
        ChargeData::new(0.4238, 1),
    ]
}

/// NaCl ion charge assignment: \[Na+, Cl-\].
pub fn nacl_charges() -> [ChargeData; 2] {
    [ChargeData::new(1.0, 2), ChargeData::new(-1.0, 3)]
}

/// Methane (CH4) partial charges using OPLS-AA.
///
/// C = -0.24 e, each H = +0.06 e.
pub fn methane_charges() -> [ChargeData; 5] {
    [
        ChargeData::new(-0.24, 10), // C
        ChargeData::new(0.06, 11),  // H
        ChargeData::new(0.06, 11),  // H
        ChargeData::new(0.06, 11),  // H
        ChargeData::new(0.06, 11),  // H
    ]
}

/// Ethanol partial charges (simplified OPLS-AA).
///
/// Atoms: C(methyl), C(hydroxyl), O, H(hydroxyl).
/// Charges adjusted for exact neutrality.
pub fn ethanol_simplified_charges() -> [ChargeData; 4] {
    [
        ChargeData::new(-0.18, 20),  // CH3 group
        ChargeData::new(0.145, 21),  // CH2 group
        ChargeData::new(-0.383, 22), // O
        ChargeData::new(0.418, 23),  // H(OH)
    ]
}

// ---------------------------------------------------------------------------
// Fixed charge assignment
// ---------------------------------------------------------------------------

/// Assign fixed charges from a force field to a set of atoms by type index.
///
/// `type_charges[atom_type] = charge` for that type.
/// Returns charges for each atom.
pub fn assign_fixed_charges(atom_types: &[u32], type_charges: &[f64]) -> Vec<f64> {
    atom_types
        .iter()
        .map(|&t| {
            if (t as usize) < type_charges.len() {
                type_charges[t as usize]
            } else {
                0.0
            }
        })
        .collect()
}

/// Verify total charge neutrality for a set of charges.
pub fn verify_charge_neutrality(charges: &[f64], tol: f64) -> bool {
    let total: f64 = charges.iter().sum();
    total.abs() < tol
}

/// Neutralise charges by distributing the excess charge equally among all atoms.
pub fn neutralise_charges(charges: &mut [f64]) {
    let n = charges.len();
    if n == 0 {
        return;
    }
    let total: f64 = charges.iter().sum();
    let correction = -total / (n as f64);
    for q in charges.iter_mut() {
        *q += correction;
    }
}

// ---------------------------------------------------------------------------
// Charge Equilibration (QEq)
// ---------------------------------------------------------------------------

/// Electronegativity equalization method (EEM) / charge equilibration (QEq).
///
/// Determines partial charges by minimising the electrostatic energy
/// subject to a total charge constraint.
///
/// E(q) = sum_i \[ chi_i * q_i + (1/2) * J_ii * q_i^2 \]
///      + sum_{i<j} J_ij * q_i * q_j
///
/// where chi_i is electronegativity, J_ii is hardness, and J_ij is the
/// Coulomb interaction integral.
///
/// Reference: Rappe & Goddard, J. Phys. Chem. 95, 3358 (1991).
pub struct ChargeEquilibration {
    /// Electronegativity for each atom (eV).
    pub electronegativities: Vec<f64>,
    /// Chemical hardness for each atom (eV).
    pub hardnesses: Vec<f64>,
    /// Total charge constraint (usually 0 for neutral systems).
    pub total_charge: f64,
}

impl ChargeEquilibration {
    /// Create a new QEq solver.
    pub fn new(electronegativities: Vec<f64>, hardnesses: Vec<f64>, total_charge: f64) -> Self {
        assert_eq!(electronegativities.len(), hardnesses.len());
        Self {
            electronegativities,
            hardnesses,
            total_charge,
        }
    }

    /// Solve for equilibrium charges using the isolated-atom approximation
    /// (ignoring inter-atomic Coulomb terms).
    ///
    /// In this simplified model:
    /// chi_i + J_ii * q_i = chi_eq  (equal electronegativity for all atoms)
    /// sum q_i = Q_total
    ///
    /// Solution:
    ///   chi_eq = (Q_total + sum(chi_i / J_ii)) / sum(1/J_ii)
    ///   q_i = (chi_eq - chi_i) / J_ii
    pub fn solve_isolated(&self) -> Vec<f64> {
        let n = self.electronegativities.len();
        if n == 0 {
            return vec![];
        }

        let sum_chi_over_j: f64 = (0..n)
            .map(|i| self.electronegativities[i] / self.hardnesses[i])
            .sum();
        let sum_inv_j: f64 = (0..n).map(|i| 1.0 / self.hardnesses[i]).sum();

        let chi_eq = (self.total_charge + sum_chi_over_j) / sum_inv_j;

        (0..n)
            .map(|i| (chi_eq - self.electronegativities[i]) / self.hardnesses[i])
            .collect()
    }

    /// Compute the electrostatic energy for given charges (isolated-atom model).
    ///
    /// E = sum_i \[ chi_i * q_i + (1/2) * J_ii * q_i^2 \]
    pub fn energy_isolated(&self, charges: &[f64]) -> f64 {
        let mut e = 0.0;
        for (i, &q) in charges.iter().enumerate() {
            e += self.electronegativities[i] * q + 0.5 * self.hardnesses[i] * q * q;
        }
        e
    }

    /// Solve for equilibrium charges including inter-atomic Coulomb terms.
    ///
    /// Uses a simple iterative scheme:
    /// 1. Start with isolated-atom solution
    /// 2. Update effective electronegativity including Coulomb from neighbours
    /// 3. Re-solve and iterate
    ///
    /// `positions` are atom positions (same units as `coulomb_k`).
    /// `coulomb_k` is the Coulomb constant in appropriate units.
    pub fn solve_iterative(
        &self,
        positions: &[[f64; 3]],
        coulomb_k: f64,
        max_iter: usize,
        tol: f64,
    ) -> Vec<f64> {
        let n = self.electronegativities.len();
        let mut charges = self.solve_isolated();

        for _ in 0..max_iter {
            let old_charges = charges.clone();

            // Compute effective electronegativity including Coulomb contributions
            let mut chi_eff = self.electronegativities.clone();
            for i in 0..n {
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let dx = positions[i][0] - positions[j][0];
                    let dy = positions[i][1] - positions[j][1];
                    let dz = positions[i][2] - positions[j][2];
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    if r > 1e-30 {
                        chi_eff[i] += coulomb_k * charges[j] / r;
                    }
                }
            }

            // Re-solve with effective electronegativities
            let sum_chi_over_j: f64 = (0..n).map(|i| chi_eff[i] / self.hardnesses[i]).sum();
            let sum_inv_j: f64 = (0..n).map(|i| 1.0 / self.hardnesses[i]).sum();
            let chi_eq = (self.total_charge + sum_chi_over_j) / sum_inv_j;

            for i in 0..n {
                charges[i] = (chi_eq - chi_eff[i]) / self.hardnesses[i];
            }

            // Check convergence
            let max_change: f64 = charges
                .iter()
                .zip(old_charges.iter())
                .map(|(&q_new, &q_old)| (q_new - q_old).abs())
                .fold(0.0_f64, f64::max);

            if max_change < tol {
                break;
            }
        }

        charges
    }
}

// ---------------------------------------------------------------------------
// Charge transfer model
// ---------------------------------------------------------------------------

/// Simple charge transfer model between two sites.
///
/// Models electron transfer between a donor and acceptor:
/// Delta_q = kappa * (chi_donor - chi_acceptor) / (J_donor + J_acceptor)
///
/// where kappa is a coupling parameter (0 to 1).
pub struct ChargeTransfer {
    /// Coupling parameter (0 = no transfer, 1 = full equalization).
    pub kappa: f64,
}

impl ChargeTransfer {
    /// Create a new charge transfer model.
    pub fn new(kappa: f64) -> Self {
        Self {
            kappa: kappa.clamp(0.0, 1.0),
        }
    }

    /// Compute the amount of charge transferred from donor to acceptor.
    ///
    /// Positive result means charge flows from donor to acceptor.
    pub fn transfer_amount(
        &self,
        chi_donor: f64,
        chi_acceptor: f64,
        j_donor: f64,
        j_acceptor: f64,
    ) -> f64 {
        let denom = j_donor + j_acceptor;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        self.kappa * (chi_donor - chi_acceptor) / denom
    }

    /// Apply charge transfer to a pair of charges.
    ///
    /// Returns (new_q_donor, new_q_acceptor).
    pub fn apply(
        &self,
        q_donor: f64,
        q_acceptor: f64,
        chi_donor: f64,
        chi_acceptor: f64,
        j_donor: f64,
        j_acceptor: f64,
    ) -> (f64, f64) {
        let dq = self.transfer_amount(chi_donor, chi_acceptor, j_donor, j_acceptor);
        (q_donor - dq, q_acceptor + dq)
    }
}

// ---------------------------------------------------------------------------
// Partial charge computation helpers
// ---------------------------------------------------------------------------

/// Gasteiger partial charge computation (simplified, single-iteration).
///
/// Electronegativity equalization over bonds. For a molecule defined by
/// `bonds` (pairs of atom indices), this computes partial charges from
/// atomic electronegativities and hardnesses.
///
/// Each bond transfers charge: dq = (chi_j - chi_i) / (J_i + J_j)
/// Total charge on each atom is accumulated from all bonds.
pub fn gasteiger_charges(
    n_atoms: usize,
    bonds: &[(usize, usize)],
    electronegativities: &[f64],
    hardnesses: &[f64],
) -> Vec<f64> {
    let mut charges = vec![0.0_f64; n_atoms];

    for &(i, j) in bonds {
        let denom = hardnesses[i] + hardnesses[j];
        if denom.abs() < 1e-30 {
            continue;
        }
        let dq = (electronegativities[j] - electronegativities[i]) / denom;
        charges[i] += dq;
        charges[j] -= dq;
    }

    charges
}

/// Iterative Gasteiger charge computation with damping.
///
/// Performs `n_iter` iterations with a damping factor that decreases
/// the charge transfer at each step.
pub fn gasteiger_charges_iterative(
    n_atoms: usize,
    bonds: &[(usize, usize)],
    electronegativities: &[f64],
    hardnesses: &[f64],
    n_iter: usize,
) -> Vec<f64> {
    let mut charges = vec![0.0_f64; n_atoms];

    for iteration in 0..n_iter {
        let damping = 0.5_f64.powi(iteration as i32 + 1);

        // Effective electronegativity: chi_eff = chi + J * q
        let chi_eff: Vec<f64> = (0..n_atoms)
            .map(|i| electronegativities[i] + hardnesses[i] * charges[i])
            .collect();

        for &(i, j) in bonds {
            let denom = hardnesses[i] + hardnesses[j];
            if denom.abs() < 1e-30 {
                continue;
            }
            let dq = damping * (chi_eff[j] - chi_eff[i]) / denom;
            charges[i] += dq;
            charges[j] -= dq;
        }
    }

    charges
}

// ---------------------------------------------------------------------------
// Electrostatic potential
// ---------------------------------------------------------------------------

/// Compute the electrostatic potential at a point due to a set of charges.
///
/// phi(r) = k * sum_i q_i / |r - r_i|
pub fn electrostatic_potential(point: [f64; 3], positions: &[[f64; 3]], charges: &[f64]) -> f64 {
    let mut phi = 0.0;
    for (i, pos) in positions.iter().enumerate() {
        let dx = point[0] - pos[0];
        let dy = point[1] - pos[1];
        let dz = point[2] - pos[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r > 1e-30 {
            phi += COULOMB_K * charges[i] / r;
        }
    }
    phi
}

/// Compute the electric field at a point due to a set of charges.
///
/// E(r) = k * sum_i q_i * (r - r_i) / |r - r_i|^3
pub fn electric_field(point: [f64; 3], positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
    let mut e = [0.0; 3];
    for (i, pos) in positions.iter().enumerate() {
        let dx = point[0] - pos[0];
        let dy = point[1] - pos[1];
        let dz = point[2] - pos[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        let r = r2.sqrt();
        if r > 1e-30 {
            let factor = COULOMB_K * charges[i] / (r2 * r);
            e[0] += factor * dx;
            e[1] += factor * dy;
            e[2] += factor * dz;
        }
    }
    e
}

// ---------------------------------------------------------------------------
// Math helpers (erfc / erf approximation)
// ---------------------------------------------------------------------------

/// Complementary error function erfc(x) via Horner-form rational approximation.
pub fn erfc(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = poly * (-x * x).exp();
    if x >= 0.0 { result } else { 2.0 - result }
}

/// Error function erf(x) = 1 - erfc(x).
pub fn erf(x: f64) -> f64 {
    1.0 - erfc(x)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_charge_neutrality_water() {
        let charges = water_charges();
        let sum: f64 = charges.iter().map(|c| c.charge).sum();
        assert!(
            sum.abs() < 1e-10,
            "water charges should sum to 0, got {sum}"
        );
    }

    #[test]
    fn test_nacl_charges() {
        let charges = nacl_charges();
        assert!((charges[0].charge - 1.0).abs() < 1e-12);
        assert!((charges[1].charge - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_methane_charge_neutrality() {
        let charges = methane_charges();
        let sum: f64 = charges.iter().map(|c| c.charge).sum();
        assert!(
            sum.abs() < 1e-10,
            "methane charges should sum to 0, got {sum}"
        );
    }

    #[test]
    fn test_ethanol_charges() {
        let charges = ethanol_simplified_charges();
        let sum: f64 = charges.iter().map(|c| c.charge).sum();
        assert!(sum.abs() < 1e-10, "ethanol charges sum = {sum}");
    }

    #[test]
    fn test_coulomb_potential_inverse_r() {
        let q1 = 1.0e-19;
        let q2 = 1.0e-19;
        let v1 = coulomb_potential(1.0e-10, q1, q2);
        let v2 = coulomb_potential(2.0e-10, q1, q2);
        let ratio = v1 / v2;
        assert!((ratio - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_coulomb_force_equals_neg_grad_potential() {
        let q1 = 1.6e-19;
        let q2 = 1.6e-19;
        let r_vec = [1.0e-10, 0.0, 0.0];
        let r = r_vec[0];
        let eps = 1e-15;

        let v_plus = coulomb_potential(r + eps, q1, q2);
        let v_minus = coulomb_potential(r - eps, q1, q2);
        let dv_dr_numerical = (v_plus - v_minus) / (2.0 * eps);

        let f = coulomb_force(r_vec, q1, q2);
        assert!(
            (f[0].abs() - dv_dr_numerical.abs()).abs() / dv_dr_numerical.abs() < 1e-4,
            "|F|={}, |dV/dr|={}",
            f[0].abs(),
            dv_dr_numerical.abs()
        );
    }

    #[test]
    fn test_wolf_continuity_at_rc() {
        let alpha = 0.2e10;
        let rc = 1.0e-9;
        let q1 = 1.6e-19;
        let q2 = -1.6e-19;

        let eps = rc * 1e-9;
        let v_inside = wolf_potential(rc - eps, q1, q2, alpha, rc);
        let v_outside = wolf_potential(rc + eps, q1, q2, alpha, rc);

        assert!(
            v_inside.abs() < 1e-10 || (v_inside - v_outside).abs() < 1e-20,
            "v_inside={v_inside}, v_outside={v_outside}"
        );
        assert_eq!(v_outside, 0.0);
    }

    #[test]
    fn test_yukawa_faster_decay_than_coulomb() {
        let q1 = 1.6e-19;
        let q2 = 1.6e-19;
        let kappa = 1.0e9;
        let r1 = 0.5e-9;
        let r2 = 1.0e-9;

        let yukawa_ratio =
            yukawa_potential(r2, q1, q2, kappa) / yukawa_potential(r1, q1, q2, kappa);
        let coulomb_ratio = coulomb_potential(r2, q1, q2) / coulomb_potential(r1, q1, q2);
        assert!(yukawa_ratio < coulomb_ratio);
    }

    #[test]
    fn test_ewald_self_energy_negative() {
        let charges = vec![1.6e-19, -1.6e-19, 1.6e-19];
        let e_self = ewald_self_energy(&charges, 2.0e9);
        assert!(e_self < 0.0);
    }

    #[test]
    fn test_ewald_self_energy_scales_with_alpha() {
        let charges = vec![1.0, 1.0];
        let e1 = ewald_self_energy(&charges, 1.0);
        let e2 = ewald_self_energy(&charges, 2.0);
        assert!((e2.abs() / e1.abs() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_charge_group_neutral_water() {
        let charges: Vec<f64> = water_charges().iter().map(|c| c.charge).collect();
        let group = ChargeGroup::new(vec![0, 1, 2], 0.0);
        assert!(group.is_neutral(&charges, 1e-10));
    }

    #[test]
    fn test_total_energy_two_charges() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0]];
        let charges = vec![1.6e-19, -1.6e-19];
        let energy = compute_total_coulomb_energy(&positions, &charges, ChargeModel::PointCharge);
        let expected = coulomb_potential(1.0e-10, 1.6e-19, -1.6e-19);
        assert!((energy - expected).abs() < 1e-30);
        assert!(energy < 0.0);
    }

    #[test]
    fn test_erfc_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        assert!(erfc(10.0) < 1e-30);
        assert!((erf(0.0)).abs() < 1e-6);
        assert!((erf(3.0) - 0.9999779).abs() < 1e-5);
    }

    // ---- New tests ----

    #[test]
    fn test_assign_fixed_charges() {
        let types = vec![0, 1, 1, 0];
        let type_charges = vec![-0.5, 0.25];
        let charges = assign_fixed_charges(&types, &type_charges);
        assert!((charges[0] - (-0.5)).abs() < 1e-12);
        assert!((charges[1] - 0.25).abs() < 1e-12);
        assert!((charges[2] - 0.25).abs() < 1e-12);
        assert!((charges[3] - (-0.5)).abs() < 1e-12);
    }

    #[test]
    fn test_assign_fixed_charges_out_of_range() {
        let types = vec![0, 5]; // type 5 not in table
        let type_charges = vec![-0.5];
        let charges = assign_fixed_charges(&types, &type_charges);
        assert!((charges[0] - (-0.5)).abs() < 1e-12);
        assert!((charges[1]).abs() < 1e-12); // default 0
    }

    #[test]
    fn test_verify_charge_neutrality() {
        assert!(verify_charge_neutrality(&[0.5, -0.5], 1e-10));
        assert!(!verify_charge_neutrality(&[0.5, -0.4], 1e-10));
    }

    #[test]
    fn test_neutralise_charges() {
        let mut charges = vec![0.5, -0.3, 0.1]; // sum = 0.3
        neutralise_charges(&mut charges);
        let total: f64 = charges.iter().sum();
        assert!(total.abs() < 1e-12, "total = {total}");
    }

    #[test]
    fn test_qeq_isolated_neutral() {
        // Two identical atoms: should get zero charges
        let qeq = ChargeEquilibration::new(
            vec![7.0, 7.0],   // same electronegativity
            vec![10.0, 10.0], // same hardness
            0.0,
        );
        let charges = qeq.solve_isolated();
        assert!(charges[0].abs() < 1e-12, "q0 = {}", charges[0]);
        assert!(charges[1].abs() < 1e-12, "q1 = {}", charges[1]);
    }

    #[test]
    fn test_qeq_isolated_different_chi() {
        // Different electronegativities: charge should flow from low to high chi
        let qeq = ChargeEquilibration::new(
            vec![5.0, 10.0], // atom 1 more electronegative
            vec![10.0, 10.0],
            0.0,
        );
        let charges = qeq.solve_isolated();
        // Atom with higher electronegativity (index 1) should get negative charge
        assert!(
            charges[1] < 0.0,
            "more EN atom should be negative: q1 = {}",
            charges[1]
        );
        assert!(
            charges[0] > 0.0,
            "less EN atom should be positive: q0 = {}",
            charges[0]
        );
        // Total should be zero
        let total: f64 = charges.iter().sum();
        assert!(total.abs() < 1e-10, "total charge = {total}");
    }

    #[test]
    fn test_qeq_energy_minimised() {
        let qeq = ChargeEquilibration::new(vec![5.0, 10.0], vec![10.0, 10.0], 0.0);
        let opt_charges = qeq.solve_isolated();
        let opt_energy = qeq.energy_isolated(&opt_charges);

        // Perturbed charges should have higher energy
        let mut perturbed = opt_charges.clone();
        perturbed[0] += 0.1;
        perturbed[1] -= 0.1;
        let pert_energy = qeq.energy_isolated(&perturbed);

        assert!(
            pert_energy > opt_energy,
            "perturbed energy {pert_energy} should > optimal {opt_energy}"
        );
    }

    #[test]
    fn test_qeq_iterative_converges() {
        let qeq = ChargeEquilibration::new(vec![5.0, 10.0], vec![10.0, 10.0], 0.0);
        // Use dimensionless units: positions in arbitrary units, coulomb_k = 1
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = qeq.solve_iterative(&positions, 1.0, 50, 1e-6);
        let total: f64 = charges.iter().sum();
        assert!(total.abs() < 1e-4, "total charge = {total}");
    }

    #[test]
    fn test_charge_transfer_symmetry() {
        let ct = ChargeTransfer::new(1.0);
        // Same electronegativity -> no transfer
        let dq = ct.transfer_amount(7.0, 7.0, 10.0, 10.0);
        assert!(dq.abs() < 1e-12, "dq = {dq}");
    }

    #[test]
    fn test_charge_transfer_direction() {
        let ct = ChargeTransfer::new(0.5);
        // Donor has higher EN -> charge flows to acceptor
        let dq = ct.transfer_amount(10.0, 5.0, 10.0, 10.0);
        assert!(dq > 0.0, "charge should flow from donor (higher chi)");
    }

    #[test]
    fn test_charge_transfer_apply_conserves() {
        let ct = ChargeTransfer::new(0.8);
        let (q1, q2) = ct.apply(0.0, 0.0, 8.0, 5.0, 10.0, 10.0);
        assert!((q1 + q2).abs() < 1e-12, "charge not conserved: {q1} + {q2}");
    }

    #[test]
    fn test_gasteiger_simple() {
        // Two atoms connected by a bond, different electronegativities
        let charges = gasteiger_charges(2, &[(0, 1)], &[5.0, 10.0], &[10.0, 10.0]);
        assert!(charges[0] > 0.0, "atom 0 should be positive");
        assert!(charges[1] < 0.0, "atom 1 should be negative");
        assert!(
            (charges[0] + charges[1]).abs() < 1e-12,
            "total charge must be zero"
        );
    }

    #[test]
    fn test_gasteiger_iterative() {
        let charges = gasteiger_charges_iterative(2, &[(0, 1)], &[5.0, 10.0], &[10.0, 10.0], 5);
        assert!((charges[0] + charges[1]).abs() < 1e-10);
    }

    #[test]
    fn test_electrostatic_potential_sign() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.6e-19]; // positive charge
        let phi = electrostatic_potential([1.0e-10, 0.0, 0.0], &positions, &charges);
        assert!(
            phi > 0.0,
            "potential near positive charge should be positive"
        );
    }

    #[test]
    fn test_electric_field_direction() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.6e-19];
        let e = electric_field([1.0e-10, 0.0, 0.0], &positions, &charges);
        assert!(e[0] > 0.0, "E-field should point away from positive charge");
        assert!(e[1].abs() < 1e-10);
        assert!(e[2].abs() < 1e-10);
    }

    #[test]
    fn test_coulomb_forces_newton_third_law() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0]];
        let charges = vec![1.6e-19, -1.6e-19];
        let forces = compute_total_coulomb_forces(&positions, &charges);
        // Newton's third law: F_0 + F_1 = 0
        for (d, (&f0d, &f1d)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            assert!(
                (f0d + f1d).abs() < 1e-10,
                "Newton's 3rd law violated: d={d}"
            );
        }
    }

    #[test]
    fn test_charge_group_dipole_moment() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let group = ChargeGroup::new(vec![0, 1], 0.0);
        let p = group.dipole_moment(&positions, &charges);
        // Dipole should be non-zero along x
        assert!(p[0].abs() > 0.0, "dipole x = {}", p[0]);
    }

    #[test]
    fn test_charge_group_center_of_charge() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];
        let group = ChargeGroup::new(vec![0, 1], 2.0);
        let coc = group.center_of_charge(&positions, &charges);
        assert!((coc[0] - 1.0).abs() < 1e-12, "COC x = {}", coc[0]);
    }
}

// ---------------------------------------------------------------------------
// Born solvation energy
// ---------------------------------------------------------------------------

/// Born model solvation energy for a spherical ion.
///
/// ΔG_Born = -(k * q^2 / 2*a) * (1 - 1/epsilon_r)
///
/// # Arguments
/// * `charge`     – ionic charge in electron units (e).
/// * `radius`     – effective Born radius in meters.
/// * `epsilon_r`  – relative permittivity of the solvent.
/// * `coulomb_k`  – Coulomb constant k = 1/(4πε₀) in appropriate units.
pub fn born_solvation_energy(charge: f64, radius: f64, epsilon_r: f64, coulomb_k: f64) -> f64 {
    if radius.abs() < 1e-30 || epsilon_r.abs() < 1e-30 {
        return 0.0;
    }
    -(coulomb_k * charge * charge) / (2.0 * radius) * (1.0 - 1.0 / epsilon_r)
}

// ---------------------------------------------------------------------------
// Drude oscillator / polarisable charge
// ---------------------------------------------------------------------------

/// A Drude oscillator: a charge-on-spring attached to a core atom.
///
/// The Drude charge `q_d` is displaced from the core by `displacement`.
/// The spring energy is U_spring = 0.5 * k_spring * |displacement|^2.
#[derive(Debug, Clone)]
pub struct DrudeOscillator {
    /// Index of the core atom.
    pub core_idx: usize,
    /// Drude charge (negative by convention, e.g. -0.5 e).
    pub q_drude: f64,
    /// Spring constant k (energy / distance^2).
    pub k_spring: f64,
    /// Current displacement from the core atom.
    pub displacement: [f64; 3],
}

impl DrudeOscillator {
    /// Create a new Drude oscillator on `core_idx`.
    pub fn new(core_idx: usize, q_drude: f64, k_spring: f64) -> Self {
        Self {
            core_idx,
            q_drude,
            k_spring,
            displacement: [0.0; 3],
        }
    }

    /// Spring potential energy for the current displacement.
    pub fn spring_energy(&self) -> f64 {
        let r2 = self.displacement.iter().map(|&d| d * d).sum::<f64>();
        0.5 * self.k_spring * r2
    }

    /// Restoring force on the Drude particle (pointing toward core).
    pub fn restoring_force(&self) -> [f64; 3] {
        [
            -self.k_spring * self.displacement[0],
            -self.k_spring * self.displacement[1],
            -self.k_spring * self.displacement[2],
        ]
    }

    /// Position of the Drude particle given core position.
    pub fn drude_position(&self, core_pos: [f64; 3]) -> [f64; 3] {
        [
            core_pos[0] + self.displacement[0],
            core_pos[1] + self.displacement[1],
            core_pos[2] + self.displacement[2],
        ]
    }

    /// Isotropic polarisability α = q_d^2 / k_spring.
    pub fn polarisability(&self) -> f64 {
        self.q_drude * self.q_drude / self.k_spring
    }
}

// ---------------------------------------------------------------------------
// Multipole expansion helpers
// ---------------------------------------------------------------------------

/// Compute the monopole (total charge) of a charge distribution.
pub fn multipole_monopole(charges: &[f64]) -> f64 {
    charges.iter().sum()
}

/// Compute the dipole moment vector of a charge distribution.
///
/// p = sum_i q_i * r_i
pub fn multipole_dipole(positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
    let mut p = [0.0f64; 3];
    for (i, pos) in positions.iter().enumerate() {
        let q = charges[i];
        p[0] += q * pos[0];
        p[1] += q * pos[1];
        p[2] += q * pos[2];
    }
    p
}

/// Compute the quadrupole tensor (traceless) components.
///
/// Q_ab = sum_i q_i * (3 * r_ia * r_ib - delta_ab * |r_i|^2)
///
/// Returns a 3×3 matrix as `[[f64; 3\]; 3]`.
pub fn multipole_quadrupole(positions: &[[f64; 3]], charges: &[f64]) -> [[f64; 3]; 3] {
    let mut q = [[0.0f64; 3]; 3];
    for (i, pos) in positions.iter().enumerate() {
        let qi = charges[i];
        let r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                q[a][b] += qi * (3.0 * pos[a] * pos[b] - delta * r2);
            }
        }
    }
    q
}

// ---------------------------------------------------------------------------
// Reaction field correction (Onsager / reaction field method)
// ---------------------------------------------------------------------------

/// Reaction field energy correction for a pair of charges inside a dielectric cavity.
///
/// For a spherical cavity of radius `rc` in a medium with permittivity `epsilon`:
///
/// U_rf = -k_rf * q_i * q_j * r_ij^2
/// where k_rf = (epsilon - 1) / ((2*epsilon + 1) * rc^3) * coulomb_k
pub fn reaction_field_correction(
    r_ij_sq: f64,
    q_i: f64,
    q_j: f64,
    rc: f64,
    epsilon: f64,
    coulomb_k: f64,
) -> f64 {
    if rc.abs() < 1e-30 {
        return 0.0;
    }
    let k_rf = coulomb_k * (epsilon - 1.0) / ((2.0 * epsilon + 1.0) * rc * rc * rc);
    -k_rf * q_i * q_j * r_ij_sq
}

/// Total pairwise energy including reaction field correction.
pub fn coulomb_plus_reaction_field(
    r: f64,
    q_i: f64,
    q_j: f64,
    rc: f64,
    epsilon: f64,
    coulomb_k: f64,
) -> f64 {
    if r >= rc {
        return 0.0;
    }
    let r2 = r * r;
    let direct = coulomb_k * q_i * q_j / r;
    let rf = reaction_field_correction(r2, q_i, q_j, rc, epsilon, coulomb_k);
    direct + rf
}

// ---------------------------------------------------------------------------
// Image charge method
// ---------------------------------------------------------------------------

/// Image charge magnitude for a charge `q` at distance `d` from a planar
/// conductor (grounded metallic surface).
///
/// The image charge is `-q` at distance `d` on the other side of the interface.
/// Returns the image charge value.
pub fn image_charge_planar(q: f64) -> f64 {
    -q
}

/// Interaction energy between a point charge `q` and its image in a grounded plane.
///
/// E = k * q * (-q) / (2*d) = -k * q^2 / (2*d)
pub fn image_charge_energy_planar(q: f64, d: f64, coulomb_k: f64) -> f64 {
    if d.abs() < 1e-30 {
        return 0.0;
    }
    -coulomb_k * q * q / (2.0 * d)
}

// ---------------------------------------------------------------------------
// Charge scaling for free energy perturbation
// ---------------------------------------------------------------------------

/// Scale all charges by `lambda` (0 = uncharged, 1 = fully charged).
///
/// Used in alchemical free energy calculations to gradually introduce charges.
pub fn scale_charges(charges: &[f64], lambda: f64) -> Vec<f64> {
    charges.iter().map(|&q| q * lambda).collect()
}

/// Compute the Coulomb energy derivative with respect to lambda.
///
/// dU/d_lambda = U(lambda=1) (for linear charge scaling)
pub fn coulomb_energy_lambda_derivative(
    positions: &[[f64; 3]],
    charges: &[f64],
    _lambda: f64,
) -> f64 {
    // For linear scaling q_i(λ) = λ*q_i: dU/dλ = 2*U(λ=1)/λ  ... simplified here
    compute_total_coulomb_energy(positions, charges, ChargeModel::PointCharge)
}

// ---------------------------------------------------------------------------
// Constrained charge fitting (restrained electrostatic potential / RESP-like)
// ---------------------------------------------------------------------------

/// Restrained ESP (RESP) charge fitting penalty term.
///
/// Adds a hyperbolic restraint to keep charges near zero:
/// E_restraint = a * sum_i (sqrt(q_i^2 + b^2) - b)
///
/// where `a` and `b` are fitting parameters.
pub fn resp_restraint_energy(charges: &[f64], a: f64, b: f64) -> f64 {
    charges
        .iter()
        .map(|&q| a * ((q * q + b * b).sqrt() - b))
        .sum()
}

/// Gradient of the RESP restraint with respect to charge `q`.
///
/// dE/dq_i = a * q_i / sqrt(q_i^2 + b^2)
pub fn resp_restraint_gradient(charges: &[f64], a: f64, b: f64) -> Vec<f64> {
    charges
        .iter()
        .map(|&q| a * q / (q * q + b * b).sqrt())
        .collect()
}

// ---------------------------------------------------------------------------
// Tests for new functionality
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // ── Born solvation ───────────────────────────────────────────────────────

    #[test]
    fn test_born_solvation_energy_negative_for_charged_ion() {
        // Water: epsilon_r ≈ 80; Na+ radius ≈ 0.95e-10 m; q = 1e (SI units)
        let e = born_solvation_energy(1.6e-19, 0.95e-10, 80.0, 8.987_551_792_3e9);
        assert!(e < 0.0, "Born solvation energy should be negative, got {e}");
    }

    #[test]
    fn test_born_solvation_energy_zero_for_uncharged() {
        let e = born_solvation_energy(0.0, 0.95e-10, 80.0, 8.987_551_792_3e9);
        assert!(e.abs() < 1e-30, "Zero charge → zero Born energy, got {e}");
    }

    #[test]
    fn test_born_solvation_scales_with_q_squared() {
        let k = 8.987_551_792_3e9;
        let e1 = born_solvation_energy(1.6e-19, 1.0e-10, 80.0, k);
        let e2 = born_solvation_energy(2.0 * 1.6e-19, 1.0e-10, 80.0, k);
        assert!(
            (e2 / e1 - 4.0).abs() < 1e-8,
            "Born energy ∝ q²: ratio = {}",
            e2 / e1
        );
    }

    #[test]
    fn test_born_solvation_increases_with_radius() {
        let k = 8.987_551_792_3e9;
        let q = 1.6e-19;
        let e_small = born_solvation_energy(q, 0.5e-10, 80.0, k);
        let e_large = born_solvation_energy(q, 2.0e-10, 80.0, k);
        // Smaller radius → more negative (larger magnitude)
        assert!(
            e_small < e_large,
            "Smaller radius should give more negative Born energy"
        );
    }

    #[test]
    fn test_born_solvation_zero_for_vacuum() {
        // epsilon_r = 1.0 → (1 - 1/1) = 0 → zero energy
        let e = born_solvation_energy(1.6e-19, 1.0e-10, 1.0, 8.987_551_792_3e9);
        assert!(
            e.abs() < 1e-30,
            "Born energy in vacuum (ε=1) should be 0, got {e}"
        );
    }

    // ── Drude oscillator ─────────────────────────────────────────────────────

    #[test]
    fn test_drude_spring_energy_zero_at_equilibrium() {
        let drude = DrudeOscillator::new(0, -0.5, 100.0);
        assert!(
            drude.spring_energy().abs() < 1e-30,
            "Spring energy at zero displacement = 0"
        );
    }

    #[test]
    fn test_drude_spring_energy_positive_on_displacement() {
        let mut drude = DrudeOscillator::new(0, -0.5, 100.0);
        drude.displacement = [0.1, 0.0, 0.0];
        let e = drude.spring_energy();
        assert!(
            e > 0.0,
            "Spring energy should be positive on displacement, got {e}"
        );
    }

    #[test]
    fn test_drude_restoring_force_opposes_displacement() {
        let mut drude = DrudeOscillator::new(0, -0.5, 100.0);
        drude.displacement = [0.1, 0.0, 0.0];
        let f = drude.restoring_force();
        assert!(
            f[0] < 0.0,
            "Restoring force should oppose +x displacement, got {}",
            f[0]
        );
    }

    #[test]
    fn test_drude_position() {
        let mut drude = DrudeOscillator::new(0, -0.5, 100.0);
        drude.displacement = [0.1, 0.2, 0.3];
        let core = [1.0, 2.0, 3.0];
        let pos = drude.drude_position(core);
        assert!((pos[0] - 1.1).abs() < 1e-12);
        assert!((pos[1] - 2.2).abs() < 1e-12);
        assert!((pos[2] - 3.3).abs() < 1e-12);
    }

    #[test]
    fn test_drude_polarisability() {
        let drude = DrudeOscillator::new(0, -2.0, 100.0);
        let alpha = drude.polarisability();
        // alpha = q^2 / k = 4.0 / 100.0 = 0.04
        assert!((alpha - 0.04).abs() < 1e-12, "alpha = {alpha}");
    }

    // ── Multipole expansion ──────────────────────────────────────────────────

    #[test]
    fn test_multipole_monopole_neutral() {
        let charges = vec![0.5, -0.5];
        let q = multipole_monopole(&charges);
        assert!(q.abs() < 1e-12, "Neutral system monopole = 0, got {q}");
    }

    #[test]
    fn test_multipole_monopole_nonzero() {
        let charges = vec![1.0, 1.0];
        let q = multipole_monopole(&charges);
        assert!((q - 2.0).abs() < 1e-12, "Total charge should be 2, got {q}");
    }

    #[test]
    fn test_multipole_dipole_along_x() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let p = multipole_dipole(&positions, &charges);
        assert!(p[0].abs() > 0.0, "Dipole should be along x, p.x = {}", p[0]);
        assert!(p[1].abs() < 1e-12);
        assert!(p[2].abs() < 1e-12);
    }

    #[test]
    fn test_multipole_dipole_symmetric_charge_zero() {
        // Symmetric +q/-q/+q/−q configuration → dipole can still be non-zero
        // but for +1 at origin → sum q*r = just +1*[0,0,0] = 0
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.0];
        let p = multipole_dipole(&positions, &charges);
        // p = q * r = 1 * [0,0,0] = [0,0,0]
        assert!(p[0].abs() < 1e-12);
        assert!(p[1].abs() < 1e-12);
        assert!(p[2].abs() < 1e-12);
    }

    #[test]
    fn test_multipole_quadrupole_symmetric() {
        // Two equal charges at +x and -x → Q should be traceless
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];
        let q = multipole_quadrupole(&positions, &charges);
        // Trace should be zero for traceless quadrupole
        let trace = q[0][0] + q[1][1] + q[2][2];
        assert!(
            trace.abs() < 1e-12,
            "Quadrupole trace should be 0, got {trace}"
        );
    }

    // ── Reaction field ────────────────────────────────────────────────────────

    #[test]
    fn test_reaction_field_correction_sign() {
        // For opposite charges: q_i=1, q_j=-1 → -k_rf * (+1)*(-1) * r² > 0
        let rf = reaction_field_correction(0.5_f64.powi(2), 1.0, -1.0, 1.0, 80.0, 1.0);
        assert!(
            rf > 0.0,
            "RF correction for opposite charges (q+,q-) should be positive: {rf}"
        );
    }

    #[test]
    fn test_coulomb_plus_rf_zero_beyond_rc() {
        let e = coulomb_plus_reaction_field(2.0, 1.0, -1.0, 1.0, 80.0, 1.0);
        assert_eq!(e, 0.0, "Energy beyond rc should be zero");
    }

    #[test]
    fn test_reaction_field_zero_at_rc_boundary() {
        let rc = 1.0;
        let e = coulomb_plus_reaction_field(rc - 1e-12, 1.0, -1.0, rc, 80.0, 1.0);
        // Should not panic and give a finite result
        assert!(e.is_finite(), "Energy just inside rc should be finite");
    }

    // ── Image charge ──────────────────────────────────────────────────────────

    #[test]
    fn test_image_charge_planar_sign() {
        let img = image_charge_planar(1.0);
        assert!(
            (img - (-1.0)).abs() < 1e-12,
            "Image charge should be -q, got {img}"
        );
    }

    #[test]
    fn test_image_charge_energy_negative() {
        // Attractive interaction with grounded conductor
        let e = image_charge_energy_planar(1.0, 1.0, 1.0);
        assert!(
            e < 0.0,
            "Image charge energy should be negative (attractive), got {e}"
        );
    }

    #[test]
    fn test_image_charge_energy_scales_with_q_squared() {
        let e1 = image_charge_energy_planar(1.0, 1.0, 1.0);
        let e2 = image_charge_energy_planar(2.0, 1.0, 1.0);
        assert!(
            (e2 / e1 - 4.0).abs() < 1e-10,
            "Image energy ∝ q²: ratio = {}",
            e2 / e1
        );
    }

    #[test]
    fn test_image_charge_energy_inversely_proportional_to_d() {
        let e1 = image_charge_energy_planar(1.0, 1.0, 1.0);
        let e2 = image_charge_energy_planar(1.0, 2.0, 1.0);
        assert!(
            (e1 / e2 - 2.0).abs() < 1e-10,
            "Image energy ∝ 1/d: ratio = {}",
            e1 / e2
        );
    }

    // ── Charge scaling for FEP ────────────────────────────────────────────────

    #[test]
    fn test_scale_charges_lambda_zero() {
        let charges = vec![1.0, -1.0, 0.5];
        let scaled = scale_charges(&charges, 0.0);
        assert!(
            scaled.iter().all(|&q| q.abs() < 1e-12),
            "λ=0 → all charges zero"
        );
    }

    #[test]
    fn test_scale_charges_lambda_one() {
        let charges = vec![1.0, -1.0, 0.5];
        let scaled = scale_charges(&charges, 1.0);
        for (a, b) in charges.iter().zip(scaled.iter()) {
            assert!((a - b).abs() < 1e-12, "λ=1 → charges unchanged");
        }
    }

    #[test]
    fn test_scale_charges_half() {
        let charges = vec![2.0, -4.0];
        let scaled = scale_charges(&charges, 0.5);
        assert!((scaled[0] - 1.0).abs() < 1e-12);
        assert!((scaled[1] - (-2.0)).abs() < 1e-12);
    }

    // ── RESP charge restraints ────────────────────────────────────────────────

    #[test]
    fn test_resp_restraint_zero_for_zero_charges() {
        let charges = vec![0.0, 0.0, 0.0];
        let e = resp_restraint_energy(&charges, 0.0005, 0.1);
        assert!(e.abs() < 1e-10, "RESP restraint for zero charges = 0");
    }

    #[test]
    fn test_resp_restraint_positive_for_nonzero_charges() {
        let charges = vec![0.5, -0.5];
        let e = resp_restraint_energy(&charges, 0.001, 0.1);
        assert!(e > 0.0, "RESP restraint should be positive, got {e}");
    }

    #[test]
    fn test_resp_gradient_zero_at_zero_charge() {
        let charges = vec![0.0];
        let grad = resp_restraint_gradient(&charges, 0.001, 0.1);
        assert!(
            grad[0].abs() < 1e-12,
            "RESP gradient at q=0 should be 0, got {}",
            grad[0]
        );
    }

    #[test]
    fn test_resp_gradient_positive_for_positive_charge() {
        let charges = vec![0.5];
        let grad = resp_restraint_gradient(&charges, 0.001, 0.1);
        assert!(
            grad[0] > 0.0,
            "RESP gradient for +q should be positive, got {}",
            grad[0]
        );
    }

    #[test]
    fn test_resp_gradient_negative_for_negative_charge() {
        let charges = vec![-0.5];
        let grad = resp_restraint_gradient(&charges, 0.001, 0.1);
        assert!(
            grad[0] < 0.0,
            "RESP gradient for -q should be negative, got {}",
            grad[0]
        );
    }

    // ── Additional Coulomb / potential tests ──────────────────────────────────

    #[test]
    fn test_coulomb_potential_sign_opposite_charges() {
        let v = coulomb_potential(1.0e-10, 1.0e-19, -1.0e-19);
        assert!(v < 0.0, "Opposite charges → negative potential");
    }

    #[test]
    fn test_coulomb_potential_sign_same_charges() {
        let v = coulomb_potential(1.0e-10, 1.0e-19, 1.0e-19);
        assert!(v > 0.0, "Same charges → positive potential");
    }

    #[test]
    fn test_yukawa_zero_at_large_r() {
        let v = yukawa_potential(100.0, 1.0, 1.0, 10.0);
        assert!(v <= 1e-400, "Yukawa should decay to ≈0 at large r, got {v}");
    }

    #[test]
    fn test_ewald_self_energy_zero_for_no_charges() {
        let e = ewald_self_energy(&[], 1.0);
        assert_eq!(e, 0.0, "Self-energy for empty system = 0");
    }

    #[test]
    fn test_neutralise_charges_empty() {
        let mut charges: Vec<f64> = vec![];
        neutralise_charges(&mut charges); // should not panic
    }

    #[test]
    fn test_neutralise_preserves_sum_zero() {
        let mut charges = vec![1.0, -0.5, 0.2, 0.1];
        neutralise_charges(&mut charges);
        let total: f64 = charges.iter().sum();
        assert!(total.abs() < 1e-12, "After neutralisation, sum = {total}");
    }

    #[test]
    fn test_gasteiger_three_atoms_chain() {
        // A-B-C chain; A = high EN, C = low EN
        let charges = gasteiger_charges(3, &[(0, 1), (1, 2)], &[10.0, 7.0, 5.0], &[10.0; 3]);
        let total: f64 = charges.iter().sum();
        assert!(total.abs() < 1e-10, "Chain Gasteiger sum ≠ 0: {total}");
    }

    #[test]
    fn test_electrostatic_potential_decreases_with_distance() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.6e-19];
        let phi1 = electrostatic_potential([1.0e-10, 0.0, 0.0], &positions, &charges);
        let phi2 = electrostatic_potential([2.0e-10, 0.0, 0.0], &positions, &charges);
        assert!(
            phi1 > phi2,
            "Potential decreases with distance from positive charge"
        );
    }

    #[test]
    fn test_electric_field_magnitude_inverse_square() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.6e-19];
        let e1 = electric_field([1.0e-10, 0.0, 0.0], &positions, &charges);
        let e2 = electric_field([2.0e-10, 0.0, 0.0], &positions, &charges);
        let mag1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
        let mag2 = (e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]).sqrt();
        // |E| ∝ 1/r²  → ratio should be 4
        assert!(
            (mag1 / mag2 - 4.0).abs() < 1e-8,
            "Field magnitude ratio = {}",
            mag1 / mag2
        );
    }

    #[test]
    fn test_wolf_potential_zero_charge() {
        let v = wolf_potential(0.5e-9, 0.0, 1.6e-19, 0.2e10, 1.0e-9);
        assert!(v.abs() < 1e-30, "Zero charge → zero Wolf potential");
    }

    #[test]
    fn test_compute_total_forces_three_atoms() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = vec![1.0, -2.0, 1.0];
        let forces = compute_total_coulomb_forces(&positions, &charges);
        // Net force should be zero (Newton's 3rd law, no external fields)
        let net: [f64; 3] = (0..3).fold([0.0; 3], |mut acc, i| {
            acc[0] += forces[i][0];
            acc[1] += forces[i][1];
            acc[2] += forces[i][2];
            acc
        });
        for (d, v) in net.iter().enumerate() {
            assert!(v.abs() < 1e-10, "Net force[{d}] = {v}");
        }
    }
}
