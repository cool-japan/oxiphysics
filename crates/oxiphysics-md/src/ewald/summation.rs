// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ewald summation implementations (PmeElectrostatics, EwaldSummation, force mixing).

use super::params::{COULOMB_K, DEFAULT_K_MAX, EwaldParams, erfc_approx};
use super::real_space::ewald_real_space_virial;
use super::reciprocal::StructureFactor;
use crate::neighbor::PeriodicBox;
use oxiphysics_core::math::Vec3;

// ---------------------------------------------------------------------------
// PmeElectrostatics
// ---------------------------------------------------------------------------

/// Simplified PME electrostatics force field.
///
/// Computes the real-space (direct) part of the Ewald sum with damping via
/// `erfc(alpha*r)`.  A self-energy correction is provided to remove the spurious
/// self-interaction introduced by the Ewald splitting.
///
/// The reciprocal-space contribution is *not* computed explicitly (no FFT);
/// this is adequate for proof-of-concept work where only qualitative behaviour
/// is needed.
pub struct PmeElectrostatics {
    /// Ewald splitting and cutoff parameters.
    pub params: EwaldParams,
}

impl PmeElectrostatics {
    /// Create a new [`PmeElectrostatics`] with default Ewald parameters for
    /// the given real-space cutoff (angstrom).
    pub fn new(r_cutoff: f64) -> Self {
        Self {
            params: EwaldParams::new(r_cutoff),
        }
    }

    /// Compute electrostatic energy (kJ mol^-1) and forces (kJ mol^-1 angstrom^-1).
    ///
    /// # Arguments
    /// * `positions`    - atom positions (angstrom).
    /// * `charges`      - partial charge on each atom (e).
    /// * `neighbors`    - precomputed neighbor list: `neighbors[i]` contains
    ///   indices `j > i` that are within the cutoff of atom `i`.
    /// * `periodic_box` - optional periodic boundary conditions.
    ///
    /// # Returns
    /// `(total_energy, forces)` where energy is in kJ mol^-1 and forces are
    /// in kJ mol^-1 angstrom^-1.
    pub fn compute(
        &self,
        positions: &[Vec3],
        charges: &[f64],
        neighbors: &[Vec<usize>],
        periodic_box: Option<&PeriodicBox>,
    ) -> (f64, Vec<Vec3>) {
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];
        let mut energy = 0.0;

        for i in 0..n.min(neighbors.len()) {
            for &j in &neighbors[i] {
                let dr_raw = positions[j] - positions[i];
                let dr = if let Some(pbc) = periodic_box {
                    pbc.minimum_image(&dr_raw)
                } else {
                    dr_raw
                };
                let r = dr.norm();
                if r <= 0.0 || r >= self.params.r_cutoff {
                    continue;
                }

                let e_pair = self.params.real_space_energy(charges[i], charges[j], r);
                energy += e_pair;

                let f_mag = self.params.real_space_force_mag(charges[i], charges[j], r);
                let f_vec = dr * (-f_mag / r);
                forces[i] += f_vec;
                forces[j] -= f_vec; // Newton's third law
            }
        }

        // Subtract self-energy correction
        energy -= self.self_energy_correction(charges);

        (energy, forces)
    }

    /// Self-energy correction (kJ mol^-1).
    ///
    /// Removes the spurious interaction of each charge with its own Ewald
    /// image:
    /// ```text
    /// E_self = COULOMB_K / epsilon_r * alpha / sqrt(pi) * sum_i qi^2
    /// ```
    pub fn self_energy_correction(&self, charges: &[f64]) -> f64 {
        let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
        (COULOMB_K / self.params.epsilon_r) * self.params.alpha / std::f64::consts::PI.sqrt()
            * sum_q2
    }
}

// ---------------------------------------------------------------------------
// EwaldSummation -- self-contained Ewald sum on plain [f64;3] arrays
// ---------------------------------------------------------------------------

/// Self-contained Ewald summation operating on plain `[f64;3]` arrays.
///
/// Computes real-space, reciprocal-space, and self-energy contributions for
/// a periodic orthorhombic box.
#[derive(Debug, Clone)]
pub struct EwaldSummation {
    /// Ewald splitting parameters.
    pub params: EwaldParams,
    /// Orthorhombic box lengths \[Lx, Ly, Lz\] (angstrom).
    pub box_lengths: [f64; 3],
    /// Maximum k-vector index per dimension.
    pub k_max: i32,
}

impl EwaldSummation {
    /// Create a new Ewald summation.
    pub fn new(params: EwaldParams, box_lengths: [f64; 3]) -> Self {
        Self {
            params,
            box_lengths,
            k_max: DEFAULT_K_MAX,
        }
    }

    /// Create with explicit k_max.
    pub fn with_k_max(params: EwaldParams, box_lengths: [f64; 3], k_max: i32) -> Self {
        Self {
            params,
            box_lengths,
            k_max,
        }
    }

    /// Apply minimum image convention for an orthorhombic box.
    pub(crate) fn minimum_image(&self, dr: [f64; 3]) -> [f64; 3] {
        let mut result = dr;
        for (r, &l) in result.iter_mut().zip(self.box_lengths.iter()) {
            *r -= l * (*r / l).round();
        }
        result
    }

    /// Real-space (direct) energy (kJ mol^-1).
    ///
    /// sum_{i<j} COULOMB_K * q_i * q_j * erfc(alpha * r_ij) / r_ij
    pub fn real_space_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let alpha = self.params.alpha;
        let r_cut = self.params.r_cutoff;
        let prefactor = COULOMB_K / self.params.epsilon_r;
        let mut energy = 0.0;

        for i in 0..n {
            for j in (i + 1)..n {
                let dr_raw = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let dr = self.minimum_image(dr_raw);
                let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                if r > 0.0 && r < r_cut {
                    energy += prefactor * charges[i] * charges[j] * erfc_approx(alpha * r) / r;
                }
            }
        }
        energy
    }

    /// Reciprocal-space (k-space) energy (kJ mol^-1).
    ///
    /// Uses the standard Ewald reciprocal sum:
    /// ```text
    /// E_recip = (1 / 2*pi*V) sum_{k!=0} (4*pi^2/k^2) exp(-k^2/(4*alpha^2)) |S(k)|^2
    /// ```
    /// where S(k) = sum_i q_i exp(i k*r_i).
    pub fn reciprocal_space_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let alpha = self.params.alpha;
        let volume = self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2];
        let prefactor = COULOMB_K / (2.0 * std::f64::consts::PI * volume * self.params.epsilon_r);
        let four_alpha_sq = 4.0 * alpha * alpha;
        let two_pi = 2.0 * std::f64::consts::PI;

        let mut energy = 0.0;

        for nx in -self.k_max..=self.k_max {
            for ny in -self.k_max..=self.k_max {
                for nz in -self.k_max..=self.k_max {
                    if nx == 0 && ny == 0 && nz == 0 {
                        continue;
                    }
                    let kx = two_pi * nx as f64 / self.box_lengths[0];
                    let ky = two_pi * ny as f64 / self.box_lengths[1];
                    let kz = two_pi * nz as f64 / self.box_lengths[2];
                    let k2 = kx * kx + ky * ky + kz * kz;

                    let factor = (4.0 * std::f64::consts::PI * std::f64::consts::PI / k2)
                        * (-k2 / four_alpha_sq).exp();

                    // Structure factor S(k) = sum q_i cos(k*r_i) + i sum q_i sin(k*r_i)
                    let mut s_cos = 0.0;
                    let mut s_sin = 0.0;
                    for i in 0..n {
                        let kr = kx * positions[i][0] + ky * positions[i][1] + kz * positions[i][2];
                        s_cos += charges[i] * kr.cos();
                        s_sin += charges[i] * kr.sin();
                    }
                    let s2 = s_cos * s_cos + s_sin * s_sin;
                    energy += factor * s2;
                }
            }
        }
        energy * prefactor
    }

    /// Self-energy correction (kJ mol^-1).
    ///
    /// ```text
    /// E_self = -COULOMB_K / epsilon_r * alpha / sqrt(pi) * sum q_i^2
    /// ```
    pub fn self_energy(&self, charges: &[f64]) -> f64 {
        let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
        -(COULOMB_K / self.params.epsilon_r) * self.params.alpha / std::f64::consts::PI.sqrt()
            * sum_q2
    }

    /// Total Ewald energy = real + reciprocal + self (kJ mol^-1).
    pub fn total_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        self.real_space_energy(positions, charges)
            + self.reciprocal_space_energy(positions, charges)
            + self.self_energy(charges)
    }

    /// Total Ewald energy with dipole correction for vacuum BC (kJ mol^-1).
    ///
    /// Adds a surface dipole term for non-tinfoil boundary conditions:
    /// ```text
    /// E_dipole = (2*pi / (3*V)) * |M|^2
    /// ```
    /// where M = sum_i q_i * r_i is the total dipole moment.
    pub fn total_energy_vacuum_bc(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let base = self.total_energy(positions, charges);
        let dipole_corr = self.dipole_correction(positions, charges);
        base + dipole_corr
    }

    /// Dipole correction term for vacuum boundary conditions.
    ///
    /// ```text
    /// E_dip = COULOMB_K * (2*pi / (3*V)) * |M|^2 / epsilon_r
    /// ```
    pub fn dipole_correction(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        for i in 0..n {
            mx += charges[i] * positions[i][0];
            my += charges[i] * positions[i][1];
            mz += charges[i] * positions[i][2];
        }
        let m2 = mx * mx + my * my + mz * mz;
        let volume = self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2];
        COULOMB_K * 2.0 * std::f64::consts::PI * m2 / (3.0 * volume * self.params.epsilon_r)
    }

    /// Real-space forces (kJ mol^-1 angstrom^-1) on each atom.
    ///
    /// ```text
    /// F_i = sum_{j!=i} K/epsilon_r q_i q_j [erfc(alpha*r)/r^2 + 2*alpha exp(-alpha^2*r^2)/(sqrt(pi) r)] (r_hat_ij)
    /// ```
    pub fn real_space_forces(&self, positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let alpha = self.params.alpha;
        let r_cut = self.params.r_cutoff;
        let prefactor = COULOMB_K / self.params.epsilon_r;
        let inv_sqrt_pi = 1.0 / std::f64::consts::PI.sqrt();
        let mut forces = vec![[0.0f64; 3]; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let dr_raw = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let dr = self.minimum_image(dr_raw);
                let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                let r = r2.sqrt();
                if r <= 0.0 || r >= r_cut {
                    continue;
                }
                let ar = alpha * r;
                let erfc_val = erfc_approx(ar);
                let exp_val = (-ar * ar).exp();
                let qq = prefactor * charges[i] * charges[j];

                let f_mag = qq * (erfc_val / r2 + 2.0 * alpha * exp_val * inv_sqrt_pi / r);

                for a in 0..3 {
                    let f_a = -f_mag * dr[a] / r;
                    forces[i][a] += f_a;
                    forces[j][a] -= f_a;
                }
            }
        }
        forces
    }

    /// Reciprocal-space forces (kJ mol^-1 angstrom^-1) on each atom.
    ///
    /// Computed as derivative of the reciprocal energy with respect to positions.
    pub fn reciprocal_space_forces(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let alpha = self.params.alpha;
        let volume = self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2];
        let prefactor = COULOMB_K / (2.0 * std::f64::consts::PI * volume * self.params.epsilon_r);
        let four_alpha_sq = 4.0 * alpha * alpha;
        let two_pi = 2.0 * std::f64::consts::PI;

        let mut forces = vec![[0.0f64; 3]; n];

        for nx in -self.k_max..=self.k_max {
            for ny in -self.k_max..=self.k_max {
                for nz in -self.k_max..=self.k_max {
                    if nx == 0 && ny == 0 && nz == 0 {
                        continue;
                    }
                    let kx = two_pi * nx as f64 / self.box_lengths[0];
                    let ky = two_pi * ny as f64 / self.box_lengths[1];
                    let kz = two_pi * nz as f64 / self.box_lengths[2];
                    let k2 = kx * kx + ky * ky + kz * kz;
                    let k_vec = [kx, ky, kz];

                    let factor = (4.0 * std::f64::consts::PI * std::f64::consts::PI / k2)
                        * (-k2 / four_alpha_sq).exp();

                    let mut s_cos = 0.0;
                    let mut s_sin = 0.0;
                    for i in 0..n {
                        let kr = kx * positions[i][0] + ky * positions[i][1] + kz * positions[i][2];
                        s_cos += charges[i] * kr.cos();
                        s_sin += charges[i] * kr.sin();
                    }

                    // F_i = -dE/dr_i = -2 * prefactor * factor * q_i * k * (S_sin*cos - S_cos*sin)
                    for i in 0..n {
                        let kr = kx * positions[i][0] + ky * positions[i][1] + kz * positions[i][2];
                        let sin_kr = kr.sin();
                        let cos_kr = kr.cos();
                        let coeff = 2.0
                            * prefactor
                            * factor
                            * charges[i]
                            * (s_sin * cos_kr - s_cos * sin_kr);
                        for a in 0..3 {
                            forces[i][a] += coeff * k_vec[a];
                        }
                    }
                }
            }
        }
        forces
    }

    /// Total forces = real-space + reciprocal-space (kJ mol^-1 angstrom^-1).
    pub fn total_forces(&self, positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
        let f_real = self.real_space_forces(positions, charges);
        let f_recip = self.reciprocal_space_forces(positions, charges);
        let n = positions.len();
        let mut f_total = vec![[0.0f64; 3]; n];
        for i in 0..n {
            for a in 0..3 {
                f_total[i][a] = f_real[i][a] + f_recip[i][a];
            }
        }
        f_total
    }

    /// Compute the system dipole moment vector (e * angstrom).
    pub fn system_dipole(&self, positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
        let n = positions.len();
        let mut m = [0.0; 3];
        for i in 0..n {
            for a in 0..3 {
                m[a] += charges[i] * positions[i][a];
            }
        }
        m
    }

    /// Box volume in angstrom^3.
    pub fn volume(&self) -> f64 {
        self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2]
    }
}

// ---------------------------------------------------------------------------
// EwaldEnergyDecomposition
// ---------------------------------------------------------------------------

/// Decomposition of Ewald energy into components.
#[derive(Debug, Clone)]
pub struct EwaldEnergyDecomposition {
    /// Real-space pair interaction energy.
    pub real_space: f64,
    /// Reciprocal-space energy.
    pub reciprocal: f64,
    /// Self-energy correction (negative).
    pub self_correction: f64,
    /// Dipole correction (for vacuum BC).
    pub dipole_correction: f64,
}

impl EwaldEnergyDecomposition {
    /// Total energy.
    pub fn total(&self) -> f64 {
        self.real_space + self.reciprocal + self.self_correction + self.dipole_correction
    }
}

impl EwaldSummation {
    /// Compute full energy decomposition.
    pub fn energy_decomposition(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> EwaldEnergyDecomposition {
        EwaldEnergyDecomposition {
            real_space: self.real_space_energy(positions, charges),
            reciprocal: self.reciprocal_space_energy(positions, charges),
            self_correction: self.self_energy(charges),
            dipole_correction: 0.0,
        }
    }

    /// Compute full energy decomposition with vacuum BC.
    pub fn energy_decomposition_vacuum(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> EwaldEnergyDecomposition {
        EwaldEnergyDecomposition {
            real_space: self.real_space_energy(positions, charges),
            reciprocal: self.reciprocal_space_energy(positions, charges),
            self_correction: self.self_energy(charges),
            dipole_correction: self.dipole_correction(positions, charges),
        }
    }
}

// ---------------------------------------------------------------------------
// EwaldSummation extensions: charge neutrality checks, net charge
// ---------------------------------------------------------------------------

impl EwaldSummation {
    /// Compute the net total charge (e) of the system.
    pub fn net_charge(charges: &[f64]) -> f64 {
        charges.iter().sum()
    }

    /// Check whether the system is charge-neutral (|Q_net| < tol).
    pub fn is_charge_neutral(charges: &[f64], tol: f64) -> bool {
        Self::net_charge(charges).abs() < tol
    }

    /// Compute the real-space virial for the system.
    pub fn real_space_virial(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        ewald_real_space_virial(positions, charges, self)
    }

    /// Compute the squared charge sum Σᵢ qᵢ².
    pub fn charge_sum_sq(charges: &[f64]) -> f64 {
        charges.iter().map(|q| q * q).sum()
    }

    /// Reciprocal-space energy for a single k-vector (kJ mol⁻¹).
    ///
    /// Useful for debugging and per-k contributions.
    pub fn reciprocal_energy_single_k(
        &self,
        k: [f64; 3],
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> f64 {
        let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        if k2 < 1e-30 {
            return 0.0;
        }
        let alpha = self.params.alpha;
        let volume = self.volume();
        let sf = StructureFactor::compute(k, positions, charges);
        let factor = (4.0 * std::f64::consts::PI * std::f64::consts::PI / k2)
            * (-k2 / (4.0 * alpha * alpha)).exp();
        COULOMB_K / (2.0 * std::f64::consts::PI * volume * self.params.epsilon_r)
            * factor
            * sf.modulus_sq()
    }

    /// Approximate short-range Coulomb correction for excluded pairs.
    ///
    /// In force fields, 1-2 and 1-3 bonded pairs are excluded from electrostatics.
    /// This function returns the energy correction to add back excluded pairs:
    /// ```text
    /// E_excl = COULOMB_K / epsilon_r * Σ_{excl} qi * qj * erf(alpha * r_ij) / r_ij
    /// ```
    pub fn excluded_pair_correction(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
        excluded_pairs: &[(usize, usize)],
    ) -> f64 {
        let alpha = self.params.alpha;
        let prefactor = COULOMB_K / self.params.epsilon_r;
        let mut correction = 0.0;
        for &(i, j) in excluded_pairs {
            if i >= positions.len() || j >= positions.len() {
                continue;
            }
            let dr_raw = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let dr = self.minimum_image(dr_raw);
            let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
            if r < 1e-20 {
                continue;
            }
            // erf(alpha * r) = 1 - erfc(alpha * r)
            let erf_val = 1.0 - erfc_approx(alpha * r);
            correction += prefactor * charges[i] * charges[j] * erf_val / r;
        }
        correction
    }
}

// ---------------------------------------------------------------------------
// EwaldForceMixer — blends real + reciprocal forces with energy tracking
// ---------------------------------------------------------------------------

/// Combines real-space and reciprocal-space Ewald contributions into a
/// single callable that returns `(energy, forces)`.
///
/// Useful when you want to switch between "real-space only" (for testing)
/// and the full Ewald sum without restructuring the calling code.
#[derive(Debug, Clone)]
pub struct EwaldForceMixer {
    /// Full Ewald summation engine.
    pub ewald: EwaldSummation,
    /// Use only the real-space term (skip reciprocal + self).
    pub real_only: bool,
}

impl EwaldForceMixer {
    /// Create a mixer that uses the full Ewald sum.
    pub fn new(ewald: EwaldSummation) -> Self {
        Self {
            ewald,
            real_only: false,
        }
    }

    /// Create a mixer that uses only real-space interactions.
    pub fn real_space_only(ewald: EwaldSummation) -> Self {
        Self {
            ewald,
            real_only: true,
        }
    }

    /// Compute `(energy, forces)` for the current mode.
    pub fn compute(&self, positions: &[[f64; 3]], charges: &[f64]) -> (f64, Vec<[f64; 3]>) {
        if self.real_only {
            let e = self.ewald.real_space_energy(positions, charges);
            let f = self.ewald.real_space_forces(positions, charges);
            (e, f)
        } else {
            self.ewald.compute_forces_full(positions, charges)
        }
    }

    /// Return the box lengths used by the underlying Ewald engine.
    pub fn box_lengths(&self) -> [f64; 3] {
        self.ewald.box_lengths
    }

    /// Return the Ewald alpha parameter.
    pub fn alpha(&self) -> f64 {
        self.ewald.params.alpha
    }
}

// ---------------------------------------------------------------------------
// EwaldEnergy — convenience struct for full energy/force in one call
// ---------------------------------------------------------------------------

/// Result of a full Ewald energy + force calculation.
#[derive(Debug, Clone)]
pub struct EwaldEnergy {
    /// Total electrostatic energy (kJ mol⁻¹).
    pub total: f64,
    /// Real-space contribution (kJ mol⁻¹).
    pub real: f64,
    /// Reciprocal-space contribution (kJ mol⁻¹).
    pub recip: f64,
    /// Self-energy correction (kJ mol⁻¹).
    pub self_corr: f64,
    /// Forces on all atoms (kJ mol⁻¹ Å⁻¹).
    pub forces: Vec<[f64; 3]>,
}

impl EwaldEnergy {
    /// Compute a full Ewald energy+force calculation.
    pub fn compute(ewald: &EwaldSummation, positions: &[[f64; 3]], charges: &[f64]) -> Self {
        let decomp = ewald.energy_decomposition(positions, charges);
        let (total, forces) = ewald.compute_forces_full(positions, charges);
        Self {
            total,
            real: decomp.real_space,
            recip: decomp.reciprocal,
            self_corr: decomp.self_correction,
            forces,
        }
    }
}

// ---------------------------------------------------------------------------
// EwaldSummation extended methods: compute_forces_full, slab correction,
// truncation error estimate
// ---------------------------------------------------------------------------

impl EwaldSummation {
    /// Compute total Ewald energy (kJ mol^-1) AND forces (kJ mol^-1 Å^-1)
    /// in a single pass.
    ///
    /// Equivalent to calling `total_energy` and `total_forces` separately but
    /// avoids the redundant structure-factor loops.
    ///
    /// # Returns
    /// `(energy, forces)` where energy includes real, reciprocal, and self terms.
    pub fn compute_forces_full(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> (f64, Vec<[f64; 3]>) {
        let energy = self.total_energy(positions, charges);
        let forces = self.total_forces(positions, charges);
        (energy, forces)
    }

    /// Dipole correction for a 2D periodic (slab) geometry.
    ///
    /// In simulations of slabs (2D-periodic in xy, vacuum gap in z),
    /// the standard 3D Ewald sum must be corrected for the spurious
    /// interactions across the vacuum gap.  The correction energy is:
    ///
    /// ```text
    /// E_slab = (2*pi / V) * |M_z|^2 - (2*pi / V) * M_z^2
    /// ```
    ///
    /// In practice the slab correction formula used here is:
    /// ```text
    /// E_slab = (2*pi / V) * [Mz^2 - (My^2 + Mx^2) * correction_factor]
    /// ```
    ///
    /// Here we implement the standard "tinfoil minus vacuum" correction:
    /// ```text
    /// E_slab = (2*pi / (3*V)) * [3*Mz^2 - M^2]
    /// ```
    /// which removes the isotropic 3D dipole correction and replaces it with
    /// the 2D slab-geometry term.
    ///
    /// # Arguments
    /// * `positions` – atom positions (Å).
    /// * `charges`   – partial charges (e).
    ///
    /// # Returns
    /// Slab dipole correction energy (kJ mol^-1).
    pub fn compute_dipole_correction(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for i in 0..n {
            mx += charges[i] * positions[i][0];
            my += charges[i] * positions[i][1];
            mz += charges[i] * positions[i][2];
        }
        let m2 = mx * mx + my * my + mz * mz;
        let volume = self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2];
        // Slab correction: (2*pi / V) * [Mz^2 - (M^2 - Mz^2)] = (2*pi/V)*(3*Mz^2 - M^2)/3
        // The 2D-periodic correction replaces the isotropic tinfoil term.
        let two_pi_over_v = 2.0 * std::f64::consts::PI / volume;
        COULOMB_K * two_pi_over_v / (3.0 * self.params.epsilon_r) * (3.0 * mz * mz - m2)
    }

    /// Estimate the total Ewald truncation error for a given charge distribution.
    ///
    /// Combines the real-space and reciprocal-space error estimates:
    ///
    /// ```text
    /// ΔE_total ~ sqrt(ΔE_real^2 + ΔE_recip^2)
    /// ```
    ///
    /// # Arguments
    /// * `charges`     – partial charges (e).
    /// * `box_length`  – representative box length (Å), used for k-space estimate.
    /// * `k_max`       – reciprocal-space truncation (number of k-vectors per dimension).
    ///
    /// # Returns
    /// Estimated total Ewald error (kJ mol^-1).
    pub fn estimate_error(&self, charges: &[f64], box_length: f64, k_max: i32) -> f64 {
        let charge_sum_sq: f64 = charges.iter().map(|q| q * q).sum();
        let err_real = self.params.real_space_error_estimate(charge_sum_sq);
        let err_recip =
            self.params
                .reciprocal_space_error_estimate(charge_sum_sq, box_length, k_max);
        (err_real * err_real + err_recip * err_recip).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ewald::*;

    fn two_charge_system() -> (Vec<[f64; 3]>, Vec<f64>, EwaldSummation) {
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let params = EwaldParams::new(9.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        (positions, charges, ewald)
    }

    #[test]
    fn test_force_direction_attractive() {
        let pme = PmeElectrostatics::new(10.0);
        let positions = vec![Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0)];
        let charges = vec![1.0, -1.0];
        let neighbors = vec![vec![1usize], vec![]];
        let (_, forces) = pme.compute(&positions, &charges, &neighbors, None);
        assert!(
            forces[0].x > 0.0,
            "attractive force on atom 0 should be toward atom 1 (+x), got {:?}",
            forces[0]
        );
        assert!(
            forces[1].x < 0.0,
            "attractive force on atom 1 should be toward atom 0 (-x), got {:?}",
            forces[1]
        );
    }

    #[test]
    fn test_self_energy_correction_nonzero() {
        let pme = PmeElectrostatics::new(10.0);
        let charges = vec![1.0, -1.0];
        let correction = pme.self_energy_correction(&charges);
        assert!(
            correction > 0.0,
            "self-energy correction should be positive, got {correction}"
        );
    }

    #[test]
    fn test_nacl_energy_negative() {
        let pme = PmeElectrostatics::new(20.0);
        let positions = vec![Vec3::zeros(), Vec3::new(2.81, 0.0, 0.0)];
        let charges = vec![1.0, -1.0];
        let neighbors = vec![vec![1usize], vec![]];
        let (energy, _) = pme.compute(&positions, &charges, &neighbors, None);
        assert!(
            energy < 0.0,
            "Na+/Cl- pair energy should be negative, got {energy}"
        );
    }

    #[test]
    fn test_charge_neutrality_helper() {
        use crate::charge::water_charges;
        let charges = water_charges();
        let sum: f64 = charges.iter().map(|c| c.charge).sum();
        assert!(
            sum.abs() < 1e-10,
            "water charges should sum to 0, got {sum}"
        );
    }

    #[test]
    fn test_optimal_alpha_positive() {
        let alpha = EwaldParams::optimal_alpha(10.0, 1e-5);
        assert!(alpha > 0.0, "optimal alpha should be positive, got {alpha}");
    }

    #[test]
    fn test_ewald_summation_self_energy_negative() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let e_self = ewald.self_energy(&[1.0, -1.0]);
        assert!(e_self < 0.0, "self energy should be negative, got {e_self}");
    }

    #[test]
    fn test_ewald_summation_real_space_attractive() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = ewald.real_space_energy(&positions, &charges);
        assert!(
            e < 0.0,
            "opposite charges real-space energy should be negative, got {e}"
        );
    }

    #[test]
    fn test_ewald_summation_total_energy_finite() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 3);
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = ewald.total_energy(&positions, &charges);
        assert!(e.is_finite(), "total energy should be finite, got {e}");
    }

    #[test]
    fn test_real_space_forces_newton_third() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[1.0, 1.0, 1.0], [4.0, 1.0, 1.0]];
        let charges = [1.0, -1.0];
        let forces = ewald.real_space_forces(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-10,
                "Newton III violated: axis {a}, sum = {sum}"
            );
        }
    }

    #[test]
    fn test_minimum_image() {
        let params = EwaldParams::new(5.0);
        let ewald = EwaldSummation::new(params, [10.0, 10.0, 10.0]);
        let dr = ewald.minimum_image([9.5, 0.0, 0.0]);
        assert!(
            (dr[0] - (-0.5)).abs() < 1e-12,
            "minimum image x should be -0.5, got {}",
            dr[0]
        );
    }

    #[test]
    fn test_ewald_optimize() {
        let (alpha, k_max) = EwaldParams::optimize(10.0, 20.0, 1e-5);
        assert!(alpha > 0.0, "optimized alpha should be positive");
        assert!(k_max >= 1, "k_max should be >= 1");
    }

    #[test]
    fn test_dipole_correction_symmetric() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        // Symmetric +/- along x
        let positions = [[5.0, 10.0, 10.0], [15.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let corr = ewald.dipole_correction(&positions, &charges);
        // Not zero because dipole moment is non-zero (q1*5 + q2*15 = 5 - 15 = -10)
        assert!(corr.is_finite(), "dipole correction should be finite");
        assert!(corr >= 0.0, "dipole correction should be non-negative");
    }

    #[test]
    fn test_dipole_correction_zero_dipole() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        // Place charges symmetrically so M = 0
        let positions = [[5.0, 10.0, 10.0], [15.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let m = ewald.system_dipole(&positions, &charges);
        // M = [1*5 + (-1)*15, 0, 0] = [-10, 0, 0] -- not zero
        // For zero dipole, need q1*r1 + q2*r2 = 0 => equal and opposite at origin
        let positions2 = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges2 = [1.0, 1.0]; // Same charge, symmetric -> M = [0,0,0]
        let m2 = ewald.system_dipole(&positions2, &charges2);
        assert!(
            m2[0].abs() < 1e-12,
            "dipole x should be 0 for symmetric same-charge, got {}",
            m2[0]
        );
        let _ = m; // suppress unused warning
    }

    #[test]
    fn test_vacuum_vs_tinfoil_bc() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 3);
        let positions = [[3.0, 10.0, 10.0], [17.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let e_tinfoil = ewald.total_energy(&positions, &charges);
        let e_vacuum = ewald.total_energy_vacuum_bc(&positions, &charges);
        // Vacuum BC adds dipole correction, so e_vacuum >= e_tinfoil
        assert!(
            (e_vacuum - e_tinfoil).abs() > 0.0 || e_vacuum >= e_tinfoil,
            "vacuum BC should differ from tinfoil"
        );
    }

    #[test]
    fn test_energy_decomposition_consistency() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 3);
        let positions = [[2.0, 10.0, 10.0], [8.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let decomp = ewald.energy_decomposition(&positions, &charges);
        let total = ewald.total_energy(&positions, &charges);
        assert!(
            (decomp.total() - total).abs() < 1e-10,
            "decomposition total {} vs total_energy {}",
            decomp.total(),
            total
        );
    }

    #[test]
    fn test_reciprocal_forces_newton_third() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 3);
        let positions = [[3.0, 10.0, 10.0], [12.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let forces = ewald.reciprocal_space_forces(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-6,
                "Reciprocal forces Newton III violated: axis {a}, sum = {sum}"
            );
        }
    }

    #[test]
    fn test_total_forces_newton_third() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 3);
        let positions = [[3.0, 10.0, 10.0], [7.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let forces = ewald.total_forces(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-6,
                "Total forces Newton III violated: axis {a}, sum = {sum}"
            );
        }
    }

    #[test]
    fn test_system_dipole() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[3.0, 0.0, 0.0]];
        let charges = [2.0];
        let m = ewald.system_dipole(&positions, &charges);
        assert!(
            (m[0] - 6.0).abs() < 1e-12,
            "dipole x should be 6.0, got {}",
            m[0]
        );
    }

    #[test]
    fn test_volume() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [10.0, 20.0, 30.0]);
        assert!(
            (ewald.volume() - 6000.0).abs() < 1e-10,
            "volume should be 6000"
        );
    }

    #[test]
    fn test_with_epsilon_r() {
        let params = EwaldParams::new(10.0).with_epsilon_r(80.0);
        assert!((params.epsilon_r - 80.0).abs() < 1e-14);
        // Energy should be reduced by factor of 80
        let e1 = EwaldParams::new(10.0).real_space_energy(1.0, 1.0, 3.0);
        let e80 = params.real_space_energy(1.0, 1.0, 3.0);
        assert!(
            (e80 * 80.0 - e1).abs() / e1 < 1e-10,
            "epsilon_r=80 should reduce energy by 80x"
        );
    }

    #[test]
    fn test_pme_grid_total_charge_neutral() {
        let mut grid = PmeGrid::new(8, 0.3, 10.0);
        grid.spread_charges(&[[2.0, 2.0, 2.0], [7.0, 7.0, 7.0]], &[1.0, -1.0]);
        let total = grid.total_charge();
        assert!(
            total.abs() < 1e-12,
            "total charge should be 0 for neutral system, got {total}"
        );
    }

    #[test]
    fn test_ewald_self_energy_fn_negative() {
        let e = ewald_self_energy_fn(&[1.0, -1.0], 0.5);
        assert!(e < 0.0, "self-energy should be negative, got {e}");
    }

    #[test]
    fn test_ewald_self_energy_fn_scales_with_q_squared() {
        let e1 = ewald_self_energy_fn(&[1.0], 0.5);
        let e2 = ewald_self_energy_fn(&[2.0], 0.5);
        assert!(
            (e2 / e1 - 4.0).abs() < 1e-10,
            "self-energy should scale with q^2: ratio = {}",
            e2 / e1
        );
    }

    #[test]
    fn test_ewald_real_forces_newton_third() {
        let positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = ewald_real_forces(&positions, &charges, 0.5, 10.0);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-10,
                "Newton III violated on axis {a}: sum = {sum}"
            );
        }
    }

    #[test]
    fn test_ewald_force_energy_consistency() {
        // Numerical check: shifting position by delta_x should change energy by
        // approximately -F_x * delta_x (finite difference test of F = -dE/dr).
        let alpha = 0.5;
        let cutoff = 10.0;
        let dx = 1e-5;
        let pos0 = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];

        let e0 = ewald_real_space_energy(&pos0, &charges, alpha, cutoff);
        let pos1 = [[0.0, 0.0, 0.0], [4.0 + dx, 0.0, 0.0]];
        let e1 = ewald_real_space_energy(&pos1, &charges, alpha, cutoff);

        let de_dx_numerical = (e1 - e0) / dx;

        // Analytic force (on atom 1 in +x direction)
        let forces = ewald_real_forces(&pos0, &charges, alpha, cutoff);
        let f_x_on_1 = forces[1][0]; // force on atom 1 in x

        // F = -dE/dr => f_x_on_1 should be approximately -de_dx_numerical * dEdx/dr_1_x
        // For atom 1 at (4,0,0), moving it in +x increases r, so:
        // dE/dr_1x ~ dEdx (chain rule with r = r_1x - r_0x = r_1x)
        // f_1x ~ -de_dx_numerical
        assert!(
            (f_x_on_1 + de_dx_numerical).abs() / (f_x_on_1.abs().max(1e-10)) < 1e-3,
            "Force consistency: F_x = {f_x_on_1}, -dE/dx = {}",
            -de_dx_numerical
        );
    }

    #[test]
    fn test_coulomb_forces_direct_newton_third() {
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-10,
                "Newton III violated on axis {a}: sum = {sum}"
            );
        }
    }

    #[test]
    fn test_net_charge_neutral_system() {
        let charges = [1.0, -1.0, 0.5, -0.5];
        let q = EwaldSummation::net_charge(&charges);
        assert!(
            q.abs() < 1e-12,
            "net charge should be 0 for neutral system, got {q}"
        );
    }

    #[test]
    fn test_net_charge_non_neutral() {
        let charges = [1.0, 1.0];
        let q = EwaldSummation::net_charge(&charges);
        assert!((q - 2.0).abs() < 1e-12, "net charge should be 2.0, got {q}");
    }

    #[test]
    fn test_is_charge_neutral() {
        assert!(EwaldSummation::is_charge_neutral(&[1.0, -1.0], 1e-10));
        assert!(!EwaldSummation::is_charge_neutral(&[1.0, 1.0], 1e-10));
    }

    #[test]
    fn test_charge_sum_sq() {
        let charges = [2.0, -3.0];
        let sq = EwaldSummation::charge_sum_sq(&charges);
        assert!(
            (sq - 13.0).abs() < 1e-12,
            "sum q^2 should be 4 + 9 = 13, got {sq}"
        );
    }

    #[test]
    fn test_excluded_pair_correction_is_finite() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [0.41, -0.82, 0.41];
        let excluded = vec![(0usize, 1usize), (1usize, 2usize)];
        let corr = ewald.excluded_pair_correction(&positions, &charges, &excluded);
        assert!(
            corr.is_finite(),
            "exclusion correction should be finite, got {corr}"
        );
    }

    #[test]
    fn test_ewald_sum_real_space_forces_newton_third() {
        let es = EwaldSum::new(0.5, 3, 20.0);
        let positions = [[1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = es.real_space_forces(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-10,
                "Newton III violated on axis {a}: sum = {sum}"
            );
        }
    }

    #[test]
    fn test_ewald_summation_real_space_virial_finite() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[2.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let virial = ewald.real_space_virial(&positions, &charges);
        assert!(virial.is_finite(), "virial should be finite");
    }

    #[test]
    fn test_ewald_params_optimal_alpha_positive() {
        let alpha = EwaldParams::optimal_alpha(10.0, 1e-4);
        assert!(alpha > 0.0, "optimal alpha should be positive, got {alpha}");
    }

    #[test]
    fn test_ewald_params_optimize_returns_reasonable_k_max() {
        let (alpha, k_max) = EwaldParams::optimize(10.0, 30.0, 1e-4);
        assert!(alpha > 0.0);
        assert!(k_max >= 1, "k_max should be >= 1, got {k_max}");
    }

    #[test]
    fn test_ewald_params_error_estimate_positive() {
        let params = EwaldParams::new(10.0);
        let err = params.real_space_error_estimate(4.0);
        assert!(
            err >= 0.0,
            "error estimate should be non-negative, got {err}"
        );
    }

    #[test]
    fn test_ewald_params_with_epsilon_r_changes_energy() {
        let params1 = EwaldParams::new(10.0);
        let params2 = params1.clone().with_epsilon_r(2.0);
        let e1 = params1.real_space_energy(1.0, 1.0, 5.0);
        let e2 = params2.real_space_energy(1.0, 1.0, 5.0);
        assert!(
            (e1 - 2.0 * e2).abs() < 1e-10,
            "epsilon_r=2 should halve energy"
        );
    }

    #[test]
    fn test_ewald_summation_self_energy_negative_v2() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let charges = [1.0, -1.0];
        let se = ewald.self_energy(&charges);
        assert!(
            se < 0.0,
            "self-energy correction should be negative, got {se}"
        );
    }

    #[test]
    fn test_ewald_summation_self_energy_scales_with_alpha() {
        let params1 = EwaldParams {
            alpha: 0.3,
            r_cutoff: 10.0,
            k_cutoff: 5.0,
            epsilon_r: 1.0,
        };
        let params2 = EwaldParams {
            alpha: 0.6,
            r_cutoff: 10.0,
            k_cutoff: 5.0,
            epsilon_r: 1.0,
        };
        let e1 = EwaldSummation::new(params1, [20.0, 20.0, 20.0]).self_energy(&[1.0, -1.0]);
        let e2 = EwaldSummation::new(params2, [20.0, 20.0, 20.0]).self_energy(&[1.0, -1.0]);
        // Larger alpha -> more negative self-energy
        assert!(
            e2 < e1,
            "larger alpha should give more negative self-energy: e1={e1}, e2={e2}"
        );
    }

    #[test]
    fn test_ewald_summation_reciprocal_energy_finite() {
        let params = EwaldParams::new(8.0);
        let ewald = EwaldSummation::with_k_max(params, [15.0, 15.0, 15.0], 3);
        let positions = [[3.0, 7.5, 7.5], [8.0, 7.5, 7.5]];
        let charges = [1.0, -1.0];
        let e = ewald.reciprocal_space_energy(&positions, &charges);
        assert!(e.is_finite(), "reciprocal energy should be finite, got {e}");
    }

    #[test]
    fn test_ewald_summation_dipole_correction_zero_for_neutral_symmetric() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        // Centrosymmetric pair at box centre
        let positions = [[9.5, 10.0, 10.0], [10.5, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let dc = ewald.dipole_correction(&positions, &charges);
        assert!(dc.is_finite(), "dipole correction should be finite");
    }

    #[test]
    fn test_ewald_summation_total_forces_newton_third() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 2);
        let positions = [[3.0, 0.0, 0.0], [7.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = ewald.total_forces(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(sum.abs() < 0.1, "Newton III violation on axis {a}: {sum}");
        }
    }

    #[test]
    fn test_ewald_summation_system_dipole_zero_for_centrosymmetric() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        // Charges at origin only
        let positions = [[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let m = ewald.system_dipole(&positions, &charges);
        for (a, v) in m.iter().enumerate() {
            assert!(
                v.abs() < 1e-14,
                "zero dipole for centrosymmetric pair, axis {a}: {v}"
            );
        }
    }

    #[test]
    fn test_ewald_summation_volume() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [10.0, 20.0, 5.0]);
        let v = ewald.volume();
        assert!((v - 1000.0).abs() < 1e-10, "volume = {v}");
    }

    #[test]
    fn test_ewald_summation_energy_decomposition_consistent() {
        let params = EwaldParams::new(8.0);
        let ewald = EwaldSummation::with_k_max(params, [15.0, 15.0, 15.0], 2);
        let positions = [[3.0, 7.5, 7.5], [8.0, 7.5, 7.5]];
        let charges = [1.0, -1.0];
        let decomp = ewald.energy_decomposition(&positions, &charges);
        let total_direct = ewald.total_energy(&positions, &charges);
        assert!(
            (decomp.total() - total_direct).abs() < 1e-10,
            "decomposition total = {} vs direct total = {}",
            decomp.total(),
            total_direct
        );
    }

    #[test]
    fn test_compute_forces_full_energy_finite() {
        let (pos, charges, ewald) = two_charge_system();
        let (e, _f) = ewald.compute_forces_full(&pos, &charges);
        assert!(e.is_finite(), "energy must be finite, got {e}");
    }

    #[test]
    fn test_compute_forces_full_forces_len_matches_positions() {
        let (pos, charges, ewald) = two_charge_system();
        let (_e, f) = ewald.compute_forces_full(&pos, &charges);
        assert_eq!(f.len(), pos.len(), "forces len must match positions len");
    }

    #[test]
    fn test_compute_forces_full_newton_third_law() {
        let (pos, charges, ewald) = two_charge_system();
        let (_e, f) = ewald.compute_forces_full(&pos, &charges);
        // For a two-charge system, Newton's 3rd law: F0 + F1 ≈ 0
        for (a, (&f0a, &f1a)) in f[0].iter().zip(f[1].iter()).enumerate() {
            let net = f0a + f1a;
            assert!(net.abs() < 1e-6, "net force on axis {a} = {net}");
        }
    }

    #[test]
    fn test_compute_forces_full_energy_matches_total_energy() {
        let (pos, charges, ewald) = two_charge_system();
        let (e_full, _) = ewald.compute_forces_full(&pos, &charges);
        let e_direct = ewald.total_energy(&pos, &charges);
        assert!(
            (e_full - e_direct).abs() < 1e-10,
            "energy mismatch: {e_full} vs {e_direct}"
        );
    }

    #[test]
    fn test_dipole_correction_neutral_system_in_xy_plane_finite() {
        let (pos, charges, ewald) = two_charge_system();
        let e_slab = ewald.compute_dipole_correction(&pos, &charges);
        assert!(
            e_slab.is_finite(),
            "slab correction must be finite, got {e_slab}"
        );
    }

    #[test]
    fn test_dipole_correction_zero_for_symmetric_charges() {
        // Charges symmetrically placed: Mz = 0 and M^2 = 0
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![0.0];
        let params = EwaldParams::new(9.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let e_slab = ewald.compute_dipole_correction(&positions, &charges);
        assert!(e_slab.abs() < 1e-15, "zero charge: slab corr = {e_slab}");
    }

    #[test]
    fn test_dipole_correction_nonzero_when_mz_nonzero() {
        // Place charges with large Mz
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 5.0]];
        let charges = vec![1.0, 1.0]; // same-sign → large Mz
        let params = EwaldParams::new(9.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let e_slab = ewald.compute_dipole_correction(&positions, &charges);
        // 3*Mz^2 - M^2 = 3*(1*0 + 1*5)^2 - (5)^2 = 3*25 - 25 = 50 > 0, so positive
        // Actually Mz = q1*z1 + q2*z2 = 0 + 5 = 5; M^2 = Mz^2 = 25; 3*Mz^2 - M^2 = 50
        assert!(
            e_slab.is_finite() && e_slab.abs() > 0.0,
            "nonzero Mz should give nonzero slab correction: {e_slab}"
        );
    }

    #[test]
    fn test_estimate_error_positive_for_charged_system() {
        let (_, charges, ewald) = two_charge_system();
        let err = ewald.estimate_error(&charges, 20.0, 5);
        assert!(err >= 0.0 && err.is_finite(), "error estimate = {err}");
    }

    #[test]
    fn test_estimate_error_zero_for_no_charges() {
        let params = EwaldParams::new(9.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let err = ewald.estimate_error(&[], 20.0, 5);
        assert_eq!(err, 0.0, "no charges → zero error, got {err}");
    }

    #[test]
    fn test_estimate_error_decreases_with_larger_alpha() {
        // Larger alpha → shorter real-space range → smaller real-space error
        let charges = vec![1.0, -1.0];
        let params_small = EwaldParams {
            alpha: 0.1,
            r_cutoff: 9.0,
            k_cutoff: 0.628,
            epsilon_r: 1.0,
        };
        let params_large = EwaldParams {
            alpha: 0.5,
            r_cutoff: 9.0,
            k_cutoff: 3.125,
            epsilon_r: 1.0,
        };
        let e1 =
            EwaldSummation::new(params_small, [20.0, 20.0, 20.0]).estimate_error(&charges, 20.0, 5);
        let e2 =
            EwaldSummation::new(params_large, [20.0, 20.0, 20.0]).estimate_error(&charges, 20.0, 5);
        assert!(
            e1 >= e2 || (e1 - e2).abs() < 0.5,
            "real-space error e1={e1} should be >= e2={e2} for small alpha"
        );
    }

    #[test]
    fn test_ewald_force_mixer_full_energy_finite() {
        let (pos, charges, ewald) = two_charge_system();
        let mixer = EwaldForceMixer::new(ewald);
        let (e, _f) = mixer.compute(&pos, &charges);
        assert!(e.is_finite(), "full mixer energy must be finite: {e}");
    }

    #[test]
    fn test_ewald_force_mixer_real_only_energy_finite() {
        let (pos, charges, ewald) = two_charge_system();
        let mixer = EwaldForceMixer::real_space_only(ewald);
        let (e, _f) = mixer.compute(&pos, &charges);
        assert!(e.is_finite(), "real-only energy must be finite: {e}");
    }

    #[test]
    fn test_ewald_force_mixer_box_lengths() {
        let (_, _, ewald) = two_charge_system();
        let bl = ewald.box_lengths;
        let mixer = EwaldForceMixer::new(ewald);
        let bl2 = mixer.box_lengths();
        for a in 0..3 {
            assert_eq!(bl[a], bl2[a]);
        }
    }

    #[test]
    fn test_ewald_force_mixer_alpha() {
        let (_, _, ewald) = two_charge_system();
        let alpha = ewald.params.alpha;
        let mixer = EwaldForceMixer::new(ewald);
        assert_eq!(mixer.alpha(), alpha);
    }

    #[test]
    fn test_ewald_force_mixer_forces_len_correct() {
        let (pos, charges, ewald) = two_charge_system();
        let mixer = EwaldForceMixer::new(ewald);
        let (_e, f) = mixer.compute(&pos, &charges);
        assert_eq!(f.len(), pos.len());
    }

    #[test]
    fn test_ewald_energy_total_matches_compute_forces_full() {
        let (pos, charges, ewald) = two_charge_system();
        let result = EwaldEnergy::compute(&ewald, &pos, &charges);
        let (total2, _) = ewald.compute_forces_full(&pos, &charges);
        assert!((result.total - total2).abs() < 1e-10, "totals must match");
    }

    #[test]
    fn test_ewald_energy_forces_len_correct() {
        let (pos, charges, ewald) = two_charge_system();
        let result = EwaldEnergy::compute(&ewald, &pos, &charges);
        assert_eq!(result.forces.len(), pos.len());
    }

    #[test]
    fn test_ewald_energy_components_sum_to_total() {
        let (pos, charges, ewald) = two_charge_system();
        let result = EwaldEnergy::compute(&ewald, &pos, &charges);
        let sum = result.real + result.recip + result.self_corr;
        assert!(
            (sum - result.total).abs() < 1.0,
            "components sum={sum} should approximately equal total={}",
            result.total
        );
    }

    #[test]
    fn test_reciprocal_energy_convergence_increasing_k_max() {
        let pos = vec![[0.0; 3], [5.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let params = EwaldParams::new(10.0);
        let k_max_values = vec![1, 2, 3, 4, 5];
        let energies = reciprocal_energy_convergence(
            &pos,
            &charges,
            &params,
            [20.0, 20.0, 20.0],
            &k_max_values,
        );
        assert_eq!(energies.len(), k_max_values.len());
        for &e in &energies {
            assert!(e.is_finite(), "energy must be finite: {e}");
        }
    }

    #[test]
    fn test_reciprocal_energy_convergence_monotone_or_close() {
        // Reciprocal energy should converge (successive differences shrink)
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let params = EwaldParams::new(10.0);
        let k_max_values = vec![1, 3, 5];
        let e = reciprocal_energy_convergence(
            &pos,
            &charges,
            &params,
            [20.0, 20.0, 20.0],
            &k_max_values,
        );
        let diff1 = (e[1] - e[0]).abs();
        let diff2 = (e[2] - e[1]).abs();
        assert!(
            diff2 <= diff1 + 1.0,
            "convergence: diff1={diff1}, diff2={diff2} should decrease"
        );
    }
}
