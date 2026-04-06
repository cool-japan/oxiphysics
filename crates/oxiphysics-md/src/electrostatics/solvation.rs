#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Solvation models: Generalized Born, Born solvation, image charges,
//! and analytical Poisson-Boltzmann.

use super::coulomb::COULOMB_K;

// ---------------------------------------------------------------------------
// GeneralizedBorn
// ---------------------------------------------------------------------------

/// Generalized Born (GB) solvation model.
///
/// Approximates the electrostatic component of the solvation free energy
/// using effective Born radii. The GB energy for a pair of charges is:
///
/// ```text
/// E_GB = -K * (1/epsilon_in - 1/epsilon_out) * q_i * q_j / f_GB(r_ij, R_i, R_j)
/// ```
///
/// where `f_GB = sqrt(r^2 + R_i * R_j * exp(-r^2 / (4 * R_i * R_j)))`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GeneralizedBorn {
    /// Interior dielectric constant (protein/solute).
    pub epsilon_in: f64,
    /// Exterior dielectric constant (solvent, typically ~80 for water).
    pub epsilon_out: f64,
    /// Effective Born radii for each atom (angstrom).
    pub born_radii: Vec<f64>,
}

impl GeneralizedBorn {
    /// Create a new GB model.
    #[allow(dead_code)]
    pub fn new(epsilon_in: f64, epsilon_out: f64, born_radii: Vec<f64>) -> Self {
        Self {
            epsilon_in,
            epsilon_out,
            born_radii,
        }
    }

    /// The f_GB function for a pair.
    ///
    /// ```text
    /// f_GB(r, R_i, R_j) = sqrt(r^2 + R_i * R_j * exp(-r^2 / (4 * R_i * R_j)))
    /// ```
    #[allow(dead_code)]
    pub fn f_gb(r: f64, r_i: f64, r_j: f64) -> f64 {
        let ri_rj = r_i * r_j;
        if ri_rj < 1e-20 {
            return r;
        }
        let r2 = r * r;
        (r2 + ri_rj * (-r2 / (4.0 * ri_rj)).exp()).sqrt()
    }

    /// GB pair energy (kJ mol^-1).
    #[allow(dead_code)]
    pub fn pair_energy(&self, q_i: f64, q_j: f64, r: f64, i: usize, j: usize) -> f64 {
        let factor = 1.0 / self.epsilon_in - 1.0 / self.epsilon_out;
        let f = Self::f_gb(r, self.born_radii[i], self.born_radii[j]);
        if f < 1e-20 {
            return 0.0;
        }
        -COULOMB_K * factor * q_i * q_j / f
    }

    /// Self-energy (Born solvation energy) for a single atom.
    ///
    /// ```text
    /// E_self = -K * (1/epsilon_in - 1/epsilon_out) * q^2 / (2*R)
    /// ```
    #[allow(dead_code)]
    pub fn self_energy(&self, q: f64, i: usize) -> f64 {
        let r_born = self.born_radii[i];
        if r_born < 1e-20 {
            return 0.0;
        }
        let factor = 1.0 / self.epsilon_in - 1.0 / self.epsilon_out;
        -COULOMB_K * factor * q * q / (2.0 * r_born)
    }

    /// Total GB energy for a system of charges.
    #[allow(dead_code)]
    pub fn total_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let mut energy = 0.0;
        // Self terms
        for i in 0..n {
            energy += self.self_energy(charges[i], i);
        }
        // Pair terms
        for i in 0..n {
            for j in (i + 1)..n {
                let dr = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                energy += self.pair_energy(charges[i], charges[j], r, i, j);
            }
        }
        energy
    }
}

// ---------------------------------------------------------------------------
// ImageChargeMethod — handle conducting boundaries
// ---------------------------------------------------------------------------

/// Electrostatic method using image charges for a conducting planar boundary.
///
/// Places an image charge `-q_i` at the mirror position across the z=0 plane
/// for each real charge `q_i` at position `(x_i, y_i, z_i)`.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ImageChargeMethod {
    /// Position of the conducting plane (z coordinate, Å).
    pub plane_z: f64,
}

#[allow(dead_code)]
impl ImageChargeMethod {
    /// Create a new image-charge method with the conducting plane at `z = plane_z`.
    pub fn new(plane_z: f64) -> Self {
        Self { plane_z }
    }

    /// Compute the image positions and charges for a set of real charges.
    pub fn image_system(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        let n = positions.len().min(charges.len());
        let img_pos: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let mut p = positions[i];
                // Reflect z across plane_z
                p[2] = 2.0 * self.plane_z - positions[i][2];
                p
            })
            .collect();
        let img_charges: Vec<f64> = charges[..n].iter().map(|&q| -q).collect();
        (img_pos, img_charges)
    }

    /// Compute the total interaction energy (kJ mol⁻¹) between real charges
    /// and their images.
    pub fn image_interaction_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let (img_pos, img_charges) = self.image_system(positions, charges);
        let n = positions.len().min(charges.len());
        let mut energy = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                let mut d2 = 0.0_f64;
                for a in 0..3 {
                    d2 += (positions[i][a] - img_pos[j][a]).powi(2);
                }
                if d2 > 1e-20 {
                    energy += COULOMB_K * charges[i] * img_charges[j] / d2.sqrt();
                }
            }
        }
        // Factor of 1/2 to avoid double counting
        energy * 0.5
    }

    /// Force on atom `idx` due to its own image (self-image force).
    ///
    /// The self-image force only has a z-component: F_z = -COULOMB_K * q² / (4 * dz²).
    pub fn self_image_force_z(&self, pos: [f64; 3], charge: f64) -> f64 {
        let dz = (pos[2] - self.plane_z).abs();
        if dz < 1e-10 {
            return 0.0;
        }
        -COULOMB_K * charge * charge / (4.0_f64 * dz * dz)
    }
}

// ---------------------------------------------------------------------------
// ChargedSphericalCavity — Born solvation model
// ---------------------------------------------------------------------------

/// Born model: the free energy of transferring a charge `q` from vacuum into
/// a spherical cavity of radius `a` embedded in a dielectric `epsilon`.
///
/// ΔG_Born = -(COULOMB_K / 2) * q² * (1 - 1/ε) / a
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct BornSolvation {
    /// Dielectric constant of the solvent.
    pub epsilon: f64,
}

#[allow(dead_code)]
impl BornSolvation {
    /// Create a Born solvation model with solvent dielectric `epsilon`.
    pub fn new(epsilon: f64) -> Self {
        Self { epsilon }
    }

    /// Born solvation free energy (kJ mol⁻¹) for a sphere of radius `a` (Å)
    /// carrying charge `q` (e).
    pub fn solvation_energy(&self, q: f64, a: f64) -> f64 {
        if a < 1e-10 || self.epsilon < 1.0 {
            return 0.0;
        }
        -(COULOMB_K * 0.5_f64) * q * q * (1.0_f64 - 1.0_f64 / self.epsilon) / a
    }

    /// Desolvation penalty (kJ mol⁻¹): energy cost to remove a charge from solvent to vacuum.
    pub fn desolvation_penalty(&self, q: f64, a: f64) -> f64 {
        -self.solvation_energy(q, a)
    }

    /// Differential Born radius effect: change in solvation energy when radius
    /// changes from `a` to `a + da`.
    pub fn d_solvation_energy_da(&self, q: f64, a: f64) -> f64 {
        if a < 1e-10 || self.epsilon < 1.0 {
            return 0.0;
        }
        (COULOMB_K * 0.5_f64) * q * q * (1.0_f64 - 1.0_f64 / self.epsilon) / (a * a)
    }
}

// ---------------------------------------------------------------------------
// Poisson-Boltzmann implicit solvent (extended)
// ---------------------------------------------------------------------------

/// Analytical linearized Poisson-Boltzmann (ALPB) model.
///
/// Computes the electrostatic solvation energy using an analytical
/// approximation to the full PB equation.
///
/// Reference: Sigalov et al., J. Chem. Phys. 122, 094511 (2005).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AnalyticalPoissonBoltzmann {
    /// Solvent dielectric.
    pub epsilon_solvent: f64,
    /// Inverse Debye length kappa (Å^-1).
    pub kappa: f64,
    /// Effective Born radii per atom (Å).
    pub born_radii: Vec<f64>,
    /// Solute dielectric constant.
    pub epsilon_solute: f64,
}

#[allow(dead_code)]
impl AnalyticalPoissonBoltzmann {
    /// Create a new ALPB model.
    pub fn new(
        epsilon_solvent: f64,
        epsilon_solute: f64,
        kappa: f64,
        born_radii: Vec<f64>,
    ) -> Self {
        Self {
            epsilon_solvent,
            kappa,
            born_radii,
            epsilon_solute,
        }
    }

    /// ALPB pairwise energy correction between atoms i and j.
    ///
    /// Uses the modified GB formula with PB-like screening:
    ///   f_ALPB = sqrt(r^2 + R_i*R_j*exp(-r^2/(4*R_i*R_j)))
    ///   E_ALPB = -(1/eps_in - exp(-kappa*f)/eps_out) * qi*qj / f_ALPB
    pub fn pair_energy(&self, qi: f64, qj: f64, r: f64, i: usize, j: usize) -> f64 {
        let ri = self.born_radii[i];
        let rj = self.born_radii[j];
        let rirj = ri * rj;
        if rirj < 1e-20 {
            return 0.0;
        }
        let r2 = r * r;
        let f = (r2 + rirj * (-r2 / (4.0 * rirj)).exp()).sqrt();
        if f < 1e-20 {
            return 0.0;
        }
        let screening = (-self.kappa * f).exp() / self.epsilon_solvent;
        let factor = 1.0 / self.epsilon_solute - screening;
        -COULOMB_K * factor * qi * qj / f
    }

    /// Total ALPB solvation energy for a system of charges.
    pub fn total_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = charges
            .len()
            .min(positions.len())
            .min(self.born_radii.len());
        let mut e = 0.0;
        for i in 0..n {
            // Self-term
            let ri = self.born_radii[i];
            if ri > 1e-20 {
                let screening = (-self.kappa * ri).exp() / self.epsilon_solvent;
                e -= COULOMB_K * (1.0 / self.epsilon_solute - screening) * charges[i] * charges[i]
                    / (2.0 * ri);
            }
            for j in (i + 1)..n {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                e += self.pair_energy(charges[i], charges[j], r, i, j);
            }
        }
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gb_f_gb_equals_r_for_zero_radii() {
        let f = GeneralizedBorn::f_gb(5.0, 0.0, 0.0);
        assert!(
            (f - 5.0).abs() < 1e-12,
            "f_GB should equal r when radii are 0"
        );
    }

    #[test]
    fn test_gb_f_gb_at_zero_distance() {
        let f = GeneralizedBorn::f_gb(0.0, 2.0, 2.0);
        // f_GB(0, R, R) = sqrt(0 + R^2 * exp(0)) = R
        let expected = 2.0;
        assert!(
            (f - expected).abs() < 1e-12,
            "f_GB(0,2,2) should be 2.0, got {f}"
        );
    }

    #[test]
    fn test_gb_self_energy_negative() {
        let gb = GeneralizedBorn::new(1.0, 80.0, vec![1.5]);
        let e = gb.self_energy(1.0, 0);
        // (1/1 - 1/80) * K * q^2 / (2*R) -- should be negative
        assert!(e < 0.0, "GB self-energy should be negative, got {e}");
    }

    #[test]
    fn test_gb_pair_energy_opposite_charges() {
        let gb = GeneralizedBorn::new(1.0, 80.0, vec![1.5, 1.5]);
        let _positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = gb.pair_energy(charges[0], charges[1], 3.0, 0, 1);
        // Opposite charges with (1/eps_in - 1/eps_out) > 0 and q*q < 0: positive (unfavorable)
        // Actually: factor = 1/1 - 1/80 ~ 0.9875, E = -K * 0.9875 * 1*(-1) / f = +
        assert!(
            e > 0.0,
            "GB pair energy for opposite charges should be positive (unfavorable in vacuum)"
        );
    }

    #[test]
    fn test_gb_total_energy_finite() {
        let gb = GeneralizedBorn::new(1.0, 80.0, vec![1.5, 1.5]);
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = gb.total_energy(&positions, &charges);
        assert!(e.is_finite(), "GB total energy should be finite, got {e}");
    }

    #[test]
    fn test_gb_self_energy_single_ion_negative() {
        // Born solvation energy for a single ion should be negative (solvation is favourable)
        // epsilon_in < epsilon_out => factor > 0 => self_energy < 0
        let gb = GeneralizedBorn::new(1.0, 80.0, vec![2.0]);
        let e = gb.self_energy(1.0, 0);
        assert!(
            e < 0.0,
            "GB self-energy for single ion should be negative, got {e}"
        );
    }

    #[test]
    fn test_gb_self_energy_zero_charge_zero() {
        let gb = GeneralizedBorn::new(1.0, 80.0, vec![2.0]);
        let e = gb.self_energy(0.0, 0);
        assert_eq!(e, 0.0, "zero charge should give zero GB self-energy");
    }

    #[test]
    fn test_gb_total_energy_finite_for_two_ions() {
        let gb = GeneralizedBorn::new(1.0, 80.0, vec![2.0, 2.0]);
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = gb.total_energy(&positions, &charges);
        assert!(e.is_finite(), "GB total energy should be finite, got {e}");
    }

    #[test]
    fn test_gb_total_energy_depends_on_epsilon_out() {
        let gb1 = GeneralizedBorn::new(1.0, 80.0, vec![2.0]);
        let gb2 = GeneralizedBorn::new(1.0, 40.0, vec![2.0]);
        let positions = [[0.0, 0.0, 0.0]];
        let charges = [1.0];
        let e1 = gb1.total_energy(&positions, &charges);
        let e2 = gb2.total_energy(&positions, &charges);
        // Higher epsilon_out -> more negative self-energy
        assert!(
            e1 < e2,
            "higher epsilon_out should give more negative GB energy: {e1} vs {e2}"
        );
    }

    #[test]
    fn test_gb_f_gb_at_zero_r() {
        let f = GeneralizedBorn::f_gb(0.0, 2.0, 2.0);
        // f_GB(0, R, R) = sqrt(R^2 * exp(0)) = R
        assert!((f - 2.0).abs() < 1e-12, "f_GB(0, R, R) = R, got {f}");
    }

    #[test]
    fn test_gb_f_gb_large_r_approaches_r() {
        // For large r, f_GB(r, R, R) -> r
        let r = 100.0;
        let f = GeneralizedBorn::f_gb(r, 2.0, 2.0);
        let rel = (f - r).abs() / r;
        assert!(
            rel < 1e-3,
            "f_GB for large r should approach r: f={f}, r={r}"
        );
    }

    #[test]
    fn test_image_charge_method_positions_mirrored() {
        let icm = ImageChargeMethod::new(0.0);
        let positions = [[0.0, 0.0, 5.0]];
        let charges = [1.0];
        let (img_pos, _) = icm.image_system(&positions, &charges);
        assert!(
            (img_pos[0][2] - (-5.0)).abs() < 1e-12,
            "image z = {}",
            img_pos[0][2]
        );
    }

    #[test]
    fn test_image_charge_method_charges_negated() {
        let icm = ImageChargeMethod::new(0.0);
        let positions = [[0.0, 0.0, 3.0]];
        let charges = [2.0];
        let (_, img_charges) = icm.image_system(&positions, &charges);
        assert!(
            (img_charges[0] - (-2.0)).abs() < 1e-12,
            "image charge = {}",
            img_charges[0]
        );
    }

    #[test]
    fn test_image_charge_method_interaction_energy_negative_for_positive_charge() {
        let icm = ImageChargeMethod::new(0.0);
        // Positive charge above plane: attracted to negative image below
        let positions = [[0.0, 0.0, 5.0]];
        let charges = [1.0];
        let e = icm.image_interaction_energy(&positions, &charges);
        assert!(
            e < 0.0,
            "positive charge + negative image → attractive energy: {e}"
        );
    }

    #[test]
    fn test_image_charge_method_self_image_force_attractive() {
        let icm = ImageChargeMethod::new(0.0);
        let f_z = icm.self_image_force_z([0.0, 0.0, 3.0], 1.0);
        // Force should be negative (toward plane)
        assert!(
            f_z < 0.0,
            "self-image force should pull charge toward plane: {f_z}"
        );
    }

    #[test]
    fn test_image_charge_method_self_image_force_decreases_with_distance() {
        let icm = ImageChargeMethod::new(0.0);
        let f1 = icm.self_image_force_z([0.0, 0.0, 1.0], 1.0).abs();
        let f2 = icm.self_image_force_z([0.0, 0.0, 5.0], 1.0).abs();
        assert!(
            f1 > f2,
            "closer charge has stronger image force: f1={f1}, f2={f2}"
        );
    }

    #[test]
    fn test_born_solvation_negative_for_water() {
        let bs = BornSolvation::new(80.0);
        let e = bs.solvation_energy(1.0, 1.5);
        assert!(e < 0.0, "Born solvation should be negative: {e}");
    }

    #[test]
    fn test_born_solvation_vacuum_gives_zero() {
        let bs = BornSolvation::new(1.0); // ε = 1 → vacuum
        let e = bs.solvation_energy(1.0, 1.5);
        assert_eq!(e, 0.0, "Born solvation in vacuum (ε=1) should be 0");
    }

    #[test]
    fn test_born_solvation_desolvation_penalty_positive() {
        let bs = BornSolvation::new(80.0);
        let penalty = bs.desolvation_penalty(1.0, 1.5);
        assert!(
            penalty > 0.0,
            "desolvation penalty must be positive: {penalty}"
        );
    }

    #[test]
    fn test_born_solvation_larger_radius_less_stabilized() {
        let bs = BornSolvation::new(80.0);
        let e1 = bs.solvation_energy(1.0, 1.0).abs();
        let e2 = bs.solvation_energy(1.0, 3.0).abs();
        assert!(e1 > e2, "smaller radius → more stabilized: {e1} vs {e2}");
    }

    #[test]
    fn test_born_solvation_d_energy_da_positive() {
        let bs = BornSolvation::new(80.0);
        let de = bs.d_solvation_energy_da(1.0, 2.0);
        assert!(
            de > 0.0,
            "dG/da should be positive (energy increases with radius): {de}"
        );
    }

    #[test]
    fn test_alpb_total_energy_finite() {
        let alpb = AnalyticalPoissonBoltzmann::new(80.0, 1.0, 0.1, vec![1.5, 1.5]);
        let positions = [[0.0; 3], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = alpb.total_energy(&positions, &charges);
        assert!(e.is_finite(), "ALPB energy = {e}");
    }

    #[test]
    fn test_alpb_opposite_charges_negative() {
        let alpb = AnalyticalPoissonBoltzmann::new(80.0, 1.0, 0.0, vec![1.5, 1.5]);
        let positions = [[0.0; 3], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = alpb.total_energy(&positions, &charges);
        assert!(e < 0.0, "opposite charges → negative energy: {e}");
    }
}
