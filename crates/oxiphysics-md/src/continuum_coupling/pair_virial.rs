// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Periodic (pair) virial stress of [`AtomisticToContin`].
//!
//! The per-atom virial `Σ_i r_i ⊗ f_i` needs absolute positions in one
//! non-periodic frame.  Inside a periodic box the stress has to be assembled from
//! pair separations under the minimum-image convention instead, which is what this
//! module provides.

use super::AtomisticToContin;
use crate::MdError;

/// Wrap one component of a separation into `[-L/2, L/2]`.
fn wrap_component(d: f64, length: f64) -> f64 {
    d - length * (d / length).round()
}

/// Minimum-image separation vector for an orthorhombic periodic box.
fn minimum_image(d: [f64; 3], box_lengths: [f64; 3]) -> [f64; 3] {
    [
        wrap_component(d[0], box_lengths[0]),
        wrap_component(d[1], box_lengths[1]),
        wrap_component(d[2], box_lengths[2]),
    ]
}

/// Volume of an orthorhombic box, checked to be finite and positive.
fn box_volume(box_lengths: [f64; 3]) -> Result<f64, MdError> {
    if box_lengths.iter().any(|l| !l.is_finite() || *l <= 0.0) {
        return Err(MdError::General(format!(
            "periodic box edge lengths must be finite and positive, got {box_lengths:?}"
        )));
    }
    Ok(box_lengths[0] * box_lengths[1] * box_lengths[2])
}

impl AtomisticToContin {
    /// Compute the potential (virial) stress of a **periodic** orthorhombic box
    /// from pair forces.
    ///
    /// σ_pot = −(1/V) Σ_{i<j} r_ij ⊗ f_ij,     V = L_x L_y L_z
    ///
    /// `pairs` lists every interacting pair once as `(i, j, f_ij)`, where `f_ij`
    /// (eV/Å) is the force on atom `i` due to atom `j` (the force on `j` is
    /// `−f_ij`), and `r_ij = r_i − r_j` is evaluated with the minimum-image
    /// convention in the box `box_lengths = [L_x, L_y, L_z]` (Å).  The sum is
    /// symmetrised, and the result is in eV/Å³ with the tension-positive sign of
    /// [`AtomisticToContin::kinetic_stress`]: a repulsive pair contributes a
    /// negative (compressive) normal stress and raises
    /// [`AtomisticToContin::pressure`].
    ///
    /// Because only separations enter, the result is unchanged by translating all
    /// atoms or by shifting any atom by a lattice vector.  The stress is also stored
    /// in [`AtomisticToContin::stress_voigt`].
    ///
    /// # Errors
    /// Returns [`MdError::General`] if an edge length is not finite and positive,
    /// a pair couples an atom with itself or an index is outside `positions`.
    pub fn pair_virial_stress(
        &mut self,
        pairs: &[(usize, usize, [f64; 3])],
        positions: &[[f64; 3]],
        box_lengths: [f64; 3],
    ) -> Result<[f64; 6], MdError> {
        let volume = box_volume(box_lengths)?;
        let mut s = [0.0_f64; 6];
        for &(i, j, f) in pairs {
            if i == j {
                return Err(MdError::General(format!(
                    "pair ({i}, {j}) couples an atom with itself"
                )));
            }
            let (Some(&r_i), Some(&r_j)) = (positions.get(i), positions.get(j)) else {
                return Err(MdError::General(format!(
                    "pair ({i}, {j}) refers to an atom outside the {} positions",
                    positions.len()
                )));
            };
            let r = minimum_image(
                [r_i[0] - r_j[0], r_i[1] - r_j[1], r_i[2] - r_j[2]],
                box_lengths,
            );
            s[0] -= r[0] * f[0];
            s[1] -= r[1] * f[1];
            s[2] -= r[2] * f[2];
            s[3] -= 0.5_f64 * (r[0] * f[1] + r[1] * f[0]);
            s[4] -= 0.5_f64 * (r[1] * f[2] + r[2] * f[1]);
            s[5] -= 0.5_f64 * (r[2] * f[0] + r[0] * f[2]);
        }
        for si in s.iter_mut() {
            *si /= volume;
        }
        self.stress_voigt = s;
        Ok(s)
    }

    /// Total Cauchy stress (eV/Å³) of a **periodic** orthorhombic box: the kinetic
    /// stress of [`AtomisticToContin::kinetic_stress`] plus the minimum-image pair
    /// virial of [`AtomisticToContin::pair_virial_stress`], both normalised by the
    /// box volume.
    ///
    /// For an ideal gas (no pairs) the pressure `−tr σ / 3` is `N k_B T / V`.
    ///
    /// # Errors
    /// See [`AtomisticToContin::pair_virial_stress`].
    pub fn total_pair_stress(
        &mut self,
        pairs: &[(usize, usize, [f64; 3])],
        positions: &[[f64; 3]],
        box_lengths: [f64; 3],
        velocities: &[[f64; 3]],
        masses: &[f64],
    ) -> Result<[f64; 6], MdError> {
        let volume = box_volume(box_lengths)?;
        let s_kin = self.kinetic_stress(velocities, masses, volume);
        let s_pot = self.pair_virial_stress(pairs, positions, box_lengths)?;
        let mut s_total = [0.0_f64; 6];
        for (total, (kin, pot)) in s_total.iter_mut().zip(s_kin.iter().zip(s_pot.iter())) {
            *total = kin + pot;
        }
        self.stress_voigt = s_total;
        Ok(s_total)
    }
}
