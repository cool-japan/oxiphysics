//! # DihedralForceField - Trait Implementations
//!
//! This module contains trait implementations for `DihedralForceField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ForceField`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::{PeriodicBox, distance_pbc};

use super::functions::ForceField;
use super::types::DihedralForceField;

impl Default for DihedralForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for DihedralForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        for dih in &self.dihedrals {
            let (r_ij, _) = distance_pbc(&atoms.positions[dih.i], &atoms.positions[dih.j], pbox);
            let (r_jk, _) = distance_pbc(&atoms.positions[dih.j], &atoms.positions[dih.k], pbox);
            let (r_kl, _) = distance_pbc(&atoms.positions[dih.k], &atoms.positions[dih.l], pbox);
            let n1 = r_ij.cross(&r_jk);
            let n2 = r_jk.cross(&r_kl);
            let n1_norm = n1.norm();
            let n2_norm = n2.norm();
            if n1_norm < 1e-15 || n2_norm < 1e-15 {
                continue;
            }
            let n1_hat = n1 / n1_norm;
            let n2_hat = n2 / n2_norm;
            let cos_phi = n1_hat.dot(&n2_hat).clamp(-1.0, 1.0);
            let m1 = n1_hat.cross(&n2_hat);
            let rjk_hat = r_jk / r_jk.norm().max(1e-15);
            let sin_phi = m1.dot(&rjk_hat);
            let phi = sin_phi.atan2(cos_phi);
            let n_f = dih.n as f64;
            energy += dih.v_n * (1.0 + (n_f * phi - dih.gamma).cos());
            let d_energy_dphi = -dih.v_n * n_f * (n_f * phi - dih.gamma).sin();
            let rjk_len = r_jk.norm().max(1e-15);
            let fi = n1 * (-d_energy_dphi / (n1_norm * n1_norm) * rjk_len);
            let fl = n2 * (d_energy_dphi / (n2_norm * n2_norm) * rjk_len);
            let r_ij_dot_rjk = r_ij.dot(&r_jk) / (rjk_len * rjk_len);
            let r_kl_dot_rjk = r_kl.dot(&r_jk) / (rjk_len * rjk_len);
            let fj = fi * (-1.0 + r_ij_dot_rjk) - fl * r_kl_dot_rjk;
            let fk = -(fi + fj + fl);
            atoms.forces[dih.i] += fi;
            atoms.forces[dih.j] += fj;
            atoms.forces[dih.k] += fk;
            atoms.forces[dih.l] += fl;
        }
        energy
    }
}
