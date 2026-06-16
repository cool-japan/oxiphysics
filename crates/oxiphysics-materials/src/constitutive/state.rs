// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Internal state types for the [`crate::constitutive`] models.

/// Internal state for J2 (von Mises) plasticity with isotropic hardening.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct J2State {
    /// Accumulated plastic strain (Voigt-6, engineering shear).
    pub plastic_strain: [f64; 6],
    /// Equivalent (accumulated) plastic strain — the scalar hardening variable.
    pub equiv_plastic_strain: f64,
}

/// Internal state for a generalized-Maxwell viscoelastic chain: one history
/// (partial) stress per Maxwell branch (each Voigt-6).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ViscoelasticState {
    /// Per-branch history stresses (length = number of Maxwell branches).
    pub branch_stress: Vec<[f64; 6]>,
}

impl ViscoelasticState {
    /// Create a virgin state sized for `n_branches` Maxwell elements (all zero).
    pub fn with_branches(n_branches: usize) -> Self {
        Self {
            branch_stress: vec![[0.0; 6]; n_branches],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j2_state_default_is_zero() {
        let s = J2State::default();
        assert_eq!(s.plastic_strain, [0.0; 6]);
        assert_eq!(s.equiv_plastic_strain, 0.0);
    }

    #[test]
    fn viscoelastic_state_default_is_empty() {
        let s = ViscoelasticState::default();
        assert!(s.branch_stress.is_empty());
    }

    #[test]
    fn viscoelastic_state_with_branches_sizes() {
        let s = ViscoelasticState::with_branches(3);
        assert_eq!(s.branch_stress.len(), 3);
        for b in &s.branch_stress {
            assert_eq!(*b, [0.0; 6]);
        }
    }
}
