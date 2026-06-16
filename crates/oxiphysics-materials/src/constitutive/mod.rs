// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unified UMAT-style constitutive interface.
//!
//! A single [`ConstitutiveModel`] trait wraps the stress-update plus
//! consistent-tangent contract of every material-point model, so a finite
//! element can drive any model uniformly (mirroring the Abaqus UMAT contract).
//!
//! Voigt-6 convention: `[σ₁₁, σ₂₂, σ₃₃, σ₂₃, σ₁₃, σ₁₂]` with engineering shear
//! strains (the factor of two is absorbed into the convention).

mod state;
pub use state::{J2State, ViscoelasticState};

/// Algorithmic response of a single stress-update step.
#[derive(Debug, Clone)]
pub struct ConstitutiveResponse<S> {
    /// Updated Cauchy stress in Voigt-6 notation `[σ₁₁,σ₂₂,σ₃₃,σ₂₃,σ₁₃,σ₁₂]`.
    pub stress: [f64; 6],
    /// Consistent algorithmic tangent ∂σ/∂ε (6×6, Voigt-6).
    pub tangent: [[f64; 6]; 6],
    /// Updated internal state.
    pub state: S,
}

/// UMAT-style material-point constitutive model.
///
/// Implementors map a total strain at t_{n+1} plus the internal state at t_n
/// (and the step `dt`) to the updated stress, the consistent tangent, and the
/// updated state.
pub trait ConstitutiveModel {
    /// Internal state (plastic strain, history stresses, etc.).
    /// [`Default`] yields the virgin / unstressed material.
    type State: Clone + Default;

    /// Stress update over one increment.
    ///
    /// * `strain` — total strain at t_{n+1} (Voigt-6, engineering shear).
    /// * `state` — internal state at t_n.
    /// * `dt` — time-step size (s); ignored by rate-independent models.
    fn stress_update(
        &self,
        strain: &[f64; 6],
        state: &Self::State,
        dt: f64,
    ) -> ConstitutiveResponse<Self::State>;

    /// Number of scalar internal state variables (for FEM SDV storage sizing).
    fn n_state_vars(&self) -> usize;
}

/// Flatten a row-major `[f64; 36]` Voigt stiffness into a nested `[[f64; 6]; 6]`.
pub(crate) fn unflatten_6x6(flat: &[f64; 36]) -> [[f64; 6]; 6] {
    let mut out = [[0.0_f64; 6]; 6];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = flat[i * 6 + j];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unflatten_6x6_round_trips() {
        let mut flat = [0.0_f64; 36];
        for (i, c) in flat.iter_mut().enumerate() {
            *c = i as f64;
        }
        let nested = unflatten_6x6(&flat);
        for (i, row) in nested.iter().enumerate() {
            for (j, &cell) in row.iter().enumerate() {
                assert_eq!(cell, (i * 6 + j) as f64);
            }
        }
    }

    #[test]
    fn constitutive_response_constructs() {
        let r = ConstitutiveResponse {
            stress: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            tangent: [[0.0_f64; 6]; 6],
            state: 42_u32,
        };
        assert_eq!(r.stress[2], 3.0);
        assert_eq!(r.tangent[0][0], 0.0);
        assert_eq!(r.state, 42);
    }
}
