//! Saint Venant-Kirchhoff tetrahedral element: the material model behind
//! [`FEMGpuKernel`](super::types::FEMGpuKernel).
//!
//! With the edge matrices `Ds = [x_b − x_a | x_c − x_a | x_d − x_a]` (current)
//! and `Dm` (rest) the deformation gradient is `F = Ds · Dm⁻¹`.  The model is
//!
//! * Green-Lagrange strain `E = ½ (FᵀF − I)`,
//! * energy density `ψ = μ E:E + (λ/2) tr(E)²`, i.e. element energy `V₀ ψ`,
//! * second Piola-Kirchhoff stress `S = 2μ E + λ tr(E) I` and first
//!   Piola-Kirchhoff stress `P = F S = ∂ψ/∂F`,
//! * nodal forces `f_b, f_c, f_d` = the columns of `−V₀ · P · Dm⁻ᵀ` and
//!   `f_a = −(f_b + f_c + f_d)`.
//!
//! `E` does not change under a rigid rotation of the element, so the rest pose
//! and every rigid motion are force free, and the forces are the exact negative
//! gradient of the energy.  Like every Saint Venant-Kirchhoff model it is only
//! trustworthy for moderate compression: the energy has a maximum at the fully
//! collapsed state, so the restoring force of an element squeezed towards zero
//! volume vanishes and an inverted element is not pushed back reliably.

use super::corotational::{self, Mat3};

/// A tetrahedron's rest-shape data and material parameters.
#[derive(Debug, Clone, Copy)]
pub(super) struct StvkTet {
    /// `Dm⁻ᵀ` of the rest tetrahedron (row-major).
    dm_inv_t: Mat3,
    /// Rest volume `V₀`, the weight of the strain energy density.
    rest_volume: f64,
    /// First Lamé parameter `λ`.
    lambda: f64,
    /// Shear modulus `μ`.
    mu: f64,
}

impl StvkTet {
    /// Element with the rest vertices `rest = [a, b, c, d]`, rest volume
    /// `rest_volume` (`V₀`) and Lamé parameters `(lambda, mu)`.
    ///
    /// Returns `None` when the element cannot carry a meaningful energy: the
    /// rest tetrahedron is degenerate (see [`corotational::rest_shape`]),
    /// `rest_volume` is not a positive finite number, or a Lamé parameter is
    /// not finite.
    pub(super) fn new(
        rest: &[[f64; 3]; 4],
        rest_volume: f64,
        lambda: f64,
        mu: f64,
    ) -> Option<Self> {
        let valid =
            rest_volume > 0.0 && rest_volume.is_finite() && lambda.is_finite() && mu.is_finite();
        if !valid {
            return None;
        }
        let (_, dm_inv_t) = corotational::rest_shape(rest)?;
        Some(Self {
            dm_inv_t,
            rest_volume,
            lambda,
            mu,
        })
    }

    /// Green-Lagrange strain `E = ½ (FᵀF − I)`.
    fn green_lagrange(f: &Mat3) -> Mat3 {
        let mut strain = corotational::mul(&corotational::transpose(f), f);
        for (i, row) in strain.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                let identity = if i == j { 1.0 } else { 0.0 };
                *entry = 0.5 * (*entry - identity);
            }
        }
        strain
    }

    /// Strain energy `V₀ [ μ E:E + (λ/2) tr(E)² ]` (J) at the current vertex
    /// positions `current = [a, b, c, d]`.
    pub(super) fn energy(&self, current: &[[f64; 3]; 4]) -> f64 {
        let f = corotational::deformation_gradient(&self.dm_inv_t, current);
        let strain = Self::green_lagrange(&f);
        let trace = corotational::trace(&strain);
        self.rest_volume
            * (self.mu * corotational::frobenius_sq(&strain) + 0.5 * self.lambda * trace * trace)
    }

    /// Nodal forces `[f_a, f_b, f_c, f_d]` (N): the negative gradient of
    /// [`energy`](Self::energy) with respect to the vertex positions.
    pub(super) fn forces(&self, current: &[[f64; 3]; 4]) -> [[f64; 3]; 4] {
        let f = corotational::deformation_gradient(&self.dm_inv_t, current);
        let strain = Self::green_lagrange(&f);
        let trace_term = self.lambda * corotational::trace(&strain);
        // Second Piola-Kirchhoff stress S = 2 mu E + lambda tr(E) I, built in place.
        let mut stress = strain;
        for (i, row) in stress.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry *= 2.0 * self.mu;
                if i == j {
                    *entry += trace_term;
                }
            }
        }
        let piola = corotational::mul(&f, &stress);
        corotational::nodal_forces(&piola, &self.dm_inv_t, self.rest_volume)
    }
}
