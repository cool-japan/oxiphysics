// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Stable Neo-Hookean XPBD constraint (Macklin & Muller 2021 / Smith et al. 2018).
//!
//! Splits the stable Neo-Hookean energy into a deviatoric distortion
//! constraint `C_dev = ||F||_F - sqrt(3)` and a hydrostatic volume constraint
//! `C_hyd = det(F) - alpha`, each projected with XPBD compliance so the
//! material is stable at near-incompressible Poisson ratios and through
//! inversion. Operates on raw `[[f64; 3]]` position arrays with a parallel
//! `inv_mass` slice so it is self-contained and easy to test.

use crate::fem_soft::math_helpers::{det3x3, edge_matrix_raw, inv3x3, mul3x3, transpose3x3};

/// Signed cofactor matrix of `f` (row-major), i.e. the derivative `dJ/dF` of the
/// determinant. Computed from 2x2 minors so it stays finite as `J -> 0`, which
/// is essential for stable behaviour through element inversion.
fn cofactor3(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            f[1][1] * f[2][2] - f[1][2] * f[2][1],
            -(f[1][0] * f[2][2] - f[1][2] * f[2][0]),
            f[1][0] * f[2][1] - f[1][1] * f[2][0],
        ],
        [
            -(f[0][1] * f[2][2] - f[0][2] * f[2][1]),
            f[0][0] * f[2][2] - f[0][2] * f[2][0],
            -(f[0][0] * f[2][1] - f[0][1] * f[2][0]),
        ],
        [
            f[0][1] * f[1][2] - f[0][2] * f[1][1],
            -(f[0][0] * f[1][2] - f[0][2] * f[1][0]),
            f[0][0] * f[1][1] - f[0][1] * f[1][0],
        ],
    ]
}

/// Frobenius-squared norm `I_C = ||F||_F^2 = sum_ij F[i][j]^2`.
fn frobenius_sq(f: [[f64; 3]; 3]) -> f64 {
    f.iter().flat_map(|row| row.iter()).map(|x| x * x).sum()
}

/// One scalar XPBD constraint to project: its current value `c_value`, the
/// gradient with respect to `F` (`grad_in_f = dC/dF`, 3x3 row-major) and the
/// XPBD compliance. Bundled so [`solve_constraint`] stays within an idiomatic
/// argument count.
struct ScalarConstraint {
    c_value: f64,
    grad_in_f: [[f64; 3]; 3],
    compliance: f64,
}

/// Apply a single XPBD scalar-constraint projection.
///
/// The per-vertex position gradients follow the standard FEM identity
/// `Hc = G * Dm_inv^T`, whose columns are the gradients of vertices 1..3 and
/// whose negated column-sum is the gradient of vertex 0. The XPBD multiplier
/// update then drives `positions` toward the constraint manifold. Implemented
/// as a free function so it borrows neither `&self` nor the accumulated
/// multiplier through `self`, keeping the borrow checker happy.
fn solve_constraint(
    indices: [usize; 4],
    rest_inv: [[f64; 3]; 3],
    positions: &mut [[f64; 3]],
    inv_mass: &[f64],
    constraint: &ScalarConstraint,
    dt: f64,
    lambda_acc: &mut f64,
) {
    // Hc = G * Dm_inv^T: columns hold the position gradients of vertices 1..3.
    let hc = mul3x3(constraint.grad_in_f, transpose3x3(rest_inv));

    let grad_x1 = [hc[0][0], hc[1][0], hc[2][0]];
    let grad_x2 = [hc[0][1], hc[1][1], hc[2][1]];
    let grad_x3 = [hc[0][2], hc[1][2], hc[2][2]];
    let grad_x0 = [
        -(grad_x1[0] + grad_x2[0] + grad_x3[0]),
        -(grad_x1[1] + grad_x2[1] + grad_x3[1]),
        -(grad_x1[2] + grad_x2[2] + grad_x3[2]),
    ];

    let grads = [grad_x0, grad_x1, grad_x2, grad_x3];

    let dt2 = dt * dt;
    if dt2 <= 1e-30 {
        return;
    }
    let alpha_tilde = constraint.compliance / dt2;

    // denom = sum_k w_k * |grad_k|^2 + alpha_tilde.
    let denom: f64 = indices
        .iter()
        .zip(grads.iter())
        .map(|(&idx, grad)| {
            let mag_sq = grad.iter().map(|g| g * g).sum::<f64>();
            inv_mass[idx] * mag_sq
        })
        .sum::<f64>()
        + alpha_tilde;

    if denom <= 1e-30 {
        return;
    }

    let d_lambda = (-constraint.c_value - alpha_tilde * *lambda_acc) / denom;
    *lambda_acc += d_lambda;

    for (&idx, grad) in indices.iter().zip(grads.iter()) {
        let w = inv_mass[idx];
        let pos = &mut positions[idx];
        for (component, &g) in pos.iter_mut().zip(grad.iter()) {
            *component += w * d_lambda * g;
        }
    }
}

/// Stable Neo-Hookean material constraint for a single tetrahedron, projected
/// via XPBD (extended position-based dynamics).
pub struct StableNeoHookeanXpbdConstraint {
    /// Vertex indices `[v0, v1, v2, v3]` into the position / inverse-mass slices.
    pub indices: [usize; 4],
    /// Inverse of the rest-configuration edge matrix `Dm^-1`.
    rest_inv: [[f64; 3]; 3],
    /// Rest volume of the tetrahedron.
    rest_volume: f64,
    /// Lame shear modulus `mu`.
    mu: f64,
    /// Lame first parameter `lambda`.
    lambda: f64,
    /// Accumulated XPBD multiplier for the deviatoric constraint.
    lambda_dev: f64,
    /// Accumulated XPBD multiplier for the hydrostatic constraint.
    lambda_hyd: f64,
}

impl StableNeoHookeanXpbdConstraint {
    /// Build a constraint from the rest configuration.
    ///
    /// `positions` provides the rest vertex positions; `young` is Young's
    /// modulus and `poisson` is Poisson's ratio, from which the Lame parameters
    /// `mu` and `lambda` are derived.
    pub fn new(indices: [usize; 4], positions: &[[f64; 3]], young: f64, poisson: f64) -> Self {
        let p0 = positions[indices[0]];
        let p1 = positions[indices[1]];
        let p2 = positions[indices[2]];
        let p3 = positions[indices[3]];

        let dm = edge_matrix_raw(p0, p1, p2, p3);
        let rest_volume = det3x3(dm).abs() / 6.0;
        let rest_inv = inv3x3(dm);

        let mu = young / (2.0 * (1.0 + poisson));
        let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));

        Self {
            indices,
            rest_inv,
            rest_volume,
            mu,
            lambda,
            lambda_dev: 0.0,
            lambda_hyd: 0.0,
        }
    }

    /// Reset both accumulated XPBD multipliers to zero. Call once per substep
    /// before the constraint-projection iterations.
    pub fn reset_lambda(&mut self) {
        self.lambda_dev = 0.0;
        self.lambda_hyd = 0.0;
    }

    /// Compute the deformation gradient `F = Ds * Dm^-1` for the given
    /// positions.
    pub fn deformation_gradient(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let p0 = positions[self.indices[0]];
        let p1 = positions[self.indices[1]];
        let p2 = positions[self.indices[2]];
        let p3 = positions[self.indices[3]];
        let ds = edge_matrix_raw(p0, p1, p2, p3);
        mul3x3(ds, self.rest_inv)
    }

    /// Compliance of the deviatoric constraint, `1 / (mu * rest_volume)`, with a
    /// large fallback when the product underflows.
    fn deviatoric_compliance(&self) -> f64 {
        if self.mu * self.rest_volume > 1e-30 {
            1.0 / (self.mu * self.rest_volume)
        } else {
            1e30
        }
    }

    /// Compliance of the hydrostatic constraint, `1 / (lambda * rest_volume)`,
    /// with a large fallback when the product underflows.
    fn hydrostatic_compliance(&self) -> f64 {
        if self.lambda * self.rest_volume > 1e-30 {
            1.0 / (self.lambda * self.rest_volume)
        } else {
            1e30
        }
    }

    /// Project only the deviatoric (distortion) constraint
    /// `C_dev = sqrt(I_C) - sqrt(3)`. Useful on its own for testing the
    /// deviatoric convergence in isolation.
    pub fn project_deviatoric_only(
        &mut self,
        positions: &mut [[f64; 3]],
        inv_mass: &[f64],
        dt: f64,
    ) {
        let f = self.deformation_gradient(positions);
        let i_c = frobenius_sq(f);
        let denom_ic = i_c.sqrt();
        if denom_ic < 1e-12 {
            return;
        }

        let c_dev = denom_ic - 3.0f64.sqrt();
        // dC_dev/dF = F / sqrt(I_C).
        let g_dev = [
            [f[0][0] / denom_ic, f[0][1] / denom_ic, f[0][2] / denom_ic],
            [f[1][0] / denom_ic, f[1][1] / denom_ic, f[1][2] / denom_ic],
            [f[2][0] / denom_ic, f[2][1] / denom_ic, f[2][2] / denom_ic],
        ];

        let constraint = ScalarConstraint {
            c_value: c_dev,
            grad_in_f: g_dev,
            compliance: self.deviatoric_compliance(),
        };
        let mut lambda_acc = self.lambda_dev;
        solve_constraint(
            self.indices,
            self.rest_inv,
            positions,
            inv_mass,
            &constraint,
            dt,
            &mut lambda_acc,
        );
        self.lambda_dev = lambda_acc;
    }

    /// Project the deviatoric constraint and then the hydrostatic volume
    /// constraint. `F` is recomputed before the hydrostatic solve because the
    /// deviatoric projection moves the vertices.
    pub fn project(&mut self, positions: &mut [[f64; 3]], inv_mass: &[f64], dt: f64) {
        // Deviatoric pass.
        self.project_deviatoric_only(positions, inv_mass, dt);

        // Hydrostatic pass: recompute F since positions moved.
        let f = self.deformation_gradient(positions);
        let alpha = 1.0 + self.mu / self.lambda - self.mu / (4.0 * self.lambda);
        let c_hyd = det3x3(f) - alpha;
        let g_hyd = cofactor3(f);

        let constraint = ScalarConstraint {
            c_value: c_hyd,
            grad_in_f: g_hyd,
            compliance: self.hydrostatic_compliance(),
        };
        let mut lambda_acc = self.lambda_hyd;
        solve_constraint(
            self.indices,
            self.rest_inv,
            positions,
            inv_mass,
            &constraint,
            dt,
            &mut lambda_acc,
        );
        self.lambda_hyd = lambda_acc;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem_soft::math_helpers::det3x3;

    #[test]
    fn test_cdev_converges() {
        let rest = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut c = StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &rest, 1.0e4, 0.3);
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        // Stretch the tet so the deviatoric constraint is violated.
        let mut pos = rest;
        pos[1][0] = 1.4;
        pos[2][1] = 1.3;
        pos[3][2] = 1.25;

        // Baseline deviatoric constraint value.
        let f0 = c.deformation_gradient(&pos);
        let ic0: f64 = f0.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
        let c_dev0 = ic0.sqrt() - 3.0f64.sqrt();
        assert!(
            c_dev0.abs() > 0.1,
            "test setup: deviatoric constraint should start violated: {c_dev0}"
        );

        let dt = 1.0e-2;
        for _ in 0..200 {
            c.reset_lambda();
            c.project_deviatoric_only(&mut pos, &inv_mass, dt);
        }
        // After many XPBD iterations the deviatoric constraint must be near zero.
        let f = c.deformation_gradient(&pos);
        let i_c: f64 = f.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
        let c_dev = i_c.sqrt() - 3.0f64.sqrt();
        assert!(c_dev.abs() < 0.05, "C_dev failed to converge: {c_dev}");
        // Volume should be finite and positive.
        assert!(det3x3(f) > 0.0);
    }

    #[test]
    fn test_full_project_stable() {
        // Running deviatoric + hydrostatic together must stay finite and reduce
        // the combined violation; volume stays positive.
        let rest = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut c = StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &rest, 1.0e4, 0.3);
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        let mut pos = rest;
        pos[1][0] = 1.3;
        pos[2][1] = 1.2;
        pos[3][2] = 1.15;
        let dt = 1.0e-2;
        for _ in 0..300 {
            c.reset_lambda();
            c.project(&mut pos, &inv_mass, dt);
        }
        for p in &pos {
            for v in p {
                assert!(v.is_finite(), "position must stay finite");
            }
        }
        let f = c.deformation_gradient(&pos);
        assert!(det3x3(f) > 0.0, "volume must stay positive: {}", det3x3(f));
    }
}
