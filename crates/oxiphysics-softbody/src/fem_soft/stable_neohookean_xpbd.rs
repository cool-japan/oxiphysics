// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Stable Neo-Hookean XPBD constraint (Macklin & Muller 2021 / Smith et al. 2018).
//!
//! Splits the stable Neo-Hookean energy into a deviatoric distortion
//! constraint `C_dev = ||F||_F - sqrt(3)` and a hydrostatic volume constraint
//! `C_hyd = det(F) - 1`, each projected with XPBD compliance so the material is
//! stable at near-incompressible Poisson ratios and through inversion.
//!
//! Both constraints vanish at `F = I` (and at every rotation), so the
//! undeformed element is a fixed point of the projection. The `-sqrt(3)` offset
//! of `C_dev` is what fixes the hydrostatic target at `1`: Macklin's
//! `gamma = 1 + mu / lambda` (or Smith's `alpha`) belongs to a deviatoric term
//! without the offset and would leave a resting tet with `C_hyd != 0`, i.e.
//! inflate it (and diverge for `lambda -> 0`).
//!
//! Operates on raw `[[f64; 3]]` position arrays with a parallel `inv_mass` slice
//! so it is self-contained and easy to test.

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
    /// with a large fallback when the product underflows (`lambda = 0`, i.e. a
    /// Poisson ratio of zero, makes the volume constraint effectively inactive).
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

    /// Project only the hydrostatic (volume) constraint `C_hyd = det(F) - 1`.
    /// Useful on its own for testing the volume convergence in isolation.
    ///
    /// The deviatoric constraint `C_dev = ||F||_F - sqrt(3)` is already zero at
    /// `F = I`, so the volume constraint must be zero there too. A target
    /// `alpha != 1` (for example `1 + mu / lambda`) would make the rest shape
    /// violate `C_hyd`, inflate an undeformed tet, and become infinite for
    /// `lambda = 0`.
    pub fn project_hydrostatic_only(
        &mut self,
        positions: &mut [[f64; 3]],
        inv_mass: &[f64],
        dt: f64,
    ) {
        let f = self.deformation_gradient(positions);
        let c_hyd = det3x3(f) - 1.0;
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

    /// Project the deviatoric constraint and then the hydrostatic volume
    /// constraint. `F` is recomputed before the hydrostatic solve because the
    /// deviatoric projection moves the vertices.
    pub fn project(&mut self, positions: &mut [[f64; 3]], inv_mass: &[f64], dt: f64) {
        self.project_deviatoric_only(positions, inv_mass, dt);
        self.project_hydrostatic_only(positions, inv_mass, dt);
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

    const UNIT_TET: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    /// Skewed tetrahedron: neither `Dm` nor `Dm^-1` is symmetric.
    const SKEWED_TET: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [0.3, 1.1, 0.0],
        [0.2, 0.4, 0.9],
    ];

    fn max_displacement(a: &[[f64; 3]; 4], b: &[[f64; 3]; 4]) -> f64 {
        a.iter()
            .flatten()
            .zip(b.iter().flatten())
            .fold(0.0_f64, |acc, (x, y)| acc.max((x - y).abs()))
    }

    fn rotate(points: &[[f64; 3]; 4], axis: [f64; 3], angle: f64) -> [[f64; 3]; 4] {
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let k = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        let (s, c) = angle.sin_cos();
        points.map(|p| {
            // Rodrigues: p c + (k x p) s + k (k . p)(1 - c)
            let kxp = [
                k[1] * p[2] - k[2] * p[1],
                k[2] * p[0] - k[0] * p[2],
                k[0] * p[1] - k[1] * p[0],
            ];
            let kdp = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
            [
                p[0] * c + kxp[0] * s + k[0] * kdp * (1.0 - c),
                p[1] * c + kxp[1] * s + k[1] * kdp * (1.0 - c),
                p[2] * c + kxp[2] * s + k[2] * kdp * (1.0 - c),
            ]
        })
    }

    #[test]
    fn rest_tet_is_a_fixed_point_of_the_projection() {
        // C_dev and C_hyd both vanish at F = I, for every Poisson ratio. The old
        // target alpha = 1 + mu/lambda - mu/(4 lambda) (1.5 for nu = 0.3) made
        // C_hyd = -0.75 mu/lambda at rest and inflated the undeformed tet.
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        for poisson in [0.0, 0.1, 0.3, 0.45, 0.49] {
            for rest in [UNIT_TET, SKEWED_TET] {
                let mut c =
                    StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &rest, 1.0e4, poisson);
                let mut pos = rest;
                for _ in 0..20 {
                    c.reset_lambda();
                    c.project(&mut pos, &inv_mass, 1.0e-2);
                }
                let moved = max_displacement(&pos, &rest);
                assert!(
                    moved < 1.0e-12,
                    "nu = {poisson}: the undeformed tet moved by {moved}"
                );
            }
        }
    }

    #[test]
    fn rigidly_rotated_tet_is_also_a_fixed_point() {
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        let rotated = rotate(&SKEWED_TET, [1.0, 2.0, 3.0], 0.9);
        let mut c = StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &SKEWED_TET, 1.0e4, 0.3);
        let mut pos = rotated;
        for _ in 0..20 {
            c.reset_lambda();
            c.project(&mut pos, &inv_mass, 1.0e-2);
        }
        assert!(
            max_displacement(&pos, &rotated) < 1.0e-12,
            "a rigid rotation carries no strain"
        );
    }

    #[test]
    fn zero_poisson_ratio_stays_finite_and_barely_moves() {
        // nu = 0 gives lambda = 0: the old target alpha was infinite and every
        // position became inf/NaN. The hydrostatic constraint now simply falls
        // back to its (practically infinite) compliance.
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        let mut c = StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &UNIT_TET, 1.0e4, 0.0);
        let mut pos = UNIT_TET;
        for _ in 0..50 {
            c.reset_lambda();
            c.project(&mut pos, &inv_mass, 1.0e-2);
        }
        for p in &pos {
            assert!(p.iter().all(|v| v.is_finite()), "position must stay finite");
        }
        assert!(max_displacement(&pos, &UNIT_TET) < 1.0e-12);

        // A strained tet still relaxes through the deviatoric constraint.
        pos[1][0] = 1.3;
        pos[2][1] = 1.2;
        pos[3][2] = 1.15;
        for _ in 0..300 {
            c.reset_lambda();
            c.project(&mut pos, &inv_mass, 1.0e-2);
        }
        for p in &pos {
            assert!(p.iter().all(|v| v.is_finite()), "position must stay finite");
        }
        let f = c.deformation_gradient(&pos);
        let i_c: f64 = f.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
        assert!(
            (i_c.sqrt() - 3.0f64.sqrt()).abs() < 0.05,
            "C_dev must converge for nu = 0: {}",
            i_c.sqrt() - 3.0f64.sqrt()
        );
        assert!(det3x3(f) > 0.0);
    }

    fn scaled_tet(scale: f64) -> [[f64; 3]; 4] {
        UNIT_TET.map(|p| [scale * p[0], scale * p[1], scale * p[2]])
    }

    #[test]
    fn hydrostatic_pass_alone_restores_unit_volume() {
        // The volume constraint on its own is a single smooth constraint, so the
        // projection converges geometrically to `det F = 1`. The old target
        // alpha = 1 + mu/lambda - mu/(4 lambda) (1.5 for nu = 0.3) froze an
        // inflated tet at det F = 1.5 instead.
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        for (scale, start_det) in [(1.2, 1.728), (0.8, 0.512)] {
            let mut c = StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &UNIT_TET, 1.0e4, 0.3);
            let mut pos = scaled_tet(scale);
            assert!((det3x3(c.deformation_gradient(&pos)) - start_det).abs() < 1.0e-12);
            for _ in 0..100 {
                c.reset_lambda();
                c.project_hydrostatic_only(&mut pos, &inv_mass, 1.0e-2);
            }
            let det = det3x3(c.deformation_gradient(&pos));
            assert!(
                (det - 1.0).abs() < 1.0e-9,
                "scale {scale}: the volume constraint must settle at det F = 1, got {det}"
            );
        }
    }

    #[test]
    fn combined_projection_converges_towards_the_rest_state() {
        // F = 1.2 I violates both constraints. At F = I the gradients of
        // C_dev = ||F|| - sqrt(3) and C_hyd = det F - 1 are both proportional to
        // I, i.e. the two constraint surfaces are tangent at the rest state, and
        // by the AM-GM inequality ||F||_F^2 >= 3 det(F)^(2/3) they only meet at
        // the rotations F^T F = I. Alternating projections onto tangent surfaces
        // converge sublinearly (the residual falls roughly like 1/k, with no
        // geometric rate), so a few hundred iterations leave a residual of a few
        // 1e-3 and the tolerances below must not be tightened. What they pin
        // down is the fixed point: the old target alpha = 1.5 stopped at
        // det F = 1.5.
        let inv_mass = [1.0, 1.0, 1.0, 1.0];
        let mut c = StableNeoHookeanXpbdConstraint::new([0, 1, 2, 3], &UNIT_TET, 1.0e4, 0.3);
        let mut pos = scaled_tet(1.2);
        assert!((det3x3(c.deformation_gradient(&pos)) - 1.728).abs() < 1.0e-12);
        for _ in 0..500 {
            c.reset_lambda();
            c.project(&mut pos, &inv_mass, 1.0e-2);
        }
        let f = c.deformation_gradient(&pos);
        let det = det3x3(f);
        assert!(
            (det - 1.0).abs() < 5.0e-3 && det < 1.01,
            "the rest volume is det F = 1, got {det}"
        );
        let i_c: f64 = f.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
        assert!(
            (i_c.sqrt() - 3.0f64.sqrt()).abs() < 1.0e-2,
            "C_dev residual {}",
            i_c.sqrt() - 3.0f64.sqrt()
        );
    }
}
