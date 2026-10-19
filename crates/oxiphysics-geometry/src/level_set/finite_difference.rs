// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Finite-difference derivatives and mean curvature of a [`LevelSetField`].
//!
//! Every derivative is built from per-axis stencils that never leave the grid:
//! the central difference in the interior and a one-sided difference on the grid
//! boundary (the same convention as [`LevelSetField::gradient`]).  Mixed second
//! derivatives are the composition of two first-derivative stencils, so
//! `∂²φ/∂x∂y` is the classic four-corner `(φ₊₊ − φ₊₋ − φ₋₊ + φ₋₋) / (4Δ²)` in the
//! interior and degrades consistently on faces, edges and corners.

use super::types::LevelSetField;

/// Smallest `|∇φ|³` for which the level-set normal counts as defined.  Below it
/// (a flat region or a critical point of φ) the curvature is reported as `0`.
const MIN_GRAD_CUBED: f64 = 1.0e-30;

/// Grid extents along x, y and z.
fn extents(field: &LevelSetField) -> [usize; 3] {
    [field.nx, field.ny, field.nz]
}

/// The node `idx` with its coordinate along `axis` replaced by `value`.
fn with_index(idx: [usize; 3], axis: usize, value: usize) -> [usize; 3] {
    let mut shifted = idx;
    if let Some(slot) = shifted.get_mut(axis) {
        *slot = value;
    }
    shifted
}

/// Field value at the node `idx`.
fn value_at(field: &LevelSetField, idx: [usize; 3]) -> f64 {
    field.get(idx[0], idx[1], idx[2])
}

/// Node indices `(lower, upper)` of the first-derivative stencil at position `i`
/// of an axis with `n` nodes.
///
/// Interior positions use the central pair `(i − 1, i + 1)`, the two end
/// positions the adjacent pair; an axis with fewer than two nodes yields
/// `(i, i)`, i.e. a zero derivative.  `i` must be smaller than `n`.
pub(super) fn first_stencil(i: usize, n: usize) -> (usize, usize) {
    if n < 2 {
        (i, i)
    } else if i == 0 {
        (0, 1)
    } else if i + 1 >= n {
        (n - 2, n - 1)
    } else {
        (i - 1, i + 1)
    }
}

/// First partial derivative `∂φ/∂x_axis` at the node `idx`.
fn first_partial(field: &LevelSetField, axis: usize, idx: [usize; 3]) -> f64 {
    let n = extents(field).get(axis).copied().unwrap_or(0);
    let here = idx.get(axis).copied().unwrap_or(0);
    let (lower, upper) = first_stencil(here, n);
    if upper == lower {
        return 0.0;
    }
    let span = (upper - lower) as f64 * field.dx;
    (value_at(field, with_index(idx, axis, upper)) - value_at(field, with_index(idx, axis, lower)))
        / span
}

/// Second partial derivative `∂²φ/∂x_axis²` at the node `idx`.
///
/// The second difference of three consecutive nodes: centred on the node in the
/// interior, the first/last three nodes on the grid boundary.  An axis with fewer
/// than three nodes has no second difference and yields `0`.
fn second_partial(field: &LevelSetField, axis: usize, idx: [usize; 3]) -> f64 {
    let n = extents(field).get(axis).copied().unwrap_or(0);
    let here = idx.get(axis).copied().unwrap_or(0);
    if n < 3 {
        return 0.0;
    }
    let first = if here == 0 {
        0
    } else if here + 1 >= n {
        n - 3
    } else {
        here - 1
    };
    let f0 = value_at(field, with_index(idx, axis, first));
    let f1 = value_at(field, with_index(idx, axis, first + 1));
    let f2 = value_at(field, with_index(idx, axis, first + 2));
    (f0 - 2.0 * f1 + f2) / (field.dx * field.dx)
}

/// Mixed second partial derivative `∂²φ/∂x_a∂x_b` (`a ≠ b`) at the node `idx`:
/// the `b`-derivative of the `a`-derivative.
fn mixed_partial(field: &LevelSetField, a: usize, b: usize, idx: [usize; 3]) -> f64 {
    let n = extents(field).get(b).copied().unwrap_or(0);
    let here = idx.get(b).copied().unwrap_or(0);
    let (lower, upper) = first_stencil(here, n);
    if upper == lower {
        return 0.0;
    }
    let span = (upper - lower) as f64 * field.dx;
    let high = first_partial(field, a, with_index(idx, b, upper));
    let low = first_partial(field, a, with_index(idx, b, lower));
    (high - low) / span
}

/// Gradient `∇φ` at the node `idx`: central differences in the interior and
/// one-sided differences on the grid boundary.  An axis with a single node
/// contributes `0`.
pub(super) fn gradient_at(field: &LevelSetField, idx: [usize; 3]) -> [f64; 3] {
    [
        first_partial(field, 0, idx),
        first_partial(field, 1, idx),
        first_partial(field, 2, idx),
    ]
}

/// Mean curvature `κ = ∇·(∇φ / |∇φ|)` at the grid node `(ix, iy, iz)`.
///
/// Evaluated in the expanded form
///
/// ```text
///        φxx (φy² + φz²) + φyy (φx² + φz²) + φzz (φx² + φy²)
///           − 2 (φx φy φxy + φx φz φxz + φy φz φyz)
///   κ = ───────────────────────────────────────────────────────
///                              |∇φ|³
/// ```
///
/// with all six second derivatives taken from the grid.  The curvature is positive
/// where φ increases away from a convex region (a sphere signed-distance field
/// `|x| − r` has `κ = 2 / |x|`).  Returns `0` for a node outside the grid and
/// where `∇φ` vanishes (the normal is undefined there).
pub(super) fn mean_curvature_at(field: &LevelSetField, ix: usize, iy: usize, iz: usize) -> f64 {
    if ix >= field.nx || iy >= field.ny || iz >= field.nz {
        return 0.0;
    }
    let idx = [ix, iy, iz];
    let [px, py, pz] = gradient_at(field, idx);
    let pxx = second_partial(field, 0, idx);
    let pyy = second_partial(field, 1, idx);
    let pzz = second_partial(field, 2, idx);
    let pxy = mixed_partial(field, 0, 1, idx);
    let pxz = mixed_partial(field, 0, 2, idx);
    let pyz = mixed_partial(field, 1, 2, idx);

    let grad_sq = px * px + py * py + pz * pz;
    let grad_cubed = grad_sq * grad_sq.sqrt();
    let normal_defined = grad_cubed.is_finite() && grad_cubed > MIN_GRAD_CUBED;
    if !normal_defined {
        return 0.0;
    }
    let numerator =
        pxx * (py * py + pz * pz) + pyy * (px * px + pz * pz) + pzz * (px * px + py * py)
            - 2.0 * (px * py * pxy + px * pz * pxz + py * pz * pyz);
    numerator / grad_cubed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Field `φ = f(x, y, z)` sampled on an `n³` grid with spacing `dx` whose
    /// lowest node sits at `origin`.
    fn sampled(
        n: usize,
        dx: f64,
        origin: [f64; 3],
        f: impl Fn(f64, f64, f64) -> f64,
    ) -> LevelSetField {
        let mut field = LevelSetField::new(n, n, n, dx, origin);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let p = field.node_pos(ix, iy, iz);
                    field.set(ix, iy, iz, f(p[0], p[1], p[2]));
                }
            }
        }
        field
    }

    /// Radius of the grid node `(ix, iy, iz)` of a field centred on the origin.
    fn node_radius(field: &LevelSetField, ix: usize, iy: usize, iz: usize) -> f64 {
        let p = field.node_pos(ix, iy, iz);
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    }

    #[test]
    fn stencils_stay_inside_the_grid() {
        assert_eq!(first_stencil(0, 5), (0, 1));
        assert_eq!(first_stencil(2, 5), (1, 3));
        assert_eq!(first_stencil(4, 5), (3, 4));
        assert_eq!(first_stencil(0, 2), (0, 1));
        assert_eq!(first_stencil(1, 2), (0, 1));
        assert_eq!(first_stencil(0, 1), (0, 0));
        assert_eq!(first_stencil(0, 0), (0, 0));
    }

    #[test]
    fn sphere_curvature_is_isotropic_on_axis_and_off_axis() {
        // Unit-sphere signed distance field on a grid whose nodes include
        // (1, 0, 0): the exact curvature of the SDF is 2 / |x| at every node.
        let dx = 0.05;
        let n = 61;
        let origin = [-1.5, -1.5, -1.5];
        let mut field = LevelSetField::new(n, n, n, dx, origin);
        field.init_sphere([0.0, 0.0, 0.0], 1.0);
        let to_node = |c: f64| ((c - origin[0]) / dx).round() as usize;
        let s2 = std::f64::consts::FRAC_1_SQRT_2;
        let s3 = 1.0 / 3.0_f64.sqrt();
        // (1,0,0), (1,1,0)/sqrt(2), (1,1,1)/sqrt(3) and a generic direction.
        let targets = [
            [1.0, 0.0, 0.0],
            [s2, s2, 0.0],
            [0.0, s2, s2],
            [s3, s3, s3],
            [0.8, 0.0, 0.6],
        ];
        for t in targets {
            let (ix, iy, iz) = (to_node(t[0]), to_node(t[1]), to_node(t[2]));
            let r = node_radius(&field, ix, iy, iz);
            let kappa = field.mean_curvature(ix, iy, iz);
            let exact = 2.0 / r;
            assert!(
                (kappa - exact).abs() <= 0.02 * exact,
                "node ({ix},{iy},{iz}) r={r}: kappa {kappa} vs exact {exact}"
            );
            // The curvature of the unit sphere itself, up to the distance of the
            // node from the surface (at most 4% for the generic nodes here).
            assert!(
                (kappa - 2.0).abs() <= 0.05 * 2.0,
                "node ({ix},{iy},{iz}): kappa {kappa} not within 5% of 2"
            );
        }
        // The on-axis node is exactly on the surface.
        let on_axis = field.mean_curvature(to_node(1.0), to_node(0.0), to_node(0.0));
        assert!((on_axis - 2.0).abs() < 0.02, "on-axis curvature {on_axis}");
    }

    #[test]
    fn sphere_shell_curvature_is_isotropic_on_the_audit_grid() {
        // init_sphere(R = 0.5) on a 41^3 grid over [-1, 1]^3: on every node within
        // one cell of the surface, along the axes, the face diagonals and the cube
        // diagonals alike, the curvature is the exact 2 / r of the distance field.
        let (n, dx) = (41, 0.05);
        let mut field = LevelSetField::new(n, n, n, dx, [-1.0, -1.0, -1.0]);
        field.init_sphere([0.0, 0.0, 0.0], 0.5);
        let mut shell_nodes = 0usize;
        let mut off_axis_nodes = 0usize;
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    if field.get(ix, iy, iz).abs() >= dx {
                        continue;
                    }
                    let r = node_radius(&field, ix, iy, iz);
                    let kappa = field.mean_curvature(ix, iy, iz);
                    let exact = 2.0 / r;
                    assert!(
                        (kappa - exact).abs() <= 0.05 * exact,
                        "node ({ix},{iy},{iz}) r={r}: kappa {kappa} vs exact {exact}"
                    );
                    shell_nodes += 1;
                    let off_axis = [ix, iy, iz].iter().filter(|&&i| i != 20).count() >= 2;
                    off_axis_nodes += usize::from(off_axis);
                }
            }
        }
        assert!(shell_nodes > 2000, "shell holds only {shell_nodes} nodes");
        assert!(
            off_axis_nodes > 1900,
            "only {off_axis_nodes} off-axis nodes"
        );
    }

    #[test]
    fn cylinder_curvature_ignores_the_axial_direction() {
        // φ = sqrt(x² + y²) − 0.6: κ = 1 / ρ, independent of z.
        let dx = 0.05;
        let n = 41;
        let field = sampled(n, dx, [-1.0, -1.0, -1.0], |x, y, _| {
            (x * x + y * y).sqrt() - 0.6
        });
        for (ix, iy, iz) in [(30, 20, 5), (28, 28, 20), (12, 31, 33)] {
            let p = field.node_pos(ix, iy, iz);
            let rho = (p[0] * p[0] + p[1] * p[1]).sqrt();
            let kappa = field.mean_curvature(ix, iy, iz);
            assert!(
                (kappa - 1.0 / rho).abs() <= 0.02 / rho,
                "rho={rho}: kappa {kappa} vs {}",
                1.0 / rho
            );
        }
    }

    #[test]
    fn mixed_second_derivatives_enter_the_curvature() {
        // Quadratic fields are differentiated exactly by the interior stencils, so
        // the curvature of φ = x y, x z and y z must match the closed form
        // κ = −2 a b / (a² + b²)^{3/2} (a, b the two coordinates) to rounding.
        let n = 9;
        let dx = 0.25;
        let origin = [-1.0, -1.0, -1.0];
        type FieldFn = fn(f64, f64, f64) -> f64;
        let cases: [(&str, FieldFn, FieldFn); 3] = [
            (
                "xy",
                |x, y, _| x * y,
                |x, y, _| -2.0 * x * y / (x * x + y * y).powf(1.5),
            ),
            (
                "xz",
                |x, _, z| x * z,
                |x, _, z| -2.0 * x * z / (x * x + z * z).powf(1.5),
            ),
            (
                "yz",
                |_, y, z| y * z,
                |_, y, z| -2.0 * y * z / (y * y + z * z).powf(1.5),
            ),
        ];
        for (name, phi, kappa_exact) in cases {
            let field = sampled(n, dx, origin, phi);
            for (ix, iy, iz) in [(6, 3, 2), (2, 6, 5), (5, 5, 7), (3, 6, 6)] {
                let p = field.node_pos(ix, iy, iz);
                let exact = kappa_exact(p[0], p[1], p[2]);
                let kappa = field.mean_curvature(ix, iy, iz);
                assert!(
                    (kappa - exact).abs() < 1e-9,
                    "phi={name} at ({ix},{iy},{iz}): {kappa} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn general_quadratic_field_matches_the_closed_form() {
        // φ = xy + yz + zx + x + 2y − z:  ∇φ = (y+z+1, x+z+2, x+y−1), constant
        // Hessian with unit off-diagonal entries, so
        // κ = −2 (φx φy + φx φz + φy φz) / |∇φ|³.
        let field = sampled(11, 0.2, [-1.0, -1.0, -1.0], |x, y, z| {
            x * y + y * z + z * x + x + 2.0 * y - z
        });
        for (ix, iy, iz) in [(5, 5, 5), (3, 7, 4), (8, 2, 6), (1, 1, 9), (9, 9, 1)] {
            let p = field.node_pos(ix, iy, iz);
            let (gx, gy, gz) = (p[1] + p[2] + 1.0, p[0] + p[2] + 2.0, p[0] + p[1] - 1.0);
            let exact =
                -2.0 * (gx * gy + gx * gz + gy * gz) / (gx * gx + gy * gy + gz * gz).powf(1.5);
            let kappa = field.mean_curvature(ix, iy, iz);
            assert!(
                (kappa - exact).abs() < 1e-9,
                "({ix},{iy},{iz}): {kappa} vs {exact}"
            );
        }
    }

    #[test]
    fn planes_have_zero_curvature_even_on_corners_and_faces() {
        let field = sampled(6, 0.3, [0.0, 0.0, 0.0], |x, y, z| x + 2.0 * y - z);
        for iz in 0..6 {
            for iy in 0..6 {
                for ix in 0..6 {
                    let kappa = field.mean_curvature(ix, iy, iz);
                    assert!(
                        kappa.abs() < 1e-10,
                        "plane curvature at ({ix},{iy},{iz}) is {kappa}"
                    );
                }
            }
        }
    }

    #[test]
    fn boundary_nodes_of_a_sphere_field_stay_finite_and_positive() {
        // The grid boundary lies outside the sphere: one-sided stencils must still
        // give a finite, positive curvature of roughly the right size.
        let mut field = LevelSetField::new(21, 21, 21, 0.1, [-1.0, -1.0, -1.0]);
        field.init_sphere([0.0, 0.0, 0.0], 0.5);
        for (ix, iy, iz) in [
            (0, 10, 10),
            (20, 10, 10),
            (0, 0, 0),
            (20, 20, 20),
            (10, 0, 20),
        ] {
            let kappa = field.mean_curvature(ix, iy, iz);
            let r = node_radius(&field, ix, iy, iz);
            assert!(kappa.is_finite() && kappa > 0.0, "kappa {kappa} at r={r}");
            assert!(
                (kappa - 2.0 / r).abs() < 0.5 * 2.0 / r,
                "boundary kappa {kappa} far from {}",
                2.0 / r
            );
        }
    }

    #[test]
    fn degenerate_inputs_return_zero_without_panicking() {
        // Out-of-range node, empty grid, single-node grid and a constant field.
        let field = LevelSetField::new(4, 4, 4, 0.1, [0.0; 3]);
        assert_eq!(field.mean_curvature(4, 0, 0), 0.0);
        assert_eq!(field.mean_curvature(0, 9, 0), 0.0);
        let empty = LevelSetField::new(0, 0, 0, 0.1, [0.0; 3]);
        assert_eq!(empty.mean_curvature(0, 0, 0), 0.0);
        let mut single = LevelSetField::new(1, 1, 1, 0.1, [0.0; 3]);
        single.set(0, 0, 0, 3.0);
        assert_eq!(single.mean_curvature(0, 0, 0), 0.0);
        let mut flat = LevelSetField::new(5, 5, 5, 0.1, [0.0; 3]);
        flat.data.iter_mut().for_each(|v| *v = 2.5);
        assert_eq!(flat.mean_curvature(2, 2, 2), 0.0);
        // A two-dimensional slab (single layer along z) is a valid field.
        let mut slab = LevelSetField::new(21, 21, 1, 0.1, [-1.0, -1.0, 0.0]);
        slab.init_sphere([0.0, 0.0, 0.0], 0.5);
        let kappa = slab.mean_curvature(15, 10, 0);
        // A circle of radius 0.5 in the slab: κ = 1 / ρ.
        let rho = 0.5;
        assert!(
            (kappa - 1.0 / rho).abs() < 0.1 * (1.0 / rho),
            "slab curvature {kappa}"
        );
    }

    #[test]
    fn mean_curvature_flow_shrinks_a_sphere_isotropically() {
        use crate::level_set::LevelSetEvolution;
        // ∂φ/∂t = κ |∇φ| moves every node of a sphere SDF by 2 / r per unit time.
        let dx = 0.1;
        let n = 31;
        let origin = [-1.5, -1.5, -1.5];
        let mut field = LevelSetField::new(n, n, n, dx, origin);
        field.init_sphere([0.0, 0.0, 0.0], 1.0);
        let before = field.data.clone();
        let evolution = LevelSetEvolution::new(0.002, 3);
        evolution.mean_curvature_flow(&mut field);
        let node = |c: f64| ((c - origin[0]) / dx).round() as usize;
        let elapsed = 0.002 * 3.0;
        let on_axis = (node(1.0), node(0.0), node(0.0));
        let diagonal = (node(0.7), node(0.7), node(0.0));
        let corner = (node(0.6), node(0.6), node(0.6));
        for (ix, iy, iz) in [on_axis, diagonal, corner] {
            let r = node_radius(&field, ix, iy, iz);
            let moved = field.get(ix, iy, iz) - before[field.idx(ix, iy, iz)];
            let expected = 2.0 / r * elapsed;
            assert!(
                (moved - expected).abs() <= 0.05 * expected,
                "node ({ix},{iy},{iz}) r={r}: moved {moved}, expected {expected}"
            );
        }
        assert!(field.data.iter().all(|v| v.is_finite()));
    }
}
