// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`SdfReinitialize`], the Eikonal redistancing of an
//! [`SdfGrid3D`].
//!
//! The old implementation took a single one-dimensional minimum over the six
//! neighbours, so an inflated field shrank by one cell and a compressed field was
//! never repaired.  The tests check the defining properties of the true Eikonal
//! solution: `|∇φ| = 1`, the correct distance to the zero level set, unchanged signs,
//! and exact reproduction of planar fronts.

use super::{SdfGrid3D, SdfReinitialize};

/// A grid of `n^3` samples spanning `[-1, 1]^3` (the samples sit at the voxel
/// centres, so the grid origin lies half a voxel before the first sample).
fn unit_cube_grid<F: Fn([f64; 3]) -> f64>(n: usize, f: F) -> SdfGrid3D {
    let h = 2.0 / (n - 1) as f64;
    SdfGrid3D::from_sdf(n, n, n, [-1.0 - 0.5 * h; 3], h, f)
}

fn norm(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

/// `scale * (|x| - radius)` sampled on `n^3` nodes of `[-1, 1]^3`.
fn scaled_sphere(n: usize, radius: f64, scale: f64) -> SdfGrid3D {
    unit_cube_grid(n, |p| scale * (norm(p) - radius))
}

/// Magnitude of the central-difference gradient at an interior node.
fn gradient_magnitude(grid: &SdfGrid3D, ix: usize, iy: usize, iz: usize) -> f64 {
    let span = 2.0 * grid.voxel_size;
    let g = [
        (grid.get(ix + 1, iy, iz) - grid.get(ix - 1, iy, iz)) / span,
        (grid.get(ix, iy + 1, iz) - grid.get(ix, iy - 1, iz)) / span,
        (grid.get(ix, iy, iz + 1) - grid.get(ix, iy, iz - 1)) / span,
    ];
    norm(g)
}

/// Interior nodes of the central box `|x|, |y|, |z| <= 0.8` that are at least
/// `interface_cells` cells from the sphere surface and `centre_cells` cells from
/// the sphere centre.
///
/// The first-order upwind scheme resolves a circular front with a gradient error of
/// about `0.5 h / r`, and the solution has a conical singularity at the centre, so
/// the unit-gradient checks leave out the neighbourhood of the centre.
fn band_nodes(
    grid: &SdfGrid3D,
    radius: f64,
    interface_cells: f64,
    centre_cells: f64,
) -> Vec<(usize, usize, usize)> {
    let h = grid.voxel_size;
    let mut nodes = Vec::new();
    for iz in 1..grid.nz - 1 {
        for iy in 1..grid.ny - 1 {
            for ix in 1..grid.nx - 1 {
                let p = grid.voxel_center(ix, iy, iz);
                let r = norm(p);
                let inside_box = p.iter().all(|c| c.abs() <= 0.8);
                if inside_box && r >= centre_cells * h && (r - radius).abs() >= interface_cells * h
                {
                    nodes.push((ix, iy, iz));
                }
            }
        }
    }
    nodes
}

/// Redistance `scale * (|x| - radius)` and check unit gradient, signed distance and
/// unchanged signs.
fn assert_redistanced_sphere(n: usize, radius: f64, scale: f64) {
    let mut grid = scaled_sphere(n, radius, scale);
    let before = grid.data.clone();
    SdfReinitialize::reinitialize(&mut grid);

    let nodes = band_nodes(&grid, radius, 2.0, 14.0);
    assert!(nodes.len() > 20_000, "band too small: {}", nodes.len());
    let h = grid.voxel_size;
    for &(ix, iy, iz) in &nodes {
        let g = gradient_magnitude(&grid, ix, iy, iz);
        let r = norm(grid.voxel_center(ix, iy, iz));
        assert!(
            (g - 1.0).abs() <= 0.05,
            "scale {scale}: |grad phi| = {g} at ({ix},{iy},{iz}), r = {r}"
        );
        let value = grid.get(ix, iy, iz);
        assert!(
            (value - (r - radius)).abs() <= h,
            "scale {scale}: phi = {value} but the signed distance is {}",
            r - radius
        );
    }
    for (i, (&old, &new)) in before.iter().zip(grid.data.iter()).enumerate() {
        assert_eq!(
            old < 0.0,
            new < 0.0,
            "sign of node {i} changed: {old} -> {new}"
        );
    }
}

#[test]
fn inflated_sphere_field_is_redistanced_to_unit_gradient() {
    // Three times too steep: the old one-cell-per-pass update left it at |grad| = 3.
    assert_redistanced_sphere(65, 0.6, 3.0);
}

#[test]
fn compressed_sphere_field_is_redistanced_to_unit_gradient() {
    // Four times too flat: the old update never lowered any value of such a field.
    assert_redistanced_sphere(65, 0.6, 0.25);
}

#[test]
fn a_true_signed_distance_field_stays_a_signed_distance_field() {
    assert_redistanced_sphere(65, 0.6, 1.0);
}

#[test]
fn planar_fronts_are_reproduced_exactly() {
    // A plane 0.3 h above the node plane k = 6 of every axis, with a wrong slope of 7.
    let n = 17;
    let h = 2.0 / (n - 1) as f64;
    let shift = -1.0 + 6.3 * h;
    for axis in 0..3 {
        let mut grid = unit_cube_grid(n, |p| 7.0 * (p[axis] - shift));
        SdfReinitialize::reinitialize(&mut grid);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let p = grid.voxel_center(ix, iy, iz);
                    let exact = p[axis] - shift;
                    let got = grid.get(ix, iy, iz);
                    assert!(
                        (got - exact).abs() <= 1e-12,
                        "axis {axis} node ({ix},{iy},{iz}): {got} vs {exact}"
                    );
                }
            }
        }
    }
}

#[test]
fn tilted_planar_front_is_reproduced_in_the_interior() {
    // The Godunov update is exact for a planar front in any direction and the
    // frozen layer carries the exact distance of a linear field.  The one-sided
    // stencils at the grid faces are not exact for an oblique front, and their
    // error decays by about a factor of three per cell into the grid, so the
    // exact solution is only recovered deep inside.
    let n = 33;
    let margin = 10;
    let normal = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let offset = 0.07;
    let plane = |p: [f64; 3]| -(normal[0] * p[0] + normal[1] * p[1] + normal[2] * p[2] - offset);
    let mut grid = unit_cube_grid(n, |p| 5.0 * plane(p));
    SdfReinitialize::reinitialize(&mut grid);
    let mut checked = 0usize;
    for iz in margin..n - margin {
        for iy in margin..n - margin {
            for ix in margin..n - margin {
                let exact = plane(grid.voxel_center(ix, iy, iz));
                let got = grid.get(ix, iy, iz);
                assert!(
                    (got - exact).abs() <= 5e-6,
                    "node ({ix},{iy},{iz}): {got} vs {exact}"
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 13 * 13 * 13);
}

#[test]
fn two_disjoint_spheres_give_the_distance_to_the_nearer_one() {
    // The fronts of both spheres meet halfway: the field needs several sweep
    // directions and the minimum of two distances.
    let (c1, r1) = ([-0.45, 0.0, 0.05], 0.25);
    let (c2, r2) = ([0.5, 0.1, -0.05], 0.3);
    let exact = |p: [f64; 3]| {
        let d1 = norm([p[0] - c1[0], p[1] - c1[1], p[2] - c1[2]]) - r1;
        let d2 = norm([p[0] - c2[0], p[1] - c2[1], p[2] - c2[2]]) - r2;
        d1.min(d2)
    };
    let n = 65;
    let mut grid = unit_cube_grid(n, |p| 4.0 * exact(p));
    let before = grid.data.clone();
    SdfReinitialize::reinitialize(&mut grid);
    let h = grid.voxel_size;
    let mut checked = 0usize;
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let p = grid.voxel_center(ix, iy, iz);
                let e = exact(p);
                // Outside, everywhere; inside, only the shell near the surface (the
                // centre of a sphere is a conical singularity).
                if e > 2.0 * h || (-0.1..-2.0 * h).contains(&e) {
                    let got = grid.get(ix, iy, iz);
                    assert!(
                        (got - e).abs() <= 1.5 * h,
                        "node ({ix},{iy},{iz}): {got} vs exact {e}"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 100_000, "only {checked} nodes checked");
    for (old, new) in before.iter().zip(grid.data.iter()) {
        assert_eq!(*old < 0.0, *new < 0.0);
    }
}

#[test]
fn a_spike_of_the_right_sign_is_replaced_by_the_distance() {
    let n = 33;
    let radius = 0.5;
    let mut grid = scaled_sphere(n, radius, 1.0);
    let h = grid.voxel_size;
    // Far outside the sphere, in a corner, and deep inside, off the exact centre.
    let outside = (30, 29, 30);
    let inside = (17, 16, 17);
    grid.set(outside.0, outside.1, outside.2, 1.0e3);
    grid.set(inside.0, inside.1, inside.2, -1.0e3);
    SdfReinitialize::reinitialize(&mut grid);
    let value_outside = grid.get(outside.0, outside.1, outside.2);
    let exact_outside = norm(grid.voxel_center(outside.0, outside.1, outside.2)) - radius;
    assert!(
        (value_outside - exact_outside).abs() <= h,
        "outside spike: {value_outside} vs {exact_outside}"
    );
    let value_inside = grid.get(inside.0, inside.1, inside.2);
    let exact_inside = norm(grid.voxel_center(inside.0, inside.1, inside.2)) - radius;
    assert!(
        (value_inside - exact_inside).abs() <= h,
        "inside spike: {value_inside} vs {exact_inside}"
    );
}

#[test]
fn a_sign_flipped_spike_becomes_a_small_bubble() {
    // An "outside" sample inside the sphere keeps its sign: it is a tiny island, so
    // its distance to the zero level set is below one voxel, not 100.
    let mut grid = scaled_sphere(33, 0.5, 1.0);
    let h = grid.voxel_size;
    grid.set(16, 16, 16, 100.0);
    SdfReinitialize::reinitialize(&mut grid);
    let v = grid.get(16, 16, 16);
    assert!(v > 0.0 && v <= 1.5 * h, "bubble value {v} (h = {h})");
    // Its neighbours (negative) are within a voxel of the bubble.
    let w = grid.get(17, 16, 16);
    assert!(w < 0.0 && w >= -h, "neighbour value {w}");
}

#[test]
fn narrow_band_grids_are_filled_from_the_band() {
    // Outside the band |phi| <= 0.3 the samples are +-infinity, as in a narrow-band
    // grid.  Redistancing recomputes them from the band.
    let n = 49;
    let radius = 0.55;
    let mut grid = unit_cube_grid(n, |p| {
        let d = norm(p) - radius;
        if d.abs() <= 0.3 {
            d
        } else {
            d.signum() * f64::INFINITY
        }
    });
    let infinite_before = grid.data.iter().filter(|v| v.is_infinite()).count();
    assert!(infinite_before > grid.data.len() / 2);
    SdfReinitialize::reinitialize(&mut grid);
    let h = grid.voxel_size;
    assert!(
        grid.data.iter().all(|v| v.is_finite()),
        "every sample must be reached"
    );
    let mut checked = 0usize;
    for iz in 1..n - 1 {
        for iy in 1..n - 1 {
            for ix in 1..n - 1 {
                let p = grid.voxel_center(ix, iy, iz);
                let r = norm(p);
                let exact = r - radius;
                // Sign is preserved.
                assert_eq!(grid.get(ix, iy, iz) < 0.0, exact < 0.0, "sign at {p:?}");
                // Outside the sphere the distance field is exact up to the scheme.
                if exact > 2.0 * h && p.iter().all(|c| c.abs() <= 0.9) {
                    assert!(
                        (grid.get(ix, iy, iz) - exact).abs() <= 1.5 * h,
                        "node ({ix},{iy},{iz}): {} vs {exact}",
                        grid.get(ix, iy, iz)
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 5000, "only {checked} nodes checked");
}

#[test]
fn a_surface_known_only_from_the_inside_is_placed_on_the_boundary_samples() {
    // The inside carries finite distances and the outside is the `+inf` a fresh grid
    // is filled with, so no sample says where between an inside and an outside
    // sample the surface lies.  Those samples end up on the zero level set and the
    // outside distances are measured from them: finite, with the signs kept and
    // less than a voxel short of the exact ones (the outer layer of samples lies up
    // to one voxel beyond the true surface).
    let n = 33;
    let radius = 0.5;
    let mut grid = unit_cube_grid(n, |p| {
        let d = norm(p) - radius;
        if d < 0.0 { d } else { f64::INFINITY }
    });
    let outside_before = grid.data.iter().filter(|v| v.is_infinite()).count();
    assert!(outside_before > grid.data.len() / 2);
    SdfReinitialize::reinitialize(&mut grid);
    let h = grid.voxel_size;
    assert!(grid.data.iter().all(|v| v.is_finite()), "NaN or infinity");
    let mut checked = 0usize;
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let exact = norm(grid.voxel_center(ix, iy, iz)) - radius;
                let value = grid.get(ix, iy, iz);
                if exact < 0.0 {
                    assert!(value <= 0.0, "inside node ({ix},{iy},{iz}) became {value}");
                } else {
                    assert!(value >= 0.0, "outside node ({ix},{iy},{iz}) is {value}");
                }
                if exact >= 2.0 * h {
                    assert!(
                        (value - exact).abs() <= h,
                        "node ({ix},{iy},{iz}): {value} vs exact {exact}"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 10_000, "only {checked} nodes checked");
}

#[test]
fn a_grid_without_a_zero_level_set_is_unchanged() {
    let mut grid = unit_cube_grid(9, |p| 2.0 + p[0] * p[0]);
    let before = grid.data.clone();
    SdfReinitialize::reinitialize(&mut grid);
    assert_eq!(grid.data, before);
    let mut infinite = SdfGrid3D::new(4, 4, 4, [0.0; 3], 0.5);
    SdfReinitialize::reinitialize(&mut infinite);
    assert!(infinite.data.iter().all(|v| *v == f64::INFINITY));
}

#[test]
fn invalid_grids_are_left_unchanged() {
    // A NaN sample.
    let mut grid = scaled_sphere(9, 0.5, 2.0);
    grid.set(4, 4, 4, f64::NAN);
    let before = grid.data.clone();
    SdfReinitialize::reinitialize(&mut grid);
    for (a, b) in before.iter().zip(grid.data.iter()) {
        assert!(a.to_bits() == b.to_bits(), "{a} -> {b}");
    }
    // A degenerate voxel size.
    for voxel_size in [0.0, -1.0, f64::NAN] {
        let mut grid = scaled_sphere(9, 0.5, 2.0);
        grid.voxel_size = voxel_size;
        let before = grid.data.clone();
        SdfReinitialize::reinitialize(&mut grid);
        assert_eq!(grid.data, before, "voxel size {voxel_size}");
    }
    // An empty grid.
    let mut empty = SdfGrid3D::new(0, 4, 4, [0.0; 3], 1.0);
    SdfReinitialize::reinitialize(&mut empty);
    assert!(empty.data.is_empty());
}

#[test]
fn a_single_layer_grid_is_redistanced_in_two_dimensions() {
    // nz = 1: the circle |x| = 0.5 in the plane, field scaled by 3.
    let (nx, ny) = (81, 81);
    let h = 2.0 / (nx - 1) as f64;
    let mut grid = SdfGrid3D::from_sdf(nx, ny, 1, [-1.0 - 0.5 * h, -1.0 - 0.5 * h, 0.0], h, |p| {
        3.0 * (p[0].hypot(p[1]) - 0.5)
    });
    SdfReinitialize::reinitialize(&mut grid);
    let mut checked = 0usize;
    for iy in 1..ny - 1 {
        for ix in 1..nx - 1 {
            let p = grid.voxel_center(ix, iy, 0);
            let r = p[0].hypot(p[1]);
            if r >= 8.0 * h && (r - 0.5).abs() >= 2.0 * h && r <= 0.9 {
                let gx = (grid.get(ix + 1, iy, 0) - grid.get(ix - 1, iy, 0)) / (2.0 * h);
                let gy = (grid.get(ix, iy + 1, 0) - grid.get(ix, iy - 1, 0)) / (2.0 * h);
                let magnitude = gx.hypot(gy);
                assert!((magnitude - 1.0).abs() <= 0.05, "|grad| = {magnitude}");
                assert!(
                    (grid.get(ix, iy, 0) - (r - 0.5)).abs() <= h,
                    "phi = {} at r = {r}",
                    grid.get(ix, iy, 0)
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 1500, "only {checked} nodes checked");
}
