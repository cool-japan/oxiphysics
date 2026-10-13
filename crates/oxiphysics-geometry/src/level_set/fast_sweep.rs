// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fast-sweeping redistancing of a [`LevelSetField`] (Zhao, *A fast sweeping method
//! for Eikonal equations*, Math. Comp. 74, 2005).
//!
//! The method solves `|∇d| = 1` for the unsigned distance `d` to the zero level
//! set with Gauss–Seidel sweeps in the eight diagonal grid orderings and the
//! upwind Godunov discretisation.  Nodes next to the interface are *frozen* with a
//! first-order estimate of their distance; every other node starts at a huge value
//! and only ever decreases, so over-compressed and under-compressed fields are
//! both repaired and the sign (inside/outside) of every node is preserved.

use super::finite_difference::gradient_at;
use super::types::{LevelSetField, LevelSetReinit};

/// Distance value of a node that no sweep has reached yet.
const FAR: f64 = 1.0e10;

/// Number of diagonal orderings of the Gauss–Seidel sweeps.
const ORDERINGS: usize = 8;

/// Upwind Godunov solution `u` of `Σ ((u − a_k)⁺)² = h²` for the neighbour
/// minima `a ≤ b ≤ c` (`f64::INFINITY` marks a missing neighbour).
fn godunov_update(sorted: [f64; 3], h: f64) -> f64 {
    let [a, b, c] = sorted;
    let mut u = a + h;
    if u > b {
        let discriminant = 2.0 * h * h - (a - b) * (a - b);
        u = 0.5 * (a + b + discriminant.max(0.0).sqrt());
        if u > c {
            // Larger root of 3u² − 2(a+b+c)u + a² + b² + c² − h² = 0.
            let sum = a + b + c;
            let rest = a * a + b * b + c * c - h * h;
            u = (sum + (sum * sum - 3.0 * rest).max(0.0).sqrt()) / 3.0;
        }
    }
    u
}

/// Smaller known distance among the (up to two) neighbours of node `idx` along
/// one axis; `f64::INFINITY` when there is none.  `stride` is the index distance
/// between neighbours, `pos` the node's coordinate along the axis and `n` the
/// extent of the axis.
fn axis_minimum(dist: &[f64], idx: usize, stride: usize, pos: usize, n: usize) -> f64 {
    let mut best = f64::INFINITY;
    let mut consider = |neighbour: usize| {
        if let Some(&d) = dist.get(neighbour)
            && d < FAR
        {
            best = best.min(d);
        }
    };
    if pos > 0 {
        consider(idx - stride);
    }
    if pos + 1 < n {
        consider(idx + stride);
    }
    best
}

/// Signed-distance estimate `|φ| / |∇φ|` of an interface-adjacent node, limited by
/// the distance to the closest zero crossing along a grid edge.
///
/// `crossing` is that smallest edge-crossing distance (linear interpolation of φ
/// between the two nodes of every edge whose end points differ in sign).
fn interface_distance(field: &LevelSetField, node: [usize; 3], phi: f64, crossing: f64) -> f64 {
    let [gx, gy, gz] = gradient_at(field, node);
    let grad = (gx * gx + gy * gy + gz * gz).sqrt();
    let estimate = phi.abs() / grad;
    if estimate.is_finite() {
        estimate.min(crossing)
    } else {
        crossing
    }
}

/// Redistance `field` in place with at most `n_sweeps` Gauss–Seidel sweeps.
///
/// See [`LevelSetReinit::fast_sweep`].
pub(super) fn fast_sweep(reinit: &LevelSetReinit, field: &mut LevelSetField, n_sweeps: usize) {
    let (nx, ny, nz) = (field.nx, field.ny, field.nz);
    let h = field.dx;
    let Some(total) = nx.checked_mul(ny).and_then(|v| v.checked_mul(nz)) else {
        return;
    };
    let valid_grid = total > 0 && field.data.len() == total && h.is_finite() && h > 0.0;
    if n_sweeps == 0 || !valid_grid || field.data.iter().any(|v| v.is_nan()) {
        return;
    }
    let (stride_y, stride_z) = (nx, nx * ny);

    // Freeze the nodes next to the interface and give them their initial distance;
    // every other node starts "far away".
    let mut frozen = vec![false; total];
    let mut dist = vec![FAR; total];
    let mut any_frozen = false;
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iz * stride_z + iy * stride_y + ix;
                let phi = field.data[i];
                if phi == 0.0 {
                    frozen[i] = true;
                    dist[i] = 0.0;
                    any_frozen = true;
                    continue;
                }
                let mut crossing = f64::INFINITY;
                let neighbours = [
                    (ix > 0).then(|| i - 1),
                    (ix + 1 < nx).then(|| i + 1),
                    (iy > 0).then(|| i - stride_y),
                    (iy + 1 < ny).then(|| i + stride_y),
                    (iz > 0).then(|| i - stride_z),
                    (iz + 1 < nz).then(|| i + stride_z),
                ];
                for j in neighbours.into_iter().flatten() {
                    let other = field.data[j];
                    if other != 0.0 && (other > 0.0) != (phi > 0.0) {
                        // Linear zero crossing between the two nodes, written to
                        // stay finite for huge magnitudes.
                        let (a, b) = (phi.abs(), other.abs());
                        let scale = a.max(b);
                        let share = (a / scale) / ((a / scale) + (b / scale));
                        crossing = crossing.min(h * share);
                    }
                }
                if crossing.is_finite() {
                    frozen[i] = true;
                    dist[i] = interface_distance(field, [ix, iy, iz], phi, crossing);
                    any_frozen = true;
                }
            }
        }
    }
    if !any_frozen {
        // No zero level set: there is nothing to measure a distance from.
        return;
    }

    // Gauss–Seidel sweeps; sweep `s` uses the diagonal ordering `s mod 8`.
    let mut cycle_change = 0.0_f64;
    for sweep in 0..n_sweeps {
        let ordering = sweep % ORDERINGS;
        if ordering == 0 {
            cycle_change = 0.0;
        }
        let ascending = [ordering & 1 == 0, ordering & 2 == 0, ordering & 4 == 0];
        for kz in 0..nz {
            let iz = if ascending[2] { kz } else { nz - 1 - kz };
            for ky in 0..ny {
                let iy = if ascending[1] { ky } else { ny - 1 - ky };
                for kx in 0..nx {
                    let ix = if ascending[0] { kx } else { nx - 1 - kx };
                    let i = iz * stride_z + iy * stride_y + ix;
                    if frozen[i] {
                        continue;
                    }
                    let mut minima = [
                        axis_minimum(&dist, i, 1, ix, nx),
                        axis_minimum(&dist, i, stride_y, iy, ny),
                        axis_minimum(&dist, i, stride_z, iz, nz),
                    ];
                    minima.sort_by(f64::total_cmp);
                    if !minima[0].is_finite() {
                        continue;
                    }
                    let updated = godunov_update(minima, h);
                    if updated < dist[i] {
                        let change = if dist[i] >= FAR {
                            f64::INFINITY
                        } else {
                            dist[i] - updated
                        };
                        cycle_change = cycle_change.max(change);
                        dist[i] = updated;
                    }
                }
            }
        }
        if ordering == ORDERINGS - 1 && cycle_change <= reinit.tol {
            break;
        }
    }

    // Restore the sign of every node; nodes no sweep reached keep their value.
    for (value, &d) in field.data.iter_mut().zip(dist.iter()) {
        if d < FAR {
            let sign = if *value > 0.0 {
                1.0
            } else if *value < 0.0 {
                -1.0
            } else {
                0.0
            };
            *value = sign * d;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sphere field `φ = scale (|x − c| − r)` on an `n³` grid over `[-1, 1]³`.
    fn scaled_sphere(n: usize, radius: f64, scale: f64) -> LevelSetField {
        let dx = 2.0 / (n - 1) as f64;
        let mut field = LevelSetField::new(n, n, n, dx, [-1.0, -1.0, -1.0]);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let p = field.node_pos(ix, iy, iz);
                    let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    field.set(ix, iy, iz, scale * (r - radius));
                }
            }
        }
        field
    }

    fn radius_of(field: &LevelSetField, ix: usize, iy: usize, iz: usize) -> f64 {
        let p = field.node_pos(ix, iy, iz);
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    }

    /// Nodes well away from the interface, the sphere centre and the grid faces,
    /// where the redistanced gradient magnitude is checked.
    fn band_nodes(field: &LevelSetField, radius: f64) -> Vec<(usize, usize, usize)> {
        let mut nodes = Vec::new();
        for iz in 1..field.nz - 1 {
            for iy in 1..field.ny - 1 {
                for ix in 1..field.nx - 1 {
                    let r = radius_of(field, ix, iy, iz);
                    let p = field.node_pos(ix, iy, iz);
                    let inside_box = p.iter().all(|c| c.abs() <= 0.8);
                    if inside_box && r > 0.15 && (r - radius).abs() >= 0.1 {
                        nodes.push((ix, iy, iz));
                    }
                }
            }
        }
        nodes
    }

    fn gradient_magnitude(field: &LevelSetField, ix: usize, iy: usize, iz: usize) -> f64 {
        let g = field.gradient(ix, iy, iz);
        (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
    }

    #[test]
    fn over_compressed_sphere_field_is_redistanced_to_unit_gradient() {
        let (n, radius) = (33, 0.5);
        let mut field = scaled_sphere(n, radius, 2.0);
        let before = field.data.clone();
        let reinit = LevelSetReinit::new(1e-12, 100);
        reinit.fast_sweep(&mut field, 24);

        let nodes = band_nodes(&field, radius);
        assert!(nodes.len() > 1000, "band too small: {}", nodes.len());
        let dx = field.dx;
        for &(ix, iy, iz) in &nodes {
            let g = gradient_magnitude(&field, ix, iy, iz);
            assert!(
                (g - 1.0).abs() <= 0.05,
                "|grad phi| = {g} at ({ix},{iy},{iz}), r = {}",
                radius_of(&field, ix, iy, iz)
            );
            let exact = radius_of(&field, ix, iy, iz) - radius;
            let value = field.get(ix, iy, iz);
            assert!(
                (value - exact).abs() <= dx,
                "phi = {value} but signed distance is {exact}"
            );
        }
        // The inside/outside classification never changes.
        for (b, a) in before.iter().zip(field.data.iter()) {
            assert_eq!(*b > 0.0, *a > 0.0);
            assert_eq!(*b < 0.0, *a < 0.0);
        }
    }

    #[test]
    fn under_compressed_sphere_field_is_redistanced_to_unit_gradient() {
        let (n, radius) = (33, 0.5);
        let mut field = scaled_sphere(n, radius, 0.25);
        LevelSetReinit::new(1e-12, 100).fast_sweep(&mut field, 24);
        for &(ix, iy, iz) in &band_nodes(&field, radius) {
            let g = gradient_magnitude(&field, ix, iy, iz);
            assert!(
                (g - 1.0).abs() <= 0.05,
                "|grad phi| = {g} at ({ix},{iy},{iz})"
            );
        }
    }

    #[test]
    fn distance_propagates_along_every_axis() {
        // A plane normal to y, 0.3 above the origin, with a wildly wrong slope: the
        // y-neighbour propagation of the Godunov update must repair it.
        let n = 17;
        let dx = 0.125;
        let mut field = LevelSetField::new(n, n, n, dx, [-1.0, -1.0, -1.0]);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let y = field.node_pos(ix, iy, iz)[1];
                    field.set(ix, iy, iz, 7.0 * (y - 0.3));
                }
            }
        }
        LevelSetReinit::new(0.0, 100).fast_sweep(&mut field, 32);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let y = field.node_pos(ix, iy, iz)[1];
                    let value = field.get(ix, iy, iz);
                    // The plane y = 0.3 lies between two node layers, so the
                    // first-order interface values are accurate to a fraction of dx.
                    assert!(
                        (value - (y - 0.3)).abs() <= 0.25 * dx + 1e-12,
                        "phi = {value} vs {} at ({ix},{iy},{iz})",
                        y - 0.3
                    );
                }
            }
        }
    }

    #[test]
    fn converged_result_does_not_depend_on_the_sweep_budget() {
        let mut short = scaled_sphere(21, 0.5, 3.0);
        let mut long = scaled_sphere(21, 0.5, 3.0);
        let reinit = LevelSetReinit::new(0.0, 100);
        reinit.fast_sweep(&mut short, 64);
        reinit.fast_sweep(&mut long, 4000);
        for (a, b) in short.data.iter().zip(long.data.iter()) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }

    #[test]
    fn already_exact_signed_distance_is_preserved_to_first_order() {
        let (n, radius) = (33, 0.5);
        let mut field = scaled_sphere(n, radius, 1.0);
        let exact = field.data.clone();
        LevelSetReinit::new(1e-12, 100).fast_sweep(&mut field, 24);
        let dx = field.dx;
        for (a, b) in field.data.iter().zip(exact.iter()) {
            assert!((a - b).abs() <= dx, "{a} vs {b}");
        }
    }

    #[test]
    fn two_dimensional_slab_is_supported() {
        // A single layer along z: the missing axis must contribute nothing.
        let (nx, ny) = (41, 41);
        let dx = 0.05;
        let mut field = LevelSetField::new(nx, ny, 1, dx, [-1.0, -1.0, 0.0]);
        for iy in 0..ny {
            for ix in 0..nx {
                let p = field.node_pos(ix, iy, 0);
                let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
                field.set(ix, iy, 0, 4.0 * (r - 0.5));
            }
        }
        LevelSetReinit::new(1e-12, 100).fast_sweep(&mut field, 16);
        for iy in 5..ny - 5 {
            for ix in 5..nx - 5 {
                let p = field.node_pos(ix, iy, 0);
                let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
                if r > 0.15 && (r - 0.5).abs() > 0.1 {
                    let g = field.gradient(ix, iy, 0);
                    // `gradient` is not defined for a one-layer axis; use x and y.
                    let magnitude = (g[0] * g[0] + g[1] * g[1]).sqrt();
                    assert!((magnitude - 1.0).abs() <= 0.05, "|grad| = {magnitude}");
                }
            }
        }
    }

    #[test]
    fn partial_sweeps_leave_unreached_nodes_untouched() {
        // One sweep cannot reach every node of a large grid.
        let mut field = scaled_sphere(25, 0.3, 5.0);
        let before = field.data.clone();
        LevelSetReinit::new(0.0, 10).fast_sweep(&mut field, 1);
        assert!(field.data.iter().all(|v| v.is_finite() && v.abs() < 1.0e9));
        let untouched = field
            .data
            .iter()
            .zip(before.iter())
            .filter(|(a, b)| a == b)
            .count();
        assert!(untouched > 0, "some node should keep its original value");
        assert!(
            untouched < field.data.len(),
            "the sweep should update some nodes"
        );
    }

    #[test]
    fn fields_without_interface_or_with_nan_are_left_alone() {
        let mut positive = LevelSetField::new(4, 4, 4, 0.1, [0.0; 3]);
        positive.data.iter_mut().for_each(|v| *v = 3.0);
        positive.set(1, 1, 1, 0.5);
        let copy = positive.data.clone();
        LevelSetReinit::new(1e-9, 10).fast_sweep(&mut positive, 16);
        assert_eq!(positive.data, copy);

        let mut with_nan = scaled_sphere(9, 0.5, 2.0);
        with_nan.set(3, 3, 3, f64::NAN);
        let copy_bits: Vec<u64> = with_nan.data.iter().map(|v| v.to_bits()).collect();
        LevelSetReinit::new(1e-9, 10).fast_sweep(&mut with_nan, 16);
        let after_bits: Vec<u64> = with_nan.data.iter().map(|v| v.to_bits()).collect();
        assert_eq!(copy_bits, after_bits);

        let mut untouched = scaled_sphere(9, 0.5, 2.0);
        let copy = untouched.data.clone();
        LevelSetReinit::new(1e-9, 10).fast_sweep(&mut untouched, 0);
        assert_eq!(untouched.data, copy);
    }

    #[test]
    fn nodes_exactly_on_the_interface_stay_zero() {
        // A node with φ = 0 is part of the interface and must remain 0.
        let n = 9;
        let mut field = LevelSetField::new(n, n, n, 0.25, [-1.0, -1.0, -1.0]);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let x = field.node_pos(ix, iy, iz)[0];
                    field.set(ix, iy, iz, 5.0 * x);
                }
            }
        }
        LevelSetReinit::new(0.0, 10).fast_sweep(&mut field, 16);
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let x = field.node_pos(ix, iy, iz)[0];
                    let value = field.get(ix, iy, iz);
                    assert!((value - x).abs() < 1e-12, "phi = {value}, x = {x}");
                }
            }
        }
    }

    #[test]
    fn godunov_update_matches_the_closed_forms() {
        let h = 0.5;
        // One neighbour: u = a + h.
        assert!((godunov_update([1.0, f64::INFINITY, f64::INFINITY], h) - 1.5).abs() < 1e-15);
        // Two equal neighbours: u = a + h / sqrt(2).
        let two = godunov_update([1.0, 1.0, f64::INFINITY], h);
        assert!((two - (1.0 + h / 2.0_f64.sqrt())).abs() < 1e-15);
        // Three equal neighbours: u = a + h / sqrt(3).
        let three = godunov_update([1.0, 1.0, 1.0], h);
        assert!((three - (1.0 + h / 3.0_f64.sqrt())).abs() < 1e-15);
        // The solution satisfies the Eikonal discretisation.
        let [a, b, c] = [1.0, 1.2, 1.3];
        let u = godunov_update([a, b, c], h);
        let residual: f64 = [a, b, c].iter().map(|&x| (u - x).max(0.0).powi(2)).sum();
        assert!((residual - h * h).abs() < 1e-12, "residual {residual}");
    }
}
