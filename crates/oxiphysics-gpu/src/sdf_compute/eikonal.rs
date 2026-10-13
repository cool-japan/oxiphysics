// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eikonal solver and signed-distance re-initialisation for [`SdfGrid`].
//!
//! The distance to a set `Γ` solves the Eikonal equation `|∇d| = 1` with `d = 0`
//! on `Γ`.  This module solves it on the uniform grid with Zhao's fast sweeping
//! method (SIAM J. Numer. Anal. 2005): a Godunov upwind discretisation that is
//! relaxed by Gauss–Seidel sweeps in all `2³` axis orderings
//! (`±i`, `±j`, `±k`), which makes every characteristic direction available to
//! some sweep.
//!
//! * [`fast_sweeping_update`] performs one iteration (the eight sweeps) on a
//!   grid of known and unknown distances,
//! * [`fast_sweeping_solve`] repeats the iteration until the field stops
//!   changing,
//! * [`sdf_redistance`] rebuilds a true signed distance field from the zero
//!   level set of an arbitrary scalar field.
//!
//! Throughout, an *unknown* distance is stored as `±f64::MAX` (the fill value of
//! [`SdfGrid::new`]); the sign of an unknown cell is its side of the surface.

use rayon::prelude::*;

use super::functions::{grid_is_usable, point_triangle_distance, polygonize_cell};
use super::types::SdfGrid;

/// Working magnitude of a cell whose distance is not known yet.
const FAR: f64 = f64::INFINITY;

/// Upper bound on the full 8-sweep passes run by [`sdf_redistance`].
///
/// The Gauss–Seidel iteration stops changing after a handful of passes even for
/// strongly non-convex interfaces; the cap only guards against pathological
/// input.
const REDISTANCE_MAX_PASSES: usize = 64;

/// Convergence threshold of [`sdf_redistance`] as a fraction of the cell size.
const CONVERGENCE_TOLERANCE: f64 = 1e-12;

/// Godunov upwind solution of `|∇d| = 1` at one cell (Zhao 2005, constant speed).
///
/// `ux`, `uy` and `uz` are the smaller of the two axis neighbours along each
/// axis ([`FAR`] when a neighbour is missing or unknown) and `h` is the cell
/// size.  With the neighbour minima sorted as `a ≤ b ≤ c` the solution is
///
/// * `a + h` when that does not exceed `b` (information from one axis only),
/// * otherwise `(a + b + √(2h² − (a − b)²)) / 2` when that does not exceed `c`
///   (two axes),
/// * otherwise `(a + b + c + √((a + b + c)² − 3(a² + b² + c² − h²))) / 3`.
///
/// The formulas are evaluated relative to `a`, which keeps them well conditioned
/// when the distances are large.  Returns [`FAR`] when no neighbour is known.
#[inline]
fn godunov_update(ux: f64, uy: f64, uz: f64, h: f64) -> f64 {
    let (mut a, mut b, mut c) = (ux, uy, uz);
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    if b > c {
        std::mem::swap(&mut b, &mut c);
    }
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    if !a.is_finite() {
        return FAR;
    }
    let one_axis = a + h;
    if one_axis <= b {
        return one_axis;
    }
    let db = b - a;
    let two_axes = a + 0.5 * (db + (2.0 * h * h - db * db).max(0.0).sqrt());
    if two_axes <= c {
        return two_axes;
    }
    let dc = c - a;
    let sum = db + dc;
    let discriminant = (sum * sum - 3.0 * (db * db + dc * dc - h * h)).max(0.0);
    a + (sum + discriminant.sqrt()) / 3.0
}

/// Smaller of the two neighbours of `idx` along one axis ([`FAR`] if absent).
#[inline]
fn axis_neighbour_min(
    u: &[f64],
    idx: usize,
    has_lower: bool,
    has_upper: bool,
    stride: usize,
) -> f64 {
    let lower = if has_lower { u[idx - stride] } else { FAR };
    let upper = if has_upper { u[idx + stride] } else { FAR };
    lower.min(upper)
}

/// One full iteration of the fast sweeping method: the eight Gauss–Seidel sweeps
/// over every combination of ascending/descending `i`, `j` and `k`.
///
/// `u` holds the non-negative working distances (in the grid's flat layout, see
/// [`SdfGrid::index`]) and `frozen` marks cells that are boundary data and must
/// not change.  Values are only ever lowered.  Returns the largest decrease of
/// any cell during the iteration ([`FAR`] if a cell went from unknown to known).
fn sweep_pass(u: &mut [f64], frozen: &[bool], dims: [usize; 3], h: f64) -> f64 {
    let [nx, ny, nz] = dims;
    let stride_i = ny * nz;
    let stride_j = nz;
    let mut max_change = 0.0_f64;
    for octant in 0..8_u8 {
        let reverse_i = octant & 1 != 0;
        let reverse_j = octant & 2 != 0;
        let reverse_k = octant & 4 != 0;
        for step_i in 0..nx {
            let i = if reverse_i { nx - 1 - step_i } else { step_i };
            for step_j in 0..ny {
                let j = if reverse_j { ny - 1 - step_j } else { step_j };
                let row = i * stride_i + j * stride_j;
                for step_k in 0..nz {
                    let k = if reverse_k { nz - 1 - step_k } else { step_k };
                    let idx = row + k;
                    if frozen[idx] {
                        continue;
                    }
                    let ux = axis_neighbour_min(u, idx, i > 0, i + 1 < nx, stride_i);
                    let uy = axis_neighbour_min(u, idx, j > 0, j + 1 < ny, stride_j);
                    let uz = axis_neighbour_min(u, idx, k > 0, k + 1 < nz, 1);
                    let candidate = godunov_update(ux, uy, uz, h);
                    let old = u[idx];
                    if candidate < old {
                        max_change = max_change.max(old - candidate);
                        u[idx] = candidate;
                    }
                }
            }
        }
    }
    max_change
}

/// Repeat [`sweep_pass`] until the largest change of an iteration is at most
/// `tolerance` or `max_passes` iterations have run.  Returns the iteration count.
fn solve_magnitudes(
    u: &mut [f64],
    frozen: &[bool],
    dims: [usize; 3],
    h: f64,
    max_passes: usize,
    tolerance: f64,
) -> usize {
    let mut passes = 0;
    while passes < max_passes {
        let change = sweep_pass(u, frozen, dims, h);
        passes += 1;
        if change <= tolerance {
            break;
        }
    }
    passes
}

/// Side (`true` = inside, value `< 0`) and working magnitude of every cell.
///
/// Unknown cells (`±f64::MAX`, infinities and NaN) get the magnitude [`FAR`].
fn split_signed(values: &[f64]) -> (Vec<bool>, Vec<f64>) {
    let negative = values.iter().map(|&v| v < 0.0).collect();
    let magnitude = values
        .iter()
        .map(|&v| {
            let m = v.abs();
            if m.is_finite() && m < f64::MAX {
                m
            } else {
                FAR
            }
        })
        .collect();
    (negative, magnitude)
}

/// Cells that are boundary data of [`fast_sweeping_update`]: known cells that
/// have an axis neighbour on the other side of the surface, and exact zeros.
fn boundary_data_mask(
    values: &[f64],
    negative: &[bool],
    magnitude: &[f64],
    dims: [usize; 3],
) -> Vec<bool> {
    let [nx, ny, nz] = dims;
    let stride_i = ny * nz;
    let stride_j = nz;
    let mut frozen = vec![false; values.len()];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = i * stride_i + j * stride_j + k;
                if !magnitude[idx].is_finite() {
                    continue;
                }
                let side = negative[idx];
                let crosses = (i > 0 && negative[idx - stride_i] != side)
                    || (i + 1 < nx && negative[idx + stride_i] != side)
                    || (j > 0 && negative[idx - stride_j] != side)
                    || (j + 1 < ny && negative[idx + stride_j] != side)
                    || (k > 0 && negative[idx - 1] != side)
                    || (k + 1 < nz && negative[idx + 1] != side);
                frozen[idx] = crosses || values[idx] == 0.0;
            }
        }
    }
    frozen
}

/// One iteration of Zhao's fast sweeping method on a signed distance grid.
///
/// The grid holds known distances and unknown cells (`±f64::MAX`, the fill value
/// of [`SdfGrid::new`]).  The sign of every cell is stored and the sweeps work on
/// the magnitudes `|φ|`: each of the eight Gauss–Seidel sweeps (all `±i`, `±j`,
/// `±k` orderings) lowers the magnitude of a cell to the Godunov upwind solution
/// of `|∇d| = 1` built from its axis neighbours.  With `a ≤ b ≤ c` the smaller
/// neighbour along each axis and `h` the cell size, that solution is `a + h` if
/// this does not exceed `b`, otherwise the two-axis root
/// `(a + b + √(2h² − (a − b)²)) / 2` if this does not exceed `c`, otherwise the
/// three-axis root of `Σ (x − ·)² = h²`.  The sign is restored afterwards.  Known
/// cells adjacent to a sign change and cells that are exactly zero are boundary
/// data and are never modified; all other values are only ever lowered, so the
/// grid acts as an upper bound that is tightened by the iteration.
///
/// A single iteration already reaches the converged discrete solution when the
/// characteristics are straight lines (a point source or a convex front); call
/// [`fast_sweeping_solve`] to iterate until the field has converged.  The scheme
/// is first order, so away from the axes it overestimates the Euclidean distance
/// (at the cell offset `(3, 4, 0)` from a point seed it yields `5.53 h` instead of
/// `5 h`, where a city-block relaxation would give `7 h`).  To rebuild a distance
/// field from an arbitrary scalar field use [`sdf_redistance`].
pub fn fast_sweeping_update(grid: &mut SdfGrid) {
    fast_sweeping_solve(grid, 1, 0.0);
}

/// Run [`fast_sweeping_update`] until the field has converged.
///
/// Iterates until the largest change of one iteration is at most `tolerance`
/// (in world units) or `max_passes` iterations have run, and returns the number
/// of iterations performed.  The Gauss–Seidel iteration decreases the values
/// monotonically, so it always terminates; `tolerance = 0.0` iterates until no
/// value changes at all.  A grid with inconsistent dimensions or a non-positive
/// cell size is left untouched and `0` is returned.
pub fn fast_sweeping_solve(grid: &mut SdfGrid, max_passes: usize, tolerance: f64) -> usize {
    if !grid_is_usable(grid) {
        return 0;
    }
    let dims = [grid.nx, grid.ny, grid.nz];
    let (negative, original) = split_signed(&grid.values);
    let frozen = boundary_data_mask(&grid.values, &negative, &original, dims);
    let mut magnitude = original.clone();
    let passes = solve_magnitudes(
        &mut magnitude,
        &frozen,
        dims,
        grid.dx,
        max_passes,
        tolerance,
    );
    for ((value, &is_negative), (&solved, &before)) in grid
        .values
        .iter_mut()
        .zip(&negative)
        .zip(magnitude.iter().zip(&original))
    {
        if solved < before {
            *value = if is_negative { -solved } else { solved };
        }
    }
    passes
}

/// Initial distance of grid node `idx` to the zero level set, and whether that
/// distance is exact boundary data.
///
/// Two kinds of nodes are next to the surface:
///
/// * *Crossing nodes* have an axis neighbour on the other side of the surface (or
///   are exactly zero).  The surface passes within one cell of them, so their
///   exact distance to the piecewise-linear surface is the minimum over the
///   marching-cubes triangles of the (up to eight) cells that share the node, and
///   at most the interpolated zero crossing on a sign-changing axis edge (which
///   is also the only estimate on grids that are one cell thick).  These values
///   are returned as *frozen* boundary data of the Eikonal solve.
/// * *Cut-cell corners* belong to a cell that the surface cuts but have no
///   sign-changing axis edge.  The surface may be as far as a cell diagonal away
///   and its nearest point need not lie in the cells that share the node, so the
///   triangle distance is only an upper bound.  It is returned as a starting value
///   that the sweeps may still lower.
///
/// Every other node gets [`FAR`] and is not frozen.
fn node_band_distance(grid: &SdfGrid, negative: &[bool], idx: usize) -> (f64, bool) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let stride_i = ny * nz;
    let k = idx % nz;
    let j = (idx / nz) % ny;
    let i = idx / stride_i;
    let h = grid.dx;
    let value = grid.values[idx];
    let own = value.abs();

    let mut neighbours = [0_usize; 6];
    let mut count = 0;
    if i > 0 {
        neighbours[count] = idx - stride_i;
        count += 1;
    }
    if i + 1 < nx {
        neighbours[count] = idx + stride_i;
        count += 1;
    }
    if j > 0 {
        neighbours[count] = idx - nz;
        count += 1;
    }
    if j + 1 < ny {
        neighbours[count] = idx + nz;
        count += 1;
    }
    if k > 0 {
        neighbours[count] = idx - 1;
        count += 1;
    }
    if k + 1 < nz {
        neighbours[count] = idx + 1;
        count += 1;
    }
    let mut crossing = if value == 0.0 { 0.0 } else { FAR };
    for &other_idx in &neighbours[..count] {
        if negative[other_idx] != negative[idx] {
            let sum = own + grid.values[other_idx].abs();
            if sum.is_finite() && sum > 0.0 {
                crossing = crossing.min(h * own / sum);
            }
        }
    }
    let is_crossing_node = crossing.is_finite();

    let mut best = crossing;
    if nx >= 2 && ny >= 2 && nz >= 2 {
        let position = grid.world_pos(i, j, k);
        for di in 0..2_usize {
            for dj in 0..2_usize {
                for dk in 0..2_usize {
                    let (Some(ci), Some(cj), Some(ck)) =
                        (i.checked_sub(di), j.checked_sub(dj), k.checked_sub(dk))
                    else {
                        continue;
                    };
                    if ci + 1 >= nx || cj + 1 >= ny || ck + 1 >= nz {
                        continue;
                    }
                    let cell = polygonize_cell(grid, ci, cj, ck, 0.0);
                    for triangle in &cell.triangles[..cell.count] {
                        let distance = point_triangle_distance(
                            &position,
                            &triangle[0],
                            &triangle[1],
                            &triangle[2],
                        );
                        if distance.is_finite() {
                            best = best.min(distance);
                        }
                    }
                }
            }
        }
    }
    (best, is_crossing_node)
}

/// Rebuild a signed distance field from the zero level set of `grid`.
///
/// Only the *sign pattern and the zero crossings* of the input are used; its
/// magnitudes away from the surface are discarded, which is what makes the
/// function suitable after operations that destroy the distance property
/// (uniform offsets, blending, CSG of a thin feature, noise):
///
/// 1. the nodes next to the zero level set are initialised with their exact
///    distance to the piecewise-linear surface that marching cubes extracts, see
///    [`marching_cubes`](super::marching_cubes); nodes with a sign-changing axis
///    neighbour keep that value as fixed boundary data, the remaining corners of
///    cells that the surface cuts only start from it;
/// 2. every other cell starts as "unknown" while the Eikonal equation `|∇d| = 1`
///    is solved outwards on both sides of the surface by fast sweeping
///    ([`fast_sweeping_solve`]) until the field stops changing;
/// 3. the original sign (`value < 0` is inside) is restored.
///
/// The result is a unit-gradient field whose zero level set coincides with the
/// input's.  If the input has no zero crossing at all (every value on one side)
/// there is no surface to measure the distance to, and every cell becomes
/// `±f64::MAX` on its side, i.e. an empty (or completely filled) shape.  A grid
/// with inconsistent dimensions or a non-positive cell size is returned
/// unchanged.
pub fn sdf_redistance(grid: &SdfGrid) -> SdfGrid {
    let rebuild = |values: Vec<f64>| SdfGrid {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        dx: grid.dx,
        origin: grid.origin,
        values,
    };
    if !grid_is_usable(grid) {
        return rebuild(grid.values.clone());
    }
    let dims = [grid.nx, grid.ny, grid.nz];
    let negative: Vec<bool> = grid.values.par_iter().map(|&v| v < 0.0).collect();
    let (mut magnitude, frozen): (Vec<f64>, Vec<bool>) = (0..grid.values.len())
        .into_par_iter()
        .map(|idx| node_band_distance(grid, &negative, idx))
        .unzip();
    if !magnitude.iter().any(|d| d.is_finite()) {
        return rebuild(
            negative
                .iter()
                .map(|&is_negative| if is_negative { -f64::MAX } else { f64::MAX })
                .collect(),
        );
    }
    solve_magnitudes(
        &mut magnitude,
        &frozen,
        dims,
        grid.dx,
        REDISTANCE_MAX_PASSES,
        CONVERGENCE_TOLERANCE * grid.dx,
    );
    rebuild(
        magnitude
            .iter()
            .zip(&negative)
            .map(|(&m, &is_negative)| {
                let m = if m.is_finite() { m } else { f64::MAX };
                if is_negative { -m } else { m }
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: f64 = 0.5;

    /// A cubic grid of unknown distances (`+f64::MAX`) with a point seed of value
    /// `0` at the centre cell.
    fn point_seed_grid(n: usize) -> (SdfGrid, usize) {
        let mut grid = SdfGrid::new(n, n, n, H, [0.0; 3]);
        let c = n / 2;
        grid.set(c, c, c, 0.0);
        (grid, c)
    }

    fn at(grid: &SdfGrid, c: usize, di: i32, dj: i32, dk: i32) -> f64 {
        let ix = |d: i32| usize::try_from(c as i32 + d).expect("offset stays inside the grid");
        grid.get(ix(di), ix(dj), ix(dk))
    }

    #[test]
    fn godunov_update_matches_the_closed_forms() {
        let h = 0.25;
        // One informative axis: a + h.
        assert!((godunov_update(1.0, FAR, FAR, h) - 1.25).abs() < 1e-15);
        assert!((godunov_update(FAR, 1.0, FAR, h) - 1.25).abs() < 1e-15);
        // Neighbour spread of exactly h: still the one-axis solution.
        assert!((godunov_update(1.0, 1.25, FAR, h) - 1.25).abs() < 1e-15);
        // Two equal axes: a + h / sqrt(2).
        let two = godunov_update(2.0, 2.0, FAR, h);
        assert!((two - (2.0 + h * std::f64::consts::FRAC_1_SQRT_2)).abs() < 1e-14);
        // Three equal axes: a + h / sqrt(3).
        let three = godunov_update(2.0, 2.0, 2.0, h);
        assert!((three - (2.0 + h / 3.0_f64.sqrt())).abs() < 1e-14);
        // The result solves the quadratic: sum of squared positive parts is h^2.
        let (a, b, c) = (1.0, 1.1, 1.15);
        let x = godunov_update(a, b, c, h);
        let residual = (x - a).powi(2) + (x - b).powi(2) + (x - c).powi(2);
        assert!((residual - h * h).abs() < 1e-12, "residual {residual}");
        // Argument order does not matter.
        assert_eq!(godunov_update(c, a, b, h), godunov_update(a, b, c, h));
    }

    #[test]
    fn godunov_update_without_known_neighbours_stays_far() {
        assert_eq!(godunov_update(FAR, FAR, FAR, 1.0), FAR);
    }

    #[test]
    fn point_seed_distances_are_exact_on_axes_and_far_below_manhattan() {
        let (mut grid, c) = point_seed_grid(21);
        fast_sweeping_update(&mut grid);
        for step in 1..=8_i32 {
            let want = f64::from(step) * H;
            assert!(
                (at(&grid, c, step, 0, 0) - want).abs() < 1e-12,
                "+x axis {step}"
            );
            assert!(
                (at(&grid, c, -step, 0, 0) - want).abs() < 1e-12,
                "-x axis {step}"
            );
            assert!(
                (at(&grid, c, 0, step, 0) - want).abs() < 1e-12,
                "+y axis {step}"
            );
            assert!(
                (at(&grid, c, 0, 0, -step) - want).abs() < 1e-12,
                "-z axis {step}"
            );
        }
        // Offset (3, 4, 0): the Euclidean distance is 5h and the city-block
        // (L1) value that the old chamfer produced is 7h.  First-order Godunov
        // overestimates the Euclidean distance near a point source, converging
        // to 5.530h here, still far from 7h.
        let d340 = at(&grid, c, 3, 4, 0);
        assert!(d340 > 5.0 * H, "never below the true distance, got {d340}");
        assert!(
            d340 < 5.6 * H,
            "must beat the L1 value 7h by far, got {d340}"
        );
        assert!((d340 - 5.530_022_893 * H).abs() < 1e-6 * H, "got {d340}");
        for (offset, want) in [
            ((1, 1, 0), 1.707_106_781),
            ((2, 1, 0), 2.545_328_925),
            ((3, 2, 0), 4.048_043_049),
            ((1, 1, 1), 2.284_457_050),
            ((2, 2, 2), 4.243_559_041),
        ] {
            let got = at(&grid, c, offset.0, offset.1, offset.2);
            assert!(
                (got - want * H).abs() < 1e-6 * H,
                "offset {offset:?}: got {got}, want {}",
                want * H
            );
        }
        // Symmetry under axis permutations and reflections.
        let mirrored = [
            at(&grid, c, 4, 3, 0),
            at(&grid, c, 0, 3, 4),
            at(&grid, c, -3, 0, -4),
            at(&grid, c, 0, -4, 3),
        ];
        for value in mirrored {
            assert!((value - d340).abs() < 1e-12, "{value} vs {d340}");
        }
    }

    #[test]
    fn one_iteration_converges_for_a_point_source() {
        let (mut grid, _) = point_seed_grid(17);
        fast_sweeping_update(&mut grid);
        let after_one = grid.values.clone();
        fast_sweeping_update(&mut grid);
        assert_eq!(
            grid.values, after_one,
            "second iteration must change nothing"
        );
    }

    #[test]
    fn solve_reports_the_iteration_count_and_respects_the_cap() {
        let (mut grid, _) = point_seed_grid(11);
        let passes = fast_sweeping_solve(&mut grid, 50, 0.0);
        assert_eq!(
            passes, 2,
            "one productive iteration plus the confirming one"
        );
        let (mut capped, _) = point_seed_grid(11);
        assert_eq!(fast_sweeping_solve(&mut capped, 1, 0.0), 1);
        assert_eq!(fast_sweeping_solve(&mut capped, 0, 0.0), 0);
        assert_eq!(capped.values, grid.values);
    }

    #[test]
    fn unreached_cells_keep_their_unknown_value() {
        // No seed at all: nothing to propagate, nothing changes.
        let mut grid = SdfGrid::new(5, 5, 5, H, [0.0; 3]);
        fast_sweeping_update(&mut grid);
        assert!(grid.values.iter().all(|&v| v == f64::MAX));
        // Unknown inside cells keep their negative sentinel.
        let mut inside = SdfGrid::new(4, 4, 4, H, [0.0; 3]);
        inside.values.iter_mut().for_each(|v| *v = -f64::MAX);
        fast_sweeping_update(&mut inside);
        assert!(inside.values.iter().all(|&v| v == -f64::MAX));
    }

    #[test]
    fn sweeping_preserves_signs_and_freezes_the_sign_change() {
        // A slab: x < 1.2 is inside.  Cells next to the sign change carry the
        // distance to the interface; cells further away hold unknown values.
        let n = 12;
        let mut grid = SdfGrid::new(n, 3, 3, H, [0.0; 3]);
        let interface = 6.0 * H;
        for i in 0..n {
            for j in 0..3 {
                for k in 0..3 {
                    let x = (i as f64 + 0.5) * H;
                    let signed = x - interface;
                    let near = i == 5 || i == 6;
                    // The frozen cells deliberately hold half the true distance: they
                    // are boundary data and must come out of the sweeps unchanged.
                    let value = if near {
                        signed * 0.5
                    } else if x < interface {
                        -f64::MAX
                    } else {
                        f64::MAX
                    };
                    grid.set(i, j, k, value);
                }
            }
        }
        let before_frozen = (grid.get(5, 1, 1), grid.get(6, 1, 1));
        fast_sweeping_solve(&mut grid, 10, 0.0);
        assert_eq!((grid.get(5, 1, 1), grid.get(6, 1, 1)), before_frozen);
        for i in 0..n {
            let inside = i <= 5;
            for j in 0..3 {
                for k in 0..3 {
                    let v = grid.get(i, j, k);
                    assert_eq!(v < 0.0, inside, "sign of cell {i}");
                    assert!(v.abs() < f64::MAX, "cell {i} must have received a distance");
                }
            }
        }
        // Moving away from the frozen band each side grows by exactly h per cell.
        assert!((grid.get(7, 1, 1) - (grid.get(6, 1, 1) + H)).abs() < 1e-12);
        assert!((grid.get(4, 1, 1) - (grid.get(5, 1, 1) - H)).abs() < 1e-12);
        assert!((grid.get(0, 1, 1) - (grid.get(5, 1, 1) - 5.0 * H)).abs() < 1e-12);
    }

    #[test]
    fn exact_zero_cells_are_seeds() {
        let mut grid = SdfGrid::new(9, 9, 9, 1.0, [0.0; 3]);
        grid.set(4, 4, 4, 0.0);
        grid.set(0, 0, 0, 0.0);
        fast_sweeping_update(&mut grid);
        assert_eq!(grid.get(4, 4, 4), 0.0);
        assert_eq!(grid.get(0, 0, 0), 0.0);
        assert!((grid.get(4, 4, 7) - 3.0).abs() < 1e-12);
        assert!(
            (grid.get(0, 0, 2) - 2.0).abs() < 1e-12,
            "nearer to the corner seed"
        );
    }

    #[test]
    fn inconsistent_or_degenerate_grids_are_left_untouched() {
        let mut broken = SdfGrid::new(3, 3, 3, 1.0, [0.0; 3]);
        broken.values.truncate(10);
        assert_eq!(fast_sweeping_solve(&mut broken, 5, 0.0), 0);
        let mut zero_spacing = SdfGrid::new(3, 3, 3, 0.0, [0.0; 3]);
        zero_spacing.set(1, 1, 1, 0.0);
        assert_eq!(fast_sweeping_solve(&mut zero_spacing, 5, 0.0), 0);
        assert_eq!(zero_spacing.get(0, 0, 0), f64::MAX);
        let mut empty = SdfGrid::new(0, 4, 4, 1.0, [0.0; 3]);
        assert_eq!(fast_sweeping_solve(&mut empty, 5, 0.0), 0);
    }

    #[test]
    fn thin_grids_sweep_along_the_available_axes() {
        // A single row: pure one-dimensional propagation.
        let mut row = SdfGrid::new(9, 1, 1, 0.25, [0.0; 3]);
        row.set(2, 0, 0, 0.0);
        fast_sweeping_update(&mut row);
        for i in 0..9 {
            let want = (i as f64 - 2.0).abs() * 0.25;
            assert!((row.get(i, 0, 0) - want).abs() < 1e-12, "cell {i}");
        }
    }

    #[test]
    fn redistance_of_a_grid_without_surface_is_empty_or_full() {
        let all_outside = SdfGrid {
            nx: 4,
            ny: 4,
            nz: 4,
            dx: 0.1,
            origin: [0.0; 3],
            values: vec![0.3; 64],
        };
        let empty = sdf_redistance(&all_outside);
        assert!(empty.values.iter().all(|&v| v == f64::MAX));
        let all_inside = SdfGrid {
            values: vec![-0.3; 64],
            ..SdfGrid::new(4, 4, 4, 0.1, [0.0; 3])
        };
        let full = sdf_redistance(&all_inside);
        assert!(full.values.iter().all(|&v| v == -f64::MAX));
    }

    #[test]
    fn redistance_of_a_planar_interface_is_exact() {
        // phi = 3 * (x - x0): wrong gradient magnitude, planar zero set.
        let (n, dx) = (14, 0.1);
        let x0 = 0.713;
        let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let x = grid.world_pos(i, j, k)[0];
                    grid.set(i, j, k, 3.0 * (x - x0));
                }
            }
        }
        let fixed = sdf_redistance(&grid);
        for i in 0..n {
            let x = fixed.world_pos(i, 0, 0)[0];
            for j in 0..n {
                for k in 0..n {
                    let got = fixed.get(i, j, k);
                    assert!(
                        (got - (x - x0)).abs() < 1e-9,
                        "cell ({i},{j},{k}): got {got}, want {}",
                        x - x0
                    );
                }
            }
        }
    }

    #[test]
    fn redistance_handles_a_single_cell_thick_slab() {
        // 2-D disc of radius 0.5 in a one-cell-thick grid: only the axis-edge
        // crossings are available, so accuracy is first order.
        let (n, dx) = (41, 0.05);
        let centre = [1.02, 1.0];
        let mut grid = SdfGrid::new(n, n, 1, dx, [0.0; 3]);
        for i in 0..n {
            for j in 0..n {
                let p = grid.world_pos(i, j, 0);
                let r = ((p[0] - centre[0]).powi(2) + (p[1] - centre[1]).powi(2)).sqrt();
                grid.set(i, j, 0, 7.0 * (r - 0.5));
            }
        }
        let fixed = sdf_redistance(&grid);
        let mut worst = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                let p = fixed.world_pos(i, j, 0);
                let r = ((p[0] - centre[0]).powi(2) + (p[1] - centre[1]).powi(2)).sqrt();
                let got = fixed.get(i, j, 0);
                assert_eq!(got < 0.0, r < 0.5, "sign at ({i},{j})");
                worst = worst.max((got - (r - 0.5)).abs());
            }
        }
        assert!(worst < 0.6 * dx, "worst error {worst} exceeds 0.6 dx");
    }
}
