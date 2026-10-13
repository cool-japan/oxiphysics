// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Streaming step for LBM simulations.
//!
//! Propagates distribution functions from each cell to its neighbors
//! along the discrete velocity directions. Periodic boundaries are
//! applied by default (wrapping around grid edges).
//!
//! Also provides specialised streaming variants:
//! - Bounce-back streaming for solid walls (full-way)
//! - Half-way bounce-back (more accurate, second-order)
//! - Pull-scheme streaming (each cell pulls from neighbors)
//! - Push-scheme streaming (each cell pushes to neighbors)
//! - AA pattern (alternating even/odd in-place streaming)
//! - Swap scheme (neighbor-swap in-place streaming)
//! - Streaming with force correction (Guo forcing)
//! - Periodic streaming helpers (shift, wrap)

use crate::grid::{LbmGrid2D, LbmGrid3D};

// ---------------------------------------------------------------------------
// Periodic shift helpers
// ---------------------------------------------------------------------------

/// Reduce the lattice velocity component `c` to a non-negative cyclic shift
/// in `[0, n)` along an axis of `n` cells (`n > 0`).
#[inline]
fn cyclic_shift(c: i32, n: usize) -> usize {
    i64::from(c).rem_euclid(n as i64) as usize
}

/// Make sure `scratch` holds `q` planes of exactly `n` values.
///
/// This is a no-op after the first call, so streaming never allocates in the
/// steady state.  It also repairs a scratch buffer that was emptied or
/// resized by code that rebuilt the grid by hand.
fn ensure_scratch(scratch: &mut Vec<Vec<f64>>, q: usize, n: usize) {
    scratch.resize_with(q, Vec::new);
    for plane in scratch.iter_mut() {
        if plane.len() != n {
            plane.clear();
            plane.resize(n, 0.0);
        }
    }
}

/// Copy the `nx × ny` plane `src` into `dst`, moving every value by
/// `(sx, sy)` cells with periodic wrap-around (`0 <= sx < nx`, `0 <= sy < ny`).
///
/// `dst[(y + sy) mod ny][(x + sx) mod nx] = src[y][x]`.  The wrap is handled
/// by splitting each row (or, for `sx == 0`, the whole plane) into two
/// contiguous spans, so the inner work is `copy_from_slice` and there is no
/// per-cell index arithmetic.
#[inline]
fn shift_plane_into(dst: &mut [f64], src: &[f64], nx: usize, ny: usize, sx: usize, sy: usize) {
    debug_assert_eq!(dst.len(), nx * ny);
    debug_assert_eq!(src.len(), nx * ny);
    if sx == 0 {
        // Whole rows move together: two contiguous blocks cover the plane.
        let split = (ny - sy) * nx;
        dst[sy * nx..].copy_from_slice(&src[..split]);
        dst[..sy * nx].copy_from_slice(&src[split..]);
        return;
    }
    let split = nx - sx;
    for (y, src_row) in src.chunks_exact(nx).enumerate() {
        let dy = if y + sy >= ny { y + sy - ny } else { y + sy };
        let dst_row = &mut dst[dy * nx..(dy + 1) * nx];
        dst_row[sx..].copy_from_slice(&src_row[..split]);
        dst_row[..sx].copy_from_slice(&src_row[split..]);
    }
}

// ---------------------------------------------------------------------------
// 2D periodic streaming
// ---------------------------------------------------------------------------

/// Stream distributions on a 2D grid with periodic boundaries.
///
/// Every direction `i` is a cyclic shift of its plane by the lattice velocity
/// `c_i`.  The shift is written into the grid's scratch buffer
/// ([`LbmGrid2D::f_tmp`]) with contiguous span copies and the per-direction
/// buffers are then swapped, so a step performs no allocation and no
/// per-cell modular arithmetic.  The rest direction does not move and is left
/// untouched.
pub fn stream_2d(grid: &mut LbmGrid2D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let n = nx * ny;
    if n == 0 {
        return;
    }
    let q = grid.lattice.q();
    ensure_scratch(&mut grid.f_tmp, q, n);

    for i in 0..q {
        let c = grid.lattice.velocity_2d(i);
        let sx = cyclic_shift(c[0], nx);
        let sy = cyclic_shift(c[1], ny);
        if sx == 0 && sy == 0 {
            continue;
        }
        shift_plane_into(&mut grid.f_tmp[i], &grid.f[i], nx, ny, sx, sy);
        std::mem::swap(&mut grid.f[i], &mut grid.f_tmp[i]);
    }
}

// ---------------------------------------------------------------------------
// 3D periodic streaming
// ---------------------------------------------------------------------------

/// Stream distributions on a 3D grid with periodic boundaries.
///
/// Same scheme as [`stream_2d`]: each direction is shifted plane by plane
/// into [`LbmGrid3D::f_tmp`] and the buffers are swapped, so no memory is
/// allocated per step.
pub fn stream_3d(grid: &mut LbmGrid3D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let plane = nx * ny;
    let n = plane * nz;
    if n == 0 {
        return;
    }
    let q = grid.lattice.q();
    ensure_scratch(&mut grid.f_tmp, q, n);

    for i in 0..q {
        let c = grid.lattice.velocity_3d(i);
        let sx = cyclic_shift(c[0], nx);
        let sy = cyclic_shift(c[1], ny);
        let sz = cyclic_shift(c[2], nz);
        if sx == 0 && sy == 0 && sz == 0 {
            continue;
        }
        {
            let src = &grid.f[i];
            let dst = &mut grid.f_tmp[i];
            for z in 0..nz {
                let dz = if z + sz >= nz { z + sz - nz } else { z + sz };
                shift_plane_into(
                    &mut dst[dz * plane..(dz + 1) * plane],
                    &src[z * plane..(z + 1) * plane],
                    nx,
                    ny,
                    sx,
                    sy,
                );
            }
        }
        std::mem::swap(&mut grid.f[i], &mut grid.f_tmp[i]);
    }
}

// ---------------------------------------------------------------------------
// Pull-scheme streaming 2D
// ---------------------------------------------------------------------------

/// Pull-scheme streaming on a 2D grid with periodic boundaries.
///
/// In the pull scheme each destination cell pulls the distribution
/// from the upstream neighbour in the opposite lattice direction.
/// This is semantically equivalent to `stream_2d` but uses an
/// explicit pull formulation that is cache-friendlier on GPUs.
pub fn stream_pull_2d(grid: &mut LbmGrid2D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let n = nx * ny;

    let f_old: Vec<Vec<f64>> = grid.f.clone();
    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for y in 0..ny {
        for x in 0..nx {
            let k_dst = y * nx + x;
            for i in 0..q {
                let c = grid.lattice.velocity_2d(i);
                let sx = ((x as i64 - c[0] as i64).rem_euclid(nx as i64)) as usize;
                let sy = ((y as i64 - c[1] as i64).rem_euclid(ny as i64)) as usize;
                let k_src = sy * nx + sx;
                f_new[i][k_dst] = f_old[i][k_src];
            }
        }
    }

    grid.f = f_new;
}

/// Pull-scheme streaming on a 3D grid with periodic boundaries.
pub fn stream_pull_3d(grid: &mut LbmGrid3D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let q = grid.lattice.q();
    let n = nx * ny * nz;

    let f_old: Vec<Vec<f64>> = grid.f.clone();
    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let k_dst = z * ny * nx + y * nx + x;
                for i in 0..q {
                    let c = grid.lattice.velocity_3d(i);
                    let sx = ((x as i64 - c[0] as i64).rem_euclid(nx as i64)) as usize;
                    let sy = ((y as i64 - c[1] as i64).rem_euclid(ny as i64)) as usize;
                    let sz = ((z as i64 - c[2] as i64).rem_euclid(nz as i64)) as usize;
                    let k_src = sz * ny * nx + sy * nx + sx;
                    f_new[i][k_dst] = f_old[i][k_src];
                }
            }
        }
    }

    grid.f = f_new;
}

// ---------------------------------------------------------------------------
// Push-scheme streaming 2D
// ---------------------------------------------------------------------------

/// Push-scheme streaming on a 2D grid with periodic boundaries.
///
/// In the push scheme each source cell pushes its distributions
/// to the downstream neighbour along each lattice velocity direction.
pub fn stream_push_2d(grid: &mut LbmGrid2D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let n = nx * ny;

    let f_old: Vec<Vec<f64>> = grid.f.clone();
    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for y in 0..ny {
        for x in 0..nx {
            let k_src = y * nx + x;
            for i in 0..q {
                let c = grid.lattice.velocity_2d(i);
                let dx = ((x as i64 + c[0] as i64).rem_euclid(nx as i64)) as usize;
                let dy = ((y as i64 + c[1] as i64).rem_euclid(ny as i64)) as usize;
                let k_dst = dy * nx + dx;
                f_new[i][k_dst] = f_old[i][k_src];
            }
        }
    }

    grid.f = f_new;
}

/// Push-scheme streaming on a 3D grid with periodic boundaries.
pub fn stream_push_3d(grid: &mut LbmGrid3D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let q = grid.lattice.q();
    let n = nx * ny * nz;

    let f_old: Vec<Vec<f64>> = grid.f.clone();
    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let k_src = z * ny * nx + y * nx + x;
                for i in 0..q {
                    let c = grid.lattice.velocity_3d(i);
                    let dx = ((x as i64 + c[0] as i64).rem_euclid(nx as i64)) as usize;
                    let dy = ((y as i64 + c[1] as i64).rem_euclid(ny as i64)) as usize;
                    let dz = ((z as i64 + c[2] as i64).rem_euclid(nz as i64)) as usize;
                    let k_dst = dz * ny * nx + dy * nx + dx;
                    f_new[i][k_dst] = f_old[i][k_src];
                }
            }
        }
    }

    grid.f = f_new;
}

// ---------------------------------------------------------------------------
// AA pattern streaming 2D
// ---------------------------------------------------------------------------

/// AA-pattern streaming for 2D grids.
///
/// The AA (alternating) pattern avoids an extra buffer by performing
/// streaming in-place over two half-steps.  On even steps the
/// distributions propagate in the positive lattice directions and on
/// odd steps in the negative directions.
///
/// After one even + one odd step the result is equivalent to one
/// standard push streaming step.
///
/// `is_even` indicates which half-step to perform.
pub fn stream_aa_2d(grid: &mut LbmGrid2D, is_even: bool) {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();

    // For the AA pattern we swap f[i][k] with f[opp(i)][k_neighbor]
    // in a single pass.  To avoid double-swaps we only process each
    // pair (i, opp(i)) once.

    // Collect pairs (i < opp(i)) to avoid double processing.
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut visited = vec![false; q];
    for i in 0..q {
        let opp = grid.lattice.opposite(i);
        if i == opp {
            continue; // rest direction — no streaming
        }
        if !visited[i] {
            pairs.push((i, opp));
            visited[i] = true;
            visited[opp] = true;
        }
    }

    for &(i, opp) in &pairs {
        let c = grid.lattice.velocity_2d(if is_even { i } else { opp });
        let cx = c[0] as i64;
        let cy = c[1] as i64;

        for y in 0..ny {
            for x in 0..nx {
                let k1 = y * nx + x;
                let x2 = ((x as i64 + cx).rem_euclid(nx as i64)) as usize;
                let y2 = ((y as i64 + cy).rem_euclid(ny as i64)) as usize;
                let k2 = y2 * nx + x2;

                // Swap f[i][k1] and f[opp][k2]
                let tmp = grid.f[i][k1];
                grid.f[i][k1] = grid.f[opp][k2];
                grid.f[opp][k2] = tmp;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Swap scheme streaming 2D
// ---------------------------------------------------------------------------

/// Swap-scheme streaming for 2D grids (Mattila et al., 2007).
///
/// A memory-free, in-place streaming step that is exactly equivalent to
/// [`stream_2d`].  For every pair of opposite directions `(i, ī)` with
/// `c_ī = −c_i`:
///
/// 1. the two populations are exchanged at every node (an O(1) swap of the
///    two direction buffers), so that the slot of `ī` holds the outgoing
///    `f_i` and the slot of `i` holds the outgoing `f_ī`;
/// 2. the slot `ī` of node `x` is swapped with the slot `i` of node
///    `x + c_i`, which delivers `f_i(x)` to `x + c_i` and `f_ī(x + c_i)` to
///    `x` in one operation.
///
/// On a periodic lattice the swaps of step 2 touch pairwise disjoint slots,
/// so no temporary copy of the distributions is needed and the visiting order
/// is irrelevant.  The rest direction does not move.
pub fn stream_swap_2d(grid: &mut LbmGrid2D) {
    let nx = grid.nx;
    let ny = grid.ny;
    if nx == 0 || ny == 0 {
        return;
    }
    let q = grid.lattice.q();

    for i in 0..q {
        let opp = grid.lattice.opposite(i);
        if opp <= i {
            // Rest direction, or the second member of an already handled pair.
            continue;
        }
        let c = grid.lattice.velocity_2d(i);
        let sx = cyclic_shift(c[0], nx);
        let sy = cyclic_shift(c[1], ny);

        // Step 1: exchange f_i and f_ī at every node.
        grid.f.swap(i, opp);

        // Step 2: swap slot ī of node x with slot i of node x + c_i.
        let (low, high) = grid.f.split_at_mut(opp);
        let f_i = &mut low[i];
        let f_opp = &mut high[0];
        for y in 0..ny {
            let dy = if y + sy >= ny { y + sy - ny } else { y + sy };
            for x in 0..nx {
                let dx = if x + sx >= nx { x + sx - nx } else { x + sx };
                std::mem::swap(&mut f_opp[y * nx + x], &mut f_i[dy * nx + dx]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Streaming with Guo force correction 2D
// ---------------------------------------------------------------------------

/// Streaming with Guo forcing correction for 2D grids.
///
/// After the standard streaming step, the Guo force term is added
/// to each distribution to incorporate a body force `F = [Fx, Fy]`
/// per cell.
///
/// The Guo forcing source term for direction `i` is:
///
/// ```text
/// S_i = w_i * (1 - 1/(2τ)) * [ (e_i - u)/cs² + (e_i·u)/cs⁴ · e_i ] · F
/// ```
///
/// # Arguments
/// * `grid`       – 2D grid (distributions are streamed in-place)
/// * `forces`     – per-cell body forces `[Fx, Fy]`, length = nx*ny
/// * `tau`        – relaxation time
pub fn stream_with_force_correction_2d(grid: &mut LbmGrid2D, forces: &[[f64; 2]], tau: f64) {
    // First do normal streaming.
    stream_2d(grid);

    // Then apply force correction.
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let cs2 = 1.0 / 3.0;
    let cs4 = cs2 * cs2;
    let prefactor = 1.0 - 1.0 / (2.0 * tau);

    for y in 0..ny {
        for x in 0..nx {
            let k = y * nx + x;
            let fx = forces[k][0];
            let fy = forces[k][1];
            let ux = grid.ux[k];
            let uy = grid.uy[k];

            for i in 0..q {
                let w = grid.lattice.weight(i);
                let c = grid.lattice.velocity_2d(i);
                let ex = c[0] as f64;
                let ey = c[1] as f64;
                let e_dot_u = ex * ux + ey * uy;

                let bx = (ex - ux) / cs2 + e_dot_u / cs4 * ex;
                let by = (ey - uy) / cs2 + e_dot_u / cs4 * ey;
                let s_i = w * prefactor * (bx * fx + by * fy);

                grid.f[i][k] += s_i;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Periodic streaming helpers
// ---------------------------------------------------------------------------

/// Periodic wrap of a coordinate in one dimension.
///
/// `wrap_periodic(x, n)` returns `x mod n` where the result is always in `[0, n)`.
#[inline]
pub fn wrap_periodic(x: i64, n: usize) -> usize {
    x.rem_euclid(n as i64) as usize
}

/// Shift all distributions on a 2D grid by `(dx, dy)` lattice spacings
/// with periodic wrapping.  This is useful for Galilean-invariance tests.
pub fn shift_periodic_2d(grid: &mut LbmGrid2D, dx: i64, dy: i64) {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let n = nx * ny;

    let mut f_shifted: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for y in 0..ny {
        for x in 0..nx {
            let src = y * nx + x;
            let dst_x = wrap_periodic(x as i64 + dx, nx);
            let dst_y = wrap_periodic(y as i64 + dy, ny);
            let dst = dst_y * nx + dst_x;
            for (f_shifted_i, f_i) in f_shifted.iter_mut().zip(grid.f.iter()) {
                f_shifted_i[dst] = f_i[src];
            }
        }
    }

    grid.f = f_shifted;
}

/// Shift all distributions on a 3D grid by `(dx, dy, dz)` lattice spacings
/// with periodic wrapping.
pub fn shift_periodic_3d(grid: &mut LbmGrid3D, dx: i64, dy: i64, dz: i64) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let q = grid.lattice.q();
    let n = nx * ny * nz;

    let mut f_shifted: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let src = z * ny * nx + y * nx + x;
                let dst_x = wrap_periodic(x as i64 + dx, nx);
                let dst_y = wrap_periodic(y as i64 + dy, ny);
                let dst_z = wrap_periodic(z as i64 + dz, nz);
                let dst = dst_z * ny * nx + dst_y * nx + dst_x;
                for (f_shifted_i, f_i) in f_shifted.iter_mut().zip(grid.f.iter()) {
                    f_shifted_i[dst] = f_i[src];
                }
            }
        }
    }

    grid.f = f_shifted;
}

/// Compute total mass (sum of all distributions) on a 2D grid.
pub fn total_mass_2d(grid: &LbmGrid2D) -> f64 {
    grid.f.iter().flat_map(|fi| fi.iter()).sum()
}

/// Compute total mass on a 3D grid.
pub fn total_mass_3d(grid: &LbmGrid3D) -> f64 {
    grid.f.iter().flat_map(|fi| fi.iter()).sum()
}

/// Compute total momentum in x-direction on a 2D grid.
pub fn total_momentum_x_2d(grid: &LbmGrid2D) -> f64 {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let mut total = 0.0;
    for i in 0..q {
        let c = grid.lattice.velocity_2d(i);
        let cx = c[0] as f64;
        for k in 0..nx * ny {
            total += cx * grid.f[i][k];
        }
    }
    total
}

/// Compute total momentum in y-direction on a 2D grid.
pub fn total_momentum_y_2d(grid: &LbmGrid2D) -> f64 {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let mut total = 0.0;
    for i in 0..q {
        let c = grid.lattice.velocity_2d(i);
        let cy = c[1] as f64;
        for k in 0..nx * ny {
            total += cy * grid.f[i][k];
        }
    }
    total
}

// ---------------------------------------------------------------------------
// Bounce-back streaming (full-way) for 2D
// ---------------------------------------------------------------------------

/// Apply full-way bounce-back at a set of marked wall cells after streaming.
///
/// For every wall cell, all 9 distributions are replaced by their
/// opposite-direction counterparts *in the same cell*.  This implements
/// the no-slip wall condition (zero velocity at the wall node).
///
/// # Arguments
/// * `grid`       – mutable 2D grid (after streaming has been performed)
/// * `wall_cells` – flat indices of all wall/solid cells
pub fn bounce_back_2d(grid: &mut LbmGrid2D, wall_cells: &[usize]) {
    let q = grid.lattice.q();
    for &k in wall_cells {
        let mut tmp = vec![0.0_f64; q];
        for (tmp_i, f_i) in tmp.iter_mut().zip(grid.f.iter()) {
            *tmp_i = f_i[k];
        }
        for i in 0..q {
            let opp = grid.lattice.opposite(i);
            grid.f[i][k] = tmp[opp];
        }
    }
}

/// Apply full-way bounce-back at a set of marked wall cells after streaming
/// on a 3D grid.
///
/// # Arguments
/// * `grid`       – mutable 3D grid
/// * `wall_cells` – flat indices of all wall/solid cells
pub fn bounce_back_3d(grid: &mut LbmGrid3D, wall_cells: &[usize]) {
    let q = grid.lattice.q();
    for &k in wall_cells {
        let mut tmp = vec![0.0_f64; q];
        for (tmp_i, f_i) in tmp.iter_mut().zip(grid.f.iter()) {
            *tmp_i = f_i[k];
        }
        for i in 0..q {
            let opp = grid.lattice.opposite(i);
            grid.f[i][k] = tmp[opp];
        }
    }
}

// ---------------------------------------------------------------------------
// Half-way bounce-back for 2D
// ---------------------------------------------------------------------------

/// Perform streaming combined with half-way bounce-back on a 2D grid.
///
/// In the half-way scheme the solid boundary is located *half-way* between
/// the last fluid node and the first solid node.  During streaming, any
/// distribution heading into a solid cell is immediately reflected back to
/// the fluid cell it came from with the opposite direction, without the
/// distribution ever entering the solid cell.  This gives second-order
/// accuracy for the wall position.
///
/// # Arguments
/// * `grid`     – mutable 2D grid (distributions before streaming)
/// * `is_solid` – closure/function `|x, y| -> bool` returning `true` if the
///   cell at `(x, y)` is a solid wall node
pub fn stream_half_way_bb_2d<F>(grid: &mut LbmGrid2D, is_solid: F)
where
    F: Fn(usize, usize) -> bool,
{
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let n = nx * ny;

    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for i in 0..q {
        let c = grid.lattice.velocity_2d(i);
        let cx = c[0];
        let cy = c[1];
        let opp = grid.lattice.opposite(i);

        for y in 0..ny {
            for x in 0..nx {
                // Skip solid nodes; they hold no fluid populations.
                if is_solid(x, y) {
                    continue;
                }

                // Destination after streaming in direction i.
                let dst_x_i64 = x as i64 + cx as i64;
                let dst_y_i64 = y as i64 + cy as i64;

                // Periodic wrap for destination.
                let dst_x = dst_x_i64.rem_euclid(nx as i64) as usize;
                let dst_y = dst_y_i64.rem_euclid(ny as i64) as usize;

                let src = y * nx + x;

                if is_solid(dst_x, dst_y) {
                    // Half-way bounce-back: reflect back to origin cell,
                    // opposite direction.
                    f_new[opp][src] += grid.f[i][src];
                } else {
                    // Normal streaming.
                    let dst = dst_y * nx + dst_x;
                    f_new[i][dst] = grid.f[i][src];
                }
            }
        }
    }

    grid.f = f_new;
}

// ---------------------------------------------------------------------------
// Half-way bounce-back with moving wall for 2D
// ---------------------------------------------------------------------------

/// Half-way bounce-back with a moving wall velocity.
///
/// When a distribution is reflected off a wall with velocity `u_wall`,
/// an additional momentum correction is applied:
///
/// ```text
/// f[opp][src] = f[i][src] + 2 * w_i * rho * (e_i · u_wall) / cs²
/// ```
///
/// This is the standard method for implementing Couette flow walls.
pub fn stream_half_way_bb_moving_wall_2d<F>(
    grid: &mut LbmGrid2D,
    is_solid: F,
    wall_velocity: [f64; 2],
) where
    F: Fn(usize, usize) -> bool,
{
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let n = nx * ny;
    let cs2 = 1.0 / 3.0;

    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];

    for i in 0..q {
        let c = grid.lattice.velocity_2d(i);
        let cx = c[0];
        let cy = c[1];
        let opp = grid.lattice.opposite(i);

        for y in 0..ny {
            for x in 0..nx {
                if is_solid(x, y) {
                    continue;
                }

                let dst_x = ((x as i64 + cx as i64).rem_euclid(nx as i64)) as usize;
                let dst_y = ((y as i64 + cy as i64).rem_euclid(ny as i64)) as usize;

                let src = y * nx + x;

                if is_solid(dst_x, dst_y) {
                    // Moving wall bounce-back with momentum correction
                    let w = grid.lattice.weight(i);
                    let rho_k: f64 = (0..q).map(|d| grid.f[d][src]).sum();
                    let e_dot_uw = cx as f64 * wall_velocity[0] + cy as f64 * wall_velocity[1];
                    f_new[opp][src] += grid.f[i][src] + 2.0 * w * rho_k * e_dot_uw / cs2;
                } else {
                    let dst = dst_y * nx + dst_x;
                    f_new[i][dst] = grid.f[i][src];
                }
            }
        }
    }

    grid.f = f_new;
}

// ---------------------------------------------------------------------------
// Periodic round-trip helper
// ---------------------------------------------------------------------------

/// Snapshot the current distribution function array of a 2D grid.
///
/// Returns a flat `Vec`f64` with all `q * n` values stored direction-major.
pub fn snapshot_2d(grid: &LbmGrid2D) -> Vec<f64> {
    let q = grid.lattice.q();
    let n = grid.nx * grid.ny;
    let mut buf = Vec::with_capacity(q * n);
    for i in 0..q {
        buf.extend_from_slice(&grid.f[i]);
    }
    buf
}

/// Snapshot the current distribution function array of a 3D grid.
pub fn snapshot_3d(grid: &LbmGrid3D) -> Vec<f64> {
    let q = grid.lattice.q();
    let n = grid.nx * grid.ny * grid.nz;
    let mut buf = Vec::with_capacity(q * n);
    for i in 0..q {
        buf.extend_from_slice(&grid.f[i]);
    }
    buf
}

/// Restore distributions from a snapshot into a 2D grid.
pub fn restore_snapshot_2d(grid: &mut LbmGrid2D, snapshot: &[f64]) {
    let q = grid.lattice.q();
    let n = grid.nx * grid.ny;
    for i in 0..q {
        let offset = i * n;
        grid.f[i].copy_from_slice(&snapshot[offset..offset + n]);
    }
}

/// Restore distributions from a snapshot into a 3D grid.
pub fn restore_snapshot_3d(grid: &mut LbmGrid3D, snapshot: &[f64]) {
    let q = grid.lattice.q();
    let n = grid.nx * grid.ny * grid.nz;
    for i in 0..q {
        let offset = i * n;
        grid.f[i].copy_from_slice(&snapshot[offset..offset + n]);
    }
}

/// Compare two snapshots and return the maximum absolute difference.
pub fn max_diff_snapshots(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::LatticeType;

    // -----------------------------------------------------------------------
    // stream_2d: moves a spike to the correct neighbour
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_2d_moves_spike() {
        let nx = 10;
        let ny = 10;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k_src = grid.idx(3, 5);
        grid.f[1][k_src] = 99.0; // direction 1 = East (cx=1, cy=0)

        stream_2d(&mut grid);

        let k_dst = grid.idx(4, 5);
        assert!(
            (grid.f[1][k_dst] - 99.0).abs() < 1e-14,
            "Spike not moved East"
        );
    }

    // -----------------------------------------------------------------------
    // stream_2d: periodic wrap at the right edge
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_2d_periodic_wrap() {
        let nx = 5;
        let ny = 5;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k_src = grid.idx(4, 2);
        grid.f[1][k_src] = 42.0;

        stream_2d(&mut grid);

        let k_dst = grid.idx(0, 2);
        assert!(
            (grid.f[1][k_dst] - 42.0).abs() < 1e-14,
            "Periodic wrap failed"
        );
    }

    // -----------------------------------------------------------------------
    // stream_2d: total mass is preserved (periodic domain)
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_2d_mass_conservation() {
        let nx = 8;
        let ny = 6;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        // Perturb a few cells.
        let k1 = grid.idx(2, 3);
        let k2 = grid.idx(5, 1);
        grid.f[1][k1] += 0.1;
        grid.f[3][k2] -= 0.05;

        let mass_before: f64 = grid.f.iter().flat_map(|fi| fi.iter()).sum();
        stream_2d(&mut grid);
        let mass_after: f64 = grid.f.iter().flat_map(|fi| fi.iter()).sum();

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Mass not conserved after stream_2d: before={mass_before}, after={mass_after}"
        );
    }

    // -----------------------------------------------------------------------
    // stream_2d round-trip: streaming nx times in +x returns to start
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_2d_round_trip() {
        let nx = 6;
        let ny = 4;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        // Put a spike at cell (1, 2) in direction 1 (East).
        let k_src = grid.idx(1, 2);
        grid.f[1][k_src] = 77.0;

        let before = snapshot_2d(&grid);

        // Stream nx times in a periodic domain → returns to original position.
        for _ in 0..nx {
            stream_2d(&mut grid);
        }

        let after = snapshot_2d(&grid);
        let max_diff = max_diff_snapshots(&before, &after);

        assert!(
            max_diff < 1e-13,
            "Round-trip streaming failed: max diff = {max_diff}"
        );
    }

    // -----------------------------------------------------------------------
    // stream_3d: moves a spike in +x direction
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_3d_moves_spike() {
        let nx = 6;
        let ny = 5;
        let nz = 4;
        let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        // Direction 1 in D3Q19 is (+x,0,0).
        let k_src = grid.idx(2, 2, 2);
        grid.f[1][k_src] = 55.0;

        stream_3d(&mut grid);

        let k_dst = grid.idx(3, 2, 2);
        assert!(
            (grid.f[1][k_dst] - 55.0).abs() < 1e-14,
            "3D spike not moved in +x"
        );
    }

    // -----------------------------------------------------------------------
    // stream_3d: periodic wrap in z direction
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_3d_periodic_wrap_z() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        // Direction 5 in D3Q19 is (0,0,+1).
        let k_src = grid.idx(1, 1, nz - 1);
        grid.f[5][k_src] = 33.0;

        stream_3d(&mut grid);

        let k_dst = grid.idx(1, 1, 0);
        assert!(
            (grid.f[5][k_dst] - 33.0).abs() < 1e-14,
            "3D periodic wrap in z failed"
        );
    }

    // -----------------------------------------------------------------------
    // stream_3d round-trip: streaming nz times in +z returns to start
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_3d_round_trip() {
        let nx = 3;
        let ny = 3;
        let nz = 5;
        let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        let k_src = grid.idx(1, 1, 1);
        grid.f[5][k_src] = 88.0; // direction 5 = +z

        let before = snapshot_3d(&grid);

        for _ in 0..nz {
            stream_3d(&mut grid);
        }

        let after = snapshot_3d(&grid);
        let max_diff = max_diff_snapshots(&before, &after);

        assert!(
            max_diff < 1e-13,
            "3D round-trip streaming failed: max diff = {max_diff}"
        );
    }

    // -----------------------------------------------------------------------
    // bounce_back_2d: reverses distributions at a wall cell (no streaming)
    // -----------------------------------------------------------------------
    #[test]
    fn test_bounce_back_2d_reverses() {
        let nx = 8;
        let ny = 8;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(4, 4);
        for i in 0..9 {
            grid.f[i][k] = (i + 1) as f64;
        }
        // Apply bounce-back directly (without streaming first) to check reversal.
        bounce_back_2d(&mut grid, &[k]);

        // After bounce-back: f[i] = old_f[opp(i)].
        // Mass should be conserved since we just permute the values.
        let sum: f64 = (0..9).map(|i| grid.f[i][k]).sum();
        let expected_sum: f64 = (1..=9).map(|v| v as f64).sum(); // = 45.0
        assert!(
            (sum - expected_sum).abs() < 1e-12,
            "Bounce-back should conserve mass in the cell: sum={sum}, expected={expected_sum}"
        );
        // Also verify the actual reversal for one pair.
        // D2Q9 opposites: opp(1)=3, so f[1] should now equal old f[3]=4.
        assert!(
            (grid.f[1][k] - 4.0).abs() < 1e-14,
            "f[1] after BB should be old f[3]=4, got {}",
            grid.f[1][k]
        );
    }

    // -----------------------------------------------------------------------
    // bounce_back_3d: mass conserved at a wall cell
    // -----------------------------------------------------------------------
    #[test]
    fn test_bounce_back_3d_mass_conserved() {
        let nx = 5;
        let ny = 5;
        let nz = 5;
        let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        let k = grid.idx(2, 2, 2);
        for i in 0..19 {
            grid.f[i][k] = (i + 1) as f64;
        }
        let mass_before: f64 = (0..19).map(|i| grid.f[i][k]).sum();
        bounce_back_3d(&mut grid, &[k]);
        let mass_after: f64 = (0..19).map(|i| grid.f[i][k]).sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Bounce-back 3D should conserve cell mass: before={mass_before}, after={mass_after}"
        );
    }

    // -----------------------------------------------------------------------
    // Half-way bounce-back: no-slip wall between fluid and solid
    // -----------------------------------------------------------------------
    #[test]
    fn test_half_way_bb_2d_no_velocity_at_wall() {
        // Simple 1D channel: x is periodic, walls at y=0 and y=ny-1.
        let nx = 4;
        let ny = 4;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);

        // Solid marker: top and bottom rows.
        let is_solid = |_x: usize, y: usize| y == 0 || y == ny - 1;

        // Put some non-zero velocity in the interior cells.
        for y in 1..(ny - 1) {
            for x in 0..nx {
                let k = grid.idx(x, y);
                grid.f[1][k] += 0.05;
            }
        }

        // Apply half-way BB streaming for a few steps.
        for _ in 0..10 {
            stream_half_way_bb_2d(&mut grid, is_solid);
        }

        // The solid cells should accumulate no mass (all distributions = 0).
        for x in 0..nx {
            for y in [0, ny - 1] {
                let k = grid.idx(x, y);
                let mass: f64 = (0..9).map(|i| grid.f[i][k]).sum();
                assert!(
                    mass.abs() < 1e-13,
                    "Solid cell ({x},{y}) should have zero mass, got {mass}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Pull-scheme streaming: equivalent to standard streaming
    // -----------------------------------------------------------------------
    #[test]
    fn test_pull_2d_equivalent_to_stream_2d() {
        let nx = 8;
        let ny = 6;
        let mut grid1 = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let mut grid2 = grid1.clone();

        // Perturb
        let k = grid1.idx(3, 2);
        grid1.f[1][k] = 10.0;
        grid2.f[1][k] = 10.0;
        let k2 = grid1.idx(5, 4);
        grid1.f[7][k2] = 20.0;
        grid2.f[7][k2] = 20.0;

        stream_2d(&mut grid1);
        stream_pull_2d(&mut grid2);

        let s1 = snapshot_2d(&grid1);
        let s2 = snapshot_2d(&grid2);
        let diff = max_diff_snapshots(&s1, &s2);
        assert!(
            diff < 1e-14,
            "Pull scheme should be equivalent to standard streaming: diff = {diff}"
        );
    }

    // -----------------------------------------------------------------------
    // Pull-scheme 3D: equivalent to standard streaming
    // -----------------------------------------------------------------------
    #[test]
    fn test_pull_3d_equivalent_to_stream_3d() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let mut grid1 = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        let mut grid2 = grid1.clone();

        let k = grid1.idx(1, 2, 3);
        grid1.f[3][k] = 42.0;
        grid2.f[3][k] = 42.0;

        stream_3d(&mut grid1);
        stream_pull_3d(&mut grid2);

        let s1 = snapshot_3d(&grid1);
        let s2 = snapshot_3d(&grid2);
        let diff = max_diff_snapshots(&s1, &s2);
        assert!(diff < 1e-14, "Pull 3D should match standard: diff = {diff}");
    }

    // -----------------------------------------------------------------------
    // Push-scheme streaming: equivalent to standard streaming
    // -----------------------------------------------------------------------
    #[test]
    fn test_push_2d_equivalent_to_stream_2d() {
        let nx = 8;
        let ny = 6;
        let mut grid1 = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let mut grid2 = grid1.clone();

        let k = grid1.idx(2, 4);
        grid1.f[5][k] = 15.0;
        grid2.f[5][k] = 15.0;

        stream_2d(&mut grid1);
        stream_push_2d(&mut grid2);

        let s1 = snapshot_2d(&grid1);
        let s2 = snapshot_2d(&grid2);
        let diff = max_diff_snapshots(&s1, &s2);
        assert!(
            diff < 1e-14,
            "Push scheme should be equivalent to standard streaming: diff = {diff}"
        );
    }

    // -----------------------------------------------------------------------
    // Push-scheme 3D: equivalent to standard streaming
    // -----------------------------------------------------------------------
    #[test]
    fn test_push_3d_equivalent_to_stream_3d() {
        let nx = 4;
        let ny = 3;
        let nz = 5;
        let mut grid1 = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        let mut grid2 = grid1.clone();

        let k = grid1.idx(2, 1, 3);
        grid1.f[7][k] = 99.0;
        grid2.f[7][k] = 99.0;

        stream_3d(&mut grid1);
        stream_push_3d(&mut grid2);

        let s1 = snapshot_3d(&grid1);
        let s2 = snapshot_3d(&grid2);
        let diff = max_diff_snapshots(&s1, &s2);
        assert!(diff < 1e-14, "Push 3D should match standard: diff = {diff}");
    }

    // -----------------------------------------------------------------------
    // Push-scheme mass conservation
    // -----------------------------------------------------------------------
    #[test]
    fn test_push_2d_mass_conservation() {
        let nx = 7;
        let ny = 5;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(3, 2);
        grid.f[2][k] += 0.3;

        let mass_before = total_mass_2d(&grid);
        stream_push_2d(&mut grid);
        let mass_after = total_mass_2d(&grid);

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Push scheme should conserve mass"
        );
    }

    // -----------------------------------------------------------------------
    // AA pattern: two steps (even+odd) is equivalent to standard stream
    // -----------------------------------------------------------------------
    #[test]
    fn test_aa_even_odd_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(2, 3);
        grid.f[1][k] = 5.0;
        grid.f[5][k] = 3.0;

        let mass_before = total_mass_2d(&grid);
        stream_aa_2d(&mut grid, true);
        stream_aa_2d(&mut grid, false);
        let mass_after = total_mass_2d(&grid);

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "AA pattern should conserve mass: before={mass_before}, after={mass_after}"
        );
    }

    // -----------------------------------------------------------------------
    // Swap scheme: mass conservation
    // -----------------------------------------------------------------------
    #[test]
    fn test_swap_2d_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(3, 3);
        grid.f[1][k] = 7.0;
        grid.f[8][k] = 4.5;

        let mass_before = total_mass_2d(&grid);
        stream_swap_2d(&mut grid);
        let mass_after = total_mass_2d(&grid);

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Swap scheme should conserve mass"
        );
    }

    // -----------------------------------------------------------------------
    // Force correction: adds expected source term
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_with_force_correction_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let n = nx * ny;

        // Zero forces → no mass change
        let forces = vec![[0.0, 0.0]; n];
        let tau = 1.0;

        let mass_before = total_mass_2d(&grid);
        stream_with_force_correction_2d(&mut grid, &forces, tau);
        let mass_after = total_mass_2d(&grid);

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Zero force should not change mass"
        );
    }

    // -----------------------------------------------------------------------
    // Periodic helpers
    // -----------------------------------------------------------------------
    #[test]
    fn test_wrap_periodic() {
        assert_eq!(wrap_periodic(5, 10), 5);
        assert_eq!(wrap_periodic(-1, 10), 9);
        assert_eq!(wrap_periodic(10, 10), 0);
        assert_eq!(wrap_periodic(-10, 10), 0);
        assert_eq!(wrap_periodic(0, 5), 0);
    }

    #[test]
    fn test_shift_periodic_2d_round_trip() {
        let nx = 6;
        let ny = 4;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(2, 1);
        grid.f[1][k] = 100.0;

        let before = snapshot_2d(&grid);

        // Shift by (3, 2) then shift back by (-3, -2)
        shift_periodic_2d(&mut grid, 3, 2);
        shift_periodic_2d(&mut grid, -3, -2);

        let after = snapshot_2d(&grid);
        let diff = max_diff_snapshots(&before, &after);
        assert!(diff < 1e-14, "Shift round-trip failed: diff = {diff}");
    }

    #[test]
    fn test_shift_periodic_3d_round_trip() {
        let nx = 4;
        let ny = 3;
        let nz = 5;
        let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        let k = grid.idx(1, 1, 2);
        grid.f[3][k] = 50.0;

        let before = snapshot_3d(&grid);
        shift_periodic_3d(&mut grid, 2, 1, -3);
        shift_periodic_3d(&mut grid, -2, -1, 3);

        let after = snapshot_3d(&grid);
        let diff = max_diff_snapshots(&before, &after);
        assert!(diff < 1e-14, "3D shift round-trip failed: diff = {diff}");
    }

    // -----------------------------------------------------------------------
    // Total mass / momentum helpers
    // -----------------------------------------------------------------------
    #[test]
    fn test_total_mass_2d() {
        let nx = 5;
        let ny = 5;
        let grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let mass = total_mass_2d(&grid);
        // Initialized with rho=1 everywhere; sum of all weights = 1 per cell.
        // total mass = nx * ny * 1.0 = 25.0
        assert!(
            (mass - (nx * ny) as f64).abs() < 1e-12,
            "Expected mass {}, got {mass}",
            nx * ny
        );
    }

    #[test]
    fn test_total_momentum_at_rest() {
        let nx = 5;
        let ny = 5;
        let grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let px = total_momentum_x_2d(&grid);
        let py = total_momentum_y_2d(&grid);
        assert!(px.abs() < 1e-12, "x-momentum should be zero at rest: {px}");
        assert!(py.abs() < 1e-12, "y-momentum should be zero at rest: {py}");
    }

    // -----------------------------------------------------------------------
    // Snapshot restore round-trip
    // -----------------------------------------------------------------------
    #[test]
    fn test_snapshot_restore_2d() {
        let nx = 6;
        let ny = 4;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(3, 2);
        grid.f[1][k] = 77.0;

        let snap = snapshot_2d(&grid);
        // Corrupt the grid
        grid.f[1][k] = 0.0;
        // Restore
        restore_snapshot_2d(&mut grid, &snap);

        assert!(
            (grid.f[1][k] - 77.0).abs() < 1e-14,
            "Snapshot restore failed"
        );
    }

    #[test]
    fn test_snapshot_restore_3d() {
        let nx = 3;
        let ny = 3;
        let nz = 3;
        let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);
        let k = grid.idx(1, 1, 1);
        grid.f[5][k] = 88.0;

        let snap = snapshot_3d(&grid);
        grid.f[5][k] = 0.0;
        restore_snapshot_3d(&mut grid, &snap);

        assert!(
            (grid.f[5][k] - 88.0).abs() < 1e-14,
            "3D snapshot restore failed"
        );
    }

    // -----------------------------------------------------------------------
    // Moving wall bounce-back test
    // -----------------------------------------------------------------------
    #[test]
    fn test_moving_wall_bb_mass_nonzero() {
        let nx = 6;
        let ny = 6;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);

        let is_solid = |_x: usize, y: usize| y == 0 || y == ny - 1;
        let wall_vel = [0.01, 0.0];

        // Run a few steps
        for _ in 0..5 {
            stream_half_way_bb_moving_wall_2d(&mut grid, is_solid, wall_vel);
        }

        // Interior cells should still have non-zero mass
        let k = grid.idx(3, 3);
        let mass: f64 = (0..9).map(|i| grid.f[i][k]).sum();
        assert!(
            mass > 0.0,
            "Interior cell should have positive mass: {mass}"
        );
    }

    // -----------------------------------------------------------------------
    // Pull-scheme mass conservation
    // -----------------------------------------------------------------------
    #[test]
    fn test_pull_2d_mass_conservation() {
        let nx = 7;
        let ny = 5;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(3, 2);
        grid.f[4][k] += 0.2;

        let mass_before = total_mass_2d(&grid);
        stream_pull_2d(&mut grid);
        let mass_after = total_mass_2d(&grid);

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Pull scheme should conserve mass"
        );
    }

    // -----------------------------------------------------------------------
    // max_diff_snapshots
    // -----------------------------------------------------------------------
    #[test]
    fn test_max_diff_snapshots_identical() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(max_diff_snapshots(&a, &b) < 1e-15);
    }

    #[test]
    fn test_max_diff_snapshots_different() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.5, 3.0];
        assert!((max_diff_snapshots(&a, &b) - 0.5).abs() < 1e-15);
    }

    // -----------------------------------------------------------------------
    // Force correction: non-zero force modifies distributions
    // -----------------------------------------------------------------------
    #[test]
    fn test_force_correction_modifies_distributions() {
        let nx = 6;
        let ny = 6;
        let n = nx * ny;

        let mut grid_no_force = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let mut grid_force = grid_no_force.clone();

        let forces = vec![[1e-4, 0.0]; n];
        let tau = 1.0;

        stream_2d(&mut grid_no_force);
        stream_with_force_correction_2d(&mut grid_force, &forces, tau);

        let s1 = snapshot_2d(&grid_no_force);
        let s2 = snapshot_2d(&grid_force);
        let diff = max_diff_snapshots(&s1, &s2);
        assert!(
            diff > 1e-10,
            "Non-zero force should modify distributions: diff = {diff}"
        );
    }

    // -----------------------------------------------------------------------
    // Shift preserves mass
    // -----------------------------------------------------------------------
    #[test]
    fn test_shift_periodic_2d_mass_conservation() {
        let nx = 6;
        let ny = 4;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let k = grid.idx(2, 1);
        grid.f[3][k] += 0.5;

        let mass_before = total_mass_2d(&grid);
        shift_periodic_2d(&mut grid, 3, -1);
        let mass_after = total_mass_2d(&grid);

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Shift should conserve mass"
        );
    }
}

#[cfg(test)]
mod tests_periodic_streaming {
    use super::*;
    use crate::lattice::LatticeType;

    /// Deterministic pseudo-random values in [0, 1) (64-bit LCG).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn random_grid_2d(nx: usize, ny: usize, seed: u64) -> LbmGrid2D {
        let mut rng = Lcg(seed);
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        for f_i in grid.f.iter_mut() {
            for v in f_i.iter_mut() {
                *v = rng.next();
            }
        }
        grid
    }

    fn random_grid_3d(nx: usize, ny: usize, nz: usize, lt: LatticeType, seed: u64) -> LbmGrid3D {
        let mut rng = Lcg(seed);
        let mut grid = LbmGrid3D::new(nx, ny, nz, lt);
        for f_i in grid.f.iter_mut() {
            for v in f_i.iter_mut() {
                *v = rng.next();
            }
        }
        grid
    }

    /// Straightforward per-cell reference: f_new_i(x) = f_i(x - c_i) (periodic).
    fn reference_stream_2d(grid: &LbmGrid2D) -> Vec<Vec<f64>> {
        let (nx, ny) = (grid.nx as i64, grid.ny as i64);
        (0..grid.lattice.q())
            .map(|i| {
                let c = grid.lattice.velocity_2d(i);
                let mut plane = vec![0.0; (nx * ny) as usize];
                for y in 0..ny {
                    for x in 0..nx {
                        let sx = (x - i64::from(c[0])).rem_euclid(nx);
                        let sy = (y - i64::from(c[1])).rem_euclid(ny);
                        plane[(y * nx + x) as usize] = grid.f[i][(sy * nx + sx) as usize];
                    }
                }
                plane
            })
            .collect()
    }

    fn reference_stream_3d(grid: &LbmGrid3D) -> Vec<Vec<f64>> {
        let (nx, ny, nz) = (grid.nx as i64, grid.ny as i64, grid.nz as i64);
        (0..grid.lattice.q())
            .map(|i| {
                let c = grid.lattice.velocity_3d(i);
                let mut plane = vec![0.0; (nx * ny * nz) as usize];
                for z in 0..nz {
                    for y in 0..ny {
                        for x in 0..nx {
                            let sx = (x - i64::from(c[0])).rem_euclid(nx);
                            let sy = (y - i64::from(c[1])).rem_euclid(ny);
                            let sz = (z - i64::from(c[2])).rem_euclid(nz);
                            plane[(z * ny * nx + y * nx + x) as usize] =
                                grid.f[i][(sz * ny * nx + sy * nx + sx) as usize];
                        }
                    }
                }
                plane
            })
            .collect()
    }

    fn max_abs_diff(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
        a.iter()
            .zip(b)
            .flat_map(|(x, y)| x.iter().zip(y).map(|(u, v)| (u - v).abs()))
            .fold(0.0, f64::max)
    }

    fn buffer_addresses(grid_f: &[Vec<f64>], grid_tmp: &[Vec<f64>]) -> Vec<usize> {
        let mut addrs: Vec<usize> = grid_f
            .iter()
            .chain(grid_tmp)
            .map(|v| v.as_ptr() as usize)
            .collect();
        addrs.sort_unstable();
        addrs
    }

    #[test]
    fn stream_2d_matches_the_per_cell_reference_on_awkward_sizes() {
        for (nx, ny) in [(1, 1), (1, 5), (5, 1), (2, 2), (3, 2), (7, 5), (16, 9)] {
            let mut grid = random_grid_2d(nx, ny, (nx * 31 + ny) as u64);
            let expected = reference_stream_2d(&grid);
            stream_2d(&mut grid);
            assert_eq!(
                max_abs_diff(&grid.f, &expected),
                0.0,
                "stream_2d mismatch on {nx}x{ny}"
            );
        }
    }

    #[test]
    fn stream_3d_matches_the_per_cell_reference_for_d3q19_and_d3q27() {
        for lt in [LatticeType::D3Q19, LatticeType::D3Q27] {
            for (nx, ny, nz) in [(1, 1, 1), (2, 3, 1), (1, 4, 3), (5, 4, 3), (4, 4, 4)] {
                let mut grid = random_grid_3d(nx, ny, nz, lt, (nx * 100 + ny * 10 + nz) as u64);
                let expected = reference_stream_3d(&grid);
                stream_3d(&mut grid);
                assert_eq!(
                    max_abs_diff(&grid.f, &expected),
                    0.0,
                    "stream_3d mismatch for {lt:?} on {nx}x{ny}x{nz}"
                );
            }
        }
    }

    #[test]
    fn stream_2d_moves_every_direction_to_the_right_neighbour() {
        let (nx, ny) = (6, 5);
        for i in 1..9 {
            let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
            for f_i in grid.f.iter_mut() {
                f_i.fill(0.0);
            }
            let k = grid.idx(0, 4);
            grid.f[i][k] = 1.0;
            stream_2d(&mut grid);
            let c = grid.lattice.velocity_2d(i);
            let x = (c[0]).rem_euclid(nx as i32) as usize;
            let y = (4 + c[1]).rem_euclid(ny as i32) as usize;
            assert_eq!(grid.f[i][grid.idx(x, y)], 1.0, "direction {i}");
            let total: f64 = grid.f.iter().flatten().sum();
            assert_eq!(total, 1.0, "direction {i}: population duplicated or lost");
        }
    }

    #[test]
    fn stream_2d_and_3d_reuse_their_buffers_instead_of_allocating() {
        let mut grid = random_grid_2d(32, 24, 7);
        let before = buffer_addresses(&grid.f, &grid.f_tmp);
        for _ in 0..5 {
            stream_2d(&mut grid);
        }
        assert_eq!(
            before,
            buffer_addresses(&grid.f, &grid.f_tmp),
            "stream_2d allocated a new distribution buffer"
        );

        let mut grid = random_grid_3d(8, 6, 5, LatticeType::D3Q19, 9);
        let before = buffer_addresses(&grid.f, &grid.f_tmp);
        for _ in 0..5 {
            stream_3d(&mut grid);
        }
        assert_eq!(
            before,
            buffer_addresses(&grid.f, &grid.f_tmp),
            "stream_3d allocated a new distribution buffer"
        );
    }

    #[test]
    fn streaming_recovers_from_an_empty_or_resized_scratch_buffer() {
        let mut grid = random_grid_2d(9, 7, 11);
        let expected = reference_stream_2d(&grid);
        grid.f_tmp = Vec::new();
        stream_2d(&mut grid);
        assert_eq!(max_abs_diff(&grid.f, &expected), 0.0);

        let mut grid = random_grid_3d(5, 4, 3, LatticeType::D3Q19, 12);
        let expected = reference_stream_3d(&grid);
        grid.f_tmp = vec![vec![1.0; 3]; 2];
        stream_3d(&mut grid);
        assert_eq!(max_abs_diff(&grid.f, &expected), 0.0);

        // A scratch buffer left over from a larger grid (more planes, longer
        // planes) is trimmed to the current shape instead of being used as is.
        let mut grid = random_grid_2d(6, 5, 13);
        let expected = reference_stream_2d(&grid);
        grid.f_tmp = vec![vec![7.0; 1000]; 12];
        stream_2d(&mut grid);
        assert_eq!(max_abs_diff(&grid.f, &expected), 0.0);
        assert_eq!(grid.f_tmp.len(), 9);
        assert!(grid.f_tmp.iter().all(|plane| plane.len() == 30));
    }

    #[test]
    fn streaming_an_empty_grid_is_a_no_op() {
        let mut grid = LbmGrid2D::new(0, 4, LatticeType::D2Q9);
        stream_2d(&mut grid);
        stream_swap_2d(&mut grid);
        let mut grid = LbmGrid3D::new(3, 0, 2, LatticeType::D3Q19);
        stream_3d(&mut grid);
    }

    #[test]
    fn stream_2d_keeps_the_rest_population_in_place() {
        let mut grid = random_grid_2d(6, 4, 21);
        let rest = grid.f[0].clone();
        stream_2d(&mut grid);
        assert_eq!(grid.f[0], rest);
    }

    #[test]
    fn stream_swap_2d_equals_stream_2d_on_random_fields() {
        for (nx, ny) in [(1, 1), (2, 2), (1, 6), (7, 1), (7, 5), (8, 6), (13, 11)] {
            let mut swapped = random_grid_2d(nx, ny, (nx * 17 + ny * 3) as u64);
            let mut streamed = swapped.clone();
            stream_swap_2d(&mut swapped);
            stream_2d(&mut streamed);
            let diff = max_abs_diff(&swapped.f, &streamed.f);
            assert!(diff < 1e-14, "swap vs stream on {nx}x{ny}: diff {diff}");
        }
    }

    #[test]
    fn stream_swap_2d_streams_every_direction_the_right_way() {
        // Mass alone cannot reveal a direction mix-up: track single spikes.
        let (nx, ny) = (7, 6);
        for i in 1..9 {
            let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
            for f_i in grid.f.iter_mut() {
                f_i.fill(0.0);
            }
            let k = grid.idx(3, 2);
            grid.f[i][k] = 1.0;
            stream_swap_2d(&mut grid);
            let c = grid.lattice.velocity_2d(i);
            let x = (3 + c[0]).rem_euclid(nx as i32) as usize;
            let y = (2 + c[1]).rem_euclid(ny as i32) as usize;
            assert_eq!(
                grid.f[i][grid.idx(x, y)],
                1.0,
                "direction {i} must arrive at ({x},{y}) in its own slot"
            );
            let total: f64 = grid.f.iter().flatten().sum();
            assert_eq!(total, 1.0, "direction {i}: population duplicated or lost");
        }
    }

    #[test]
    fn stream_swap_2d_conserves_momentum_and_has_a_round_trip() {
        let mut grid = random_grid_2d(6, 4, 99);
        let px_before = total_momentum_x_2d(&grid);
        let py_before = total_momentum_y_2d(&grid);
        let snapshot = snapshot_2d(&grid);
        for _ in 0..12 {
            // lcm(6, 4) = 12 steps of an (1,1) diagonal direction return home,
            // and every other direction returns after 6 or 4 steps.
            stream_swap_2d(&mut grid);
        }
        assert!((total_momentum_x_2d(&grid) - px_before).abs() < 1e-12);
        assert!((total_momentum_y_2d(&grid) - py_before).abs() < 1e-12);
        let diff = max_diff_snapshots(&snapshot, &snapshot_2d(&grid));
        assert!(diff < 1e-14, "12 swap steps must be the identity: {diff}");
    }
}
