//! Streaming, equilibrium, and macro utility functions for LBM lattices.

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    CS2, D2Q9_OPPOSITES, D2Q9_VELOCITIES, D2Q9_WEIGHTS, D3Q19_OPPOSITES, D3Q19_VELOCITIES,
    D3Q19_WEIGHTS,
};

/// Perform a full periodic streaming step for a D2Q9 population array.
///
/// `f` is stored as `f[y * nx + x][q]` (cell-major).
/// After streaming, each population `q` moves from cell `(x,y)` to
/// `(x + cx_q, y + cy_q)` with periodic wrapping.
#[allow(dead_code)]
pub fn stream_d2q9_periodic(f: &mut Vec<[f64; 9]>, nx: usize, ny: usize) {
    let cx: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
    let cy: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
    let src = f.clone();
    for y in 0..ny {
        for x in 0..nx {
            for q in 0..9usize {
                let src_x = ((x as i32 - cx[q]).rem_euclid(nx as i32)) as usize;
                let src_y = ((y as i32 - cy[q]).rem_euclid(ny as i32)) as usize;
                f[y * nx + x][q] = src[src_y * nx + src_x][q];
            }
        }
    }
}
/// Perform a full periodic streaming step for a D3Q19 population array.
///
/// `f` is stored as `f[z * nx*ny + y * nx + x][q]` (cell-major).
#[allow(dead_code)]
pub fn stream_d3q19_periodic(f: &mut Vec<[f64; 19]>, nx: usize, ny: usize, nz: usize) {
    let (cx, cy, cz) = d3q19_vel_arrays();
    let src = f.clone();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let dst = z * nx * ny + y * nx + x;
                for q in 0..19usize {
                    let sx = ((x as i32 - cx[q]).rem_euclid(nx as i32)) as usize;
                    let sy = ((y as i32 - cy[q]).rem_euclid(ny as i32)) as usize;
                    let sz = ((z as i32 - cz[q]).rem_euclid(nz as i32)) as usize;
                    f[dst][q] = src[sz * nx * ny + sy * nx + sx][q];
                }
            }
        }
    }
}
/// Helper: return separate cx/cy/cz arrays for D3Q19.
#[inline]
#[allow(dead_code)]
pub fn d3q19_vel_arrays() -> ([i32; 19], [i32; 19], [i32; 19]) {
    let mut cx = [0i32; 19];
    let mut cy = [0i32; 19];
    let mut cz = [0i32; 19];
    for (q, v) in D3Q19_VELOCITIES.iter().enumerate() {
        cx[q] = v[0];
        cy[q] = v[1];
        cz[q] = v[2];
    }
    (cx, cy, cz)
}
/// Compute the full D2Q9 equilibrium distribution for given `rho`, `ux`, `uy`.
///
/// Returns `[f64; 9]` with `feq_i = w_i * rho * (1 + eu/cs² + eu²/(2cs⁴) − u²/(2cs²))`.
#[allow(dead_code)]
pub fn equilibrium_d2q9(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let mut feq = [0.0f64; 9];
    let u2 = ux * ux + uy * uy;
    for q in 0..9 {
        let cx = D2Q9_VELOCITIES[q][0] as f64;
        let cy = D2Q9_VELOCITIES[q][1] as f64;
        let eu = cx * ux + cy * uy;
        feq[q] = D2Q9_WEIGHTS[q]
            * rho
            * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}
/// Compute the full D3Q19 equilibrium distribution for given `rho`, `ux`, `uy`, `uz`.
#[allow(dead_code)]
pub fn equilibrium_d3q19(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; 19] {
    let mut feq = [0.0f64; 19];
    let u2 = ux * ux + uy * uy + uz * uz;
    for q in 0..19 {
        let cx = D3Q19_VELOCITIES[q][0] as f64;
        let cy = D3Q19_VELOCITIES[q][1] as f64;
        let cz = D3Q19_VELOCITIES[q][2] as f64;
        let eu = cx * ux + cy * uy + cz * uz;
        feq[q] = D3Q19_WEIGHTS[q]
            * rho
            * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}
/// Compute density `rho` and velocity `(ux, uy)` from a D2Q9 population `[f64; 9]`.
#[allow(dead_code)]
pub fn macros_from_d2q9(f: &[f64; 9]) -> (f64, f64, f64) {
    let rho: f64 = f.iter().sum();
    if rho < 1e-30 {
        return (0.0, 0.0, 0.0);
    }
    let mut ux = 0.0f64;
    let mut uy = 0.0f64;
    for q in 0..9 {
        ux += D2Q9_VELOCITIES[q][0] as f64 * f[q];
        uy += D2Q9_VELOCITIES[q][1] as f64 * f[q];
    }
    (rho, ux / rho, uy / rho)
}
/// Compute density and velocity from a D3Q19 population `[f64; 19]`.
#[allow(dead_code)]
pub fn macros_from_d3q19(f: &[f64; 19]) -> (f64, f64, f64, f64) {
    let rho: f64 = f.iter().sum();
    if rho < 1e-30 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut ux = 0.0f64;
    let mut uy = 0.0f64;
    let mut uz = 0.0f64;
    for q in 0..19 {
        ux += D3Q19_VELOCITIES[q][0] as f64 * f[q];
        uy += D3Q19_VELOCITIES[q][1] as f64 * f[q];
        uz += D3Q19_VELOCITIES[q][2] as f64 * f[q];
    }
    (rho, ux / rho, uy / rho, uz / rho)
}
/// Apply BGK collision to a D2Q9 node in-place.
///
/// `omega = 1/tau` is the relaxation frequency.
#[allow(dead_code)]
pub fn bgk_d2q9(f: &mut [f64; 9], rho: f64, ux: f64, uy: f64, omega: f64) {
    let feq = equilibrium_d2q9(rho, ux, uy);
    for q in 0..9 {
        f[q] += omega * (feq[q] - f[q]);
    }
}
/// Apply BGK collision to a D3Q19 node in-place.
#[allow(dead_code)]
pub fn bgk_d3q19(f: &mut [f64; 19], rho: f64, ux: f64, uy: f64, uz: f64, omega: f64) {
    let feq = equilibrium_d3q19(rho, ux, uy, uz);
    for q in 0..19 {
        f[q] += omega * (feq[q] - f[q]);
    }
}
/// Apply half-way bounce-back at a D2Q9 node (no-slip wall).
/// Swaps each population with its opposite direction.
#[allow(dead_code)]
pub fn bounce_back_node_d2q9(f: &mut [f64; 9]) {
    for q in 0..9 {
        let opp = D2Q9_OPPOSITES[q];
        if opp > q {
            f.swap(q, opp);
        }
    }
}
/// Apply half-way bounce-back at a D3Q19 node.
#[allow(dead_code)]
pub fn bounce_back_node_d3q19(f: &mut [f64; 19]) {
    for q in 0..19 {
        let opp = D3Q19_OPPOSITES[q];
        if opp > q {
            f.swap(q, opp);
        }
    }
}
/// Compute the non-equilibrium distribution `fneq_i = f_i - feq_i` for D2Q9.
#[allow(dead_code)]
pub fn non_equilibrium_d2q9(f: &[f64; 9], feq: &[f64; 9]) -> [f64; 9] {
    let mut fneq = [0.0f64; 9];
    for q in 0..9 {
        fneq[q] = f[q] - feq[q];
    }
    fneq
}
