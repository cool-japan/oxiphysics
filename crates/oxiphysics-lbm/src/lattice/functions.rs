//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::D2Q9;
#[cfg(test)]
use super::types::{
    D2Q9Grid, D3Q19, D3Q27, GuoBodyForce, Lattice, LatticeD2Q9Grid, LatticeDimensions, LatticeType,
    MrtCollision, SmagorinskyModel, TrtCollision,
};

/// D2Q9 weights.
pub const D2Q9_WEIGHTS: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// D2Q9 velocity vectors `[cx, cy]`.
///
/// Order: rest, E, N, W, S, NE, NW, SW, SE.
pub const D2Q9_VELOCITIES: [[i32; 2]; 9] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
    [1, 1],
    [-1, 1],
    [-1, -1],
    [1, -1],
];
/// Opposite direction indices for D2Q9.
pub const D2Q9_OPPOSITES: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
/// D3Q19 weights.
pub const D3Q19_WEIGHTS: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// D3Q19 velocity vectors `[cx, cy, cz]`.
pub const D3Q19_VELOCITIES: [[i32; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, 1],
    [1, 0, -1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, 1],
    [0, 1, -1],
    [0, -1, -1],
];
/// Opposite direction indices for D3Q19.
pub const D3Q19_OPPOSITES: [usize; 19] = [
    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15,
];
/// D3Q27 weights.
pub const D3Q27_WEIGHTS: [f64; 27] = [
    8.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
];
/// D3Q27 velocity vectors `[cx, cy, cz]`.
pub const D3Q27_VELOCITIES: [[i32; 3]; 27] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, 1],
    [1, 0, -1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, 1],
    [0, 1, -1],
    [0, -1, -1],
    [1, 1, 1],
    [-1, 1, 1],
    [1, -1, 1],
    [-1, -1, 1],
    [1, 1, -1],
    [-1, 1, -1],
    [1, -1, -1],
    [-1, -1, -1],
];
/// Opposite direction indices for D3Q27.
pub const D3Q27_OPPOSITES: [usize; 27] = [
    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20,
    19,
];
/// Speed of sound squared in lattice units: cs² = 1/3.
pub const CS2: f64 = 1.0 / 3.0;
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_lattice_dimensions_total_cells() {
        let dims = LatticeDimensions::new(4, 5, 6);
        assert_eq!(dims.total_cells(), 120);
    }
    #[test]
    fn test_lattice_dimensions_idx_coords_roundtrip() {
        let dims = LatticeDimensions::new(7, 3, 5);
        for z in 0..dims.nz {
            for y in 0..dims.ny {
                for x in 0..dims.nx {
                    let idx = dims.idx_3d(x, y, z);
                    let (rx, ry, rz) = dims.coords(idx);
                    assert_eq!((x, y, z), (rx, ry, rz));
                }
            }
        }
    }
    #[test]
    fn test_d2q9grid_feq_sums_to_rho() {
        let rho = 1.7;
        let ux = 0.08;
        let uy = -0.04;
        let sum: f64 = (0..9).map(|a| D2Q9Grid::feq(rho, ux, uy, a)).sum();
        assert!((sum - rho).abs() < 1e-14, "feq sum = {sum}, expected {rho}");
    }
    #[test]
    fn test_d2q9grid_feq_zero_velocity() {
        let rho = 2.0;
        for (alpha, &w) in D2Q9_WEIGHTS.iter().enumerate() {
            let feq = D2Q9Grid::feq(rho, 0.0, 0.0, alpha);
            let expected = w * rho;
            assert!(
                (feq - expected).abs() < 1e-14,
                "feq[{alpha}] = {feq}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_d2q9grid_stream_conserves_mass() {
        let nx = 8;
        let ny = 6;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        let base = (3 * nx + 2) * 9;
        grid.f[base] += 0.5;
        grid.f[base + 1] -= 0.2;
        let mass_before = grid.total_mass();
        grid.stream();
        let mass_after = grid.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Stream did not conserve mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d2q9grid_equilibrium_unchanged_by_collision() {
        let nx = 4;
        let ny = 4;
        let rho0 = 1.0;
        let grid_before = D2Q9Grid::new(nx, ny, rho0);
        let mut grid = grid_before.clone();
        grid.collide_bgk(1.5);
        for i in 0..grid.f.len() {
            assert!(
                (grid.f[i] - grid_before.f[i]).abs() < 1e-14,
                "Equilibrium changed by collision at index {i}"
            );
        }
    }
    #[test]
    fn test_d2q9grid_collision_conserves_mass() {
        let nx = 6;
        let ny = 6;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        grid.f[0] += 0.3;
        grid.f[1] -= 0.15;
        grid.f[2] -= 0.15;
        let mass_before = grid.total_mass();
        grid.collide_bgk(1.2);
        let mass_after = grid.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Collision did not conserve mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d2q9grid_density_velocity_fields() {
        let nx = 5;
        let ny = 5;
        let rho0 = 1.3;
        let grid = D2Q9Grid::new(nx, ny, rho0);
        let rho = grid.density_field();
        let vel = grid.velocity_field();
        for idx in 0..(nx * ny) {
            assert!(
                (rho[idx] - rho0).abs() < 1e-14,
                "density[{idx}] = {}, expected {rho0}",
                rho[idx]
            );
            assert!(
                vel[idx][0].abs() < 1e-14 && vel[idx][1].abs() < 1e-14,
                "velocity[{idx}] should be zero"
            );
        }
    }
    #[test]
    fn test_d2q9grid_multi_step_conservation() {
        let nx = 10;
        let ny = 10;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        for alpha in 0..9 {
            grid.f[(3 * nx + 5) * 9 + alpha] = D2Q9Grid::feq(1.1, 0.05, 0.02, alpha);
        }
        let mass_before = grid.total_mass();
        for _ in 0..50 {
            grid.collide_bgk(1.0);
            grid.stream();
        }
        let mass_after = grid.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Multi-step mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
}
/// D2Q9 MRT transformation matrix M (row vectors are the 9 modes).
///
/// Based on Lallemand & Luo (2000).
pub const MRT_M: [[f64; 9]; 9] = [
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0],
    [4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0],
    [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],
    [0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0],
    [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0],
    [0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0],
    [0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0],
];
#[cfg(test)]
mod expanded_tests {
    use super::*;
    #[test]
    fn test_trt_collision_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        grid.f[0] += 0.3;
        grid.f[1] -= 0.15;
        grid.f[2] -= 0.15;
        let mass_before = grid.total_mass();
        let trt = TrtCollision::new(1.0, 1.0);
        trt.collide(&mut grid.f, nx * ny);
        let mass_after = grid.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "TRT mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_trt_viscosity() {
        let trt = TrtCollision::new(1.0, 1.0);
        let nu = trt.viscosity();
        assert!(nu > 0.0, "viscosity should be positive");
    }
    #[test]
    fn test_trt_magic_parameter() {
        let trt = TrtCollision::with_magic(1.0);
        assert!(
            trt.omega_minus > 0.0 && trt.omega_minus <= 2.0,
            "omega_minus = {}",
            trt.omega_minus
        );
    }
    #[test]
    fn test_mrt_moment_roundtrip() {
        let f = [0.4, 0.1, 0.1, 0.1, 0.1, 0.036, 0.036, 0.036, 0.036];
        let m = MrtCollision::to_moment_space(&f);
        let f2 = MrtCollision::from_moment_space(&m);
        for i in 0..9 {
            assert!(
                (f[i] - f2[i]).abs() < 1e-10,
                "roundtrip error at {i}: {} vs {}",
                f[i],
                f2[i]
            );
        }
    }
    #[test]
    fn test_mrt_collision_mass_conservation() {
        let nx = 4;
        let ny = 4;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        grid.f[0] += 0.1;
        let mass_before = grid.total_mass();
        let mrt = MrtCollision::uniform(1.0);
        mrt.collide(&mut grid.f, nx * ny);
        let mass_after = grid.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-11,
            "MRT mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_mrt_equilibrium_unchanged() {
        let nx = 4;
        let ny = 4;
        let grid_before = D2Q9Grid::new(nx, ny, 1.0);
        let mut grid = grid_before.clone();
        let mrt = MrtCollision::uniform(1.0);
        mrt.collide(&mut grid.f, nx * ny);
        for i in 0..grid.f.len() {
            assert!(
                (grid.f[i] - grid_before.f[i]).abs() < 1e-13,
                "Equilibrium changed at {i}"
            );
        }
    }
    #[test]
    fn test_smagorinsky_effective_omega_positive() {
        let model = SmagorinskyModel::new(0.18, 1.0, 1.0 / 6.0);
        let nx = 4;
        let ny = 4;
        let grid = D2Q9Grid::new(nx, ny, 1.0);
        let omega = model.effective_omega(&grid.f, 0);
        assert!(omega > 0.0, "effective omega should be positive");
    }
    #[test]
    fn test_guo_forcing_term_zero_force() {
        let guo = GuoBodyForce::new(0.0, 0.0);
        for a in 0..9 {
            let ft = guo.forcing_term(0.05, 0.0, a);
            assert!(
                ft.abs() < 1e-60,
                "forcing term with zero g should be zero at {a}: {ft}"
            );
        }
    }
    #[test]
    fn test_guo_collide_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        let mass_before = grid.total_mass();
        let guo = GuoBodyForce::new(1e-5, 0.0);
        guo.collide_with_force(&mut grid.f, nx * ny, 1.0);
        let mass_after = grid.total_mass();
        assert!(
            (mass_after - mass_before).abs() < 0.01 * mass_before,
            "mass change with body force: {:.2e}",
            (mass_after - mass_before).abs()
        );
    }
    #[test]
    fn test_guo_collide_velocity_increase() {
        let nx = 4;
        let ny = 4;
        let mut grid = D2Q9Grid::new(nx, ny, 1.0);
        let guo = GuoBodyForce::new(1e-4, 0.0);
        for _ in 0..10 {
            guo.collide_with_force(&mut grid.f, nx * ny, 1.0);
            grid.stream();
        }
        let vel = grid.velocity_field();
        let ux_avg: f64 = vel.iter().map(|v| v[0]).sum::<f64>() / vel.len() as f64;
        assert!(
            ux_avg >= 0.0,
            "x-velocity should increase with positive x body force"
        );
    }
}
/// Compute macroscopic density from a D2Q9 distribution.
pub fn compute_rho(f: &[f64; 9]) -> f64 {
    f.iter().sum()
}
/// Compute macroscopic velocity `[ux, uy]` from a D2Q9 distribution.
///
/// Requires `rho > 0`.
pub fn compute_velocity(f: &[f64; 9], rho: f64) -> [f64; 2] {
    let mut mx = 0.0_f64;
    let mut my = 0.0_f64;
    for i in 0..9 {
        mx += f[i] * D2Q9_VELOCITIES[i][0] as f64;
        my += f[i] * D2Q9_VELOCITIES[i][1] as f64;
    }
    if rho.abs() > 1e-15 {
        [mx / rho, my / rho]
    } else {
        [0.0, 0.0]
    }
}
/// Apply BGK collision in-place to a single D2Q9 node.
///
/// Relaxes `f` toward equilibrium at `(rho, ux, uy)` with frequency `omega`.
pub fn bgk_collision(f: &mut [f64; 9], rho: f64, ux: f64, uy: f64, omega: f64) {
    let feq = D2Q9::equilibrium(rho, ux, uy);
    for i in 0..9 {
        f[i] -= omega * (f[i] - feq[i]);
    }
}
/// D2Q9 full periodic pull-scheme streaming step.
///
/// `f` is laid out as `f[y * nx * 9 + x * 9 + alpha]` (flat, row-major).
/// After streaming, distribution `alpha` at `(x, y)` is pulled from
/// `(x - cx_alpha, y - cy_alpha)` (periodic).
pub fn stream_d2q9(f: &mut [f64], nx: usize, ny: usize) {
    let f_old = f.to_owned();
    for y in 0..ny {
        for x in 0..nx {
            let dst = (y * nx + x) * 9;
            for alpha in 0..9 {
                let cx = D2Q9_VELOCITIES[alpha][0];
                let cy = D2Q9_VELOCITIES[alpha][1];
                let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                let src = (sy * nx + sx) * 9 + alpha;
                f[dst + alpha] = f_old[src];
            }
        }
    }
}
/// D3Q19 full periodic pull-scheme streaming step.
///
/// `f` is laid out as `f[(z * ny * nx + y * nx + x) * 19 + alpha]`.
pub fn stream_d3q19(f: &mut [f64], nx: usize, ny: usize, nz: usize) {
    let f_old = f.to_owned();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let dst = (z * ny * nx + y * nx + x) * 19;
                for alpha in 0..19 {
                    let cx = D3Q19_VELOCITIES[alpha][0];
                    let cy = D3Q19_VELOCITIES[alpha][1];
                    let cz = D3Q19_VELOCITIES[alpha][2];
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let sz = (z as isize - cz as isize).rem_euclid(nz as isize) as usize;
                    let src = (sz * ny * nx + sy * nx + sx) * 19 + alpha;
                    f[dst + alpha] = f_old[src];
                }
            }
        }
    }
}
/// Compute the D2Q9 equilibrium distribution for density `rho` and velocity `(ux, uy)`.
///
/// Returns `[feq_0, ..., feq_8]`.
pub fn feq_d2q9(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u_sq = ux * ux + uy * uy;
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        let eu = cx * ux + cy * uy;
        feq[i] = D2Q9_WEIGHTS[i]
            * rho
            * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
    }
    feq
}
/// Compute the D3Q19 equilibrium distribution for density `rho` and velocity `(ux, uy, uz)`.
///
/// Returns `[feq_0, ..., feq_18]`.
pub fn feq_d3q19(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; 19] {
    let u_sq = ux * ux + uy * uy + uz * uz;
    let mut feq = [0.0_f64; 19];
    for i in 0..19 {
        let cx = D3Q19_VELOCITIES[i][0] as f64;
        let cy = D3Q19_VELOCITIES[i][1] as f64;
        let cz = D3Q19_VELOCITIES[i][2] as f64;
        let eu = cx * ux + cy * uy + cz * uz;
        feq[i] = D3Q19_WEIGHTS[i]
            * rho
            * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
    }
    feq
}
/// Extract macroscopic variables (rho, ux, uy) from a D2Q9 distribution slice.
///
/// `f` must have length 9.
/// Returns `(rho, ux, uy)`.
pub fn macros_d2q9(f: &[f64]) -> (f64, f64, f64) {
    debug_assert_eq!(f.len(), 9, "macros_d2q9: f must have 9 elements");
    let mut rho = 0.0_f64;
    let mut mx = 0.0_f64;
    let mut my = 0.0_f64;
    for i in 0..9 {
        rho += f[i];
        mx += f[i] * D2Q9_VELOCITIES[i][0] as f64;
        my += f[i] * D2Q9_VELOCITIES[i][1] as f64;
    }
    if rho.abs() > 1e-15 {
        (rho, mx / rho, my / rho)
    } else {
        (rho, 0.0, 0.0)
    }
}
/// Apply full-way bounce-back for a wall node in D2Q9.
///
/// Reverses all distribution functions in-place: `f[alpha] ↔ f[opposite(alpha)]`.
/// This implements a no-slip boundary at the wall node.
pub fn bounce_back_d2q9(f: &mut [f64; 9]) {
    let tmp = *f;
    for a in 0..9 {
        f[a] = tmp[D2Q9_OPPOSITES[a]];
    }
}
/// Apply full-way bounce-back for a wall node in D3Q19.
///
/// Reverses all distribution functions in-place: `f[alpha] ↔ f[opposite(alpha)]`.
pub fn bounce_back_d3q19(f: &mut [f64; 19]) {
    let tmp = *f;
    for a in 0..19 {
        f[a] = tmp[D3Q19_OPPOSITES[a]];
    }
}
#[cfg(test)]
mod lattice_additions_tests {
    use super::*;
    #[test]
    fn test_d2q9_equilibrium_sums_to_rho() {
        let rho = 1.3;
        let ux = 0.05;
        let uy = -0.03;
        let feq = D2Q9::equilibrium(rho, ux, uy);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "feq sum={sum}, expected {rho}");
    }
    #[test]
    fn test_d3q19_equilibrium_sums_to_rho() {
        let rho = 2.1;
        let u = [0.04, -0.02, 0.01];
        let feq = D3Q19::equilibrium(rho, u);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-12,
            "D3Q19 feq sum={sum}, expected {rho}"
        );
    }
    #[test]
    fn test_compute_rho_and_velocity_at_rest() {
        let f = D2Q9::equilibrium(1.0, 0.0, 0.0);
        let rho = compute_rho(&f);
        assert!((rho - 1.0).abs() < 1e-14, "rho={rho}");
        let vel = compute_velocity(&f, rho);
        assert!(vel[0].abs() < 1e-14 && vel[1].abs() < 1e-14, "vel={vel:?}");
    }
    #[test]
    fn test_bgk_collision_conserves_mass() {
        let mut f = D2Q9::equilibrium(1.5, 0.05, 0.0);
        f[0] += 0.1;
        f[1] -= 0.05;
        f[2] -= 0.05;
        let mass_before: f64 = f.iter().sum();
        let rho = compute_rho(&f);
        let vel = compute_velocity(&f, rho);
        bgk_collision(&mut f, rho, vel[0], vel[1], 1.0);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-13,
            "BGK mass not conserved"
        );
    }
    #[test]
    fn test_lattice_d2q9_grid_new_uniform_density() {
        let grid = LatticeD2Q9Grid::new(5, 5);
        let rho_field = grid.density_field();
        for (k, &r) in rho_field.iter().enumerate() {
            assert!((r - 1.0).abs() < 1e-14, "density[{k}]={r}");
        }
    }
    #[test]
    fn test_lattice_d2q9_grid_step_mass_conservation() {
        let mut grid = LatticeD2Q9Grid::new(6, 6);
        let mass_before = grid.total_mass();
        for _ in 0..20 {
            grid.step(1.0);
        }
        let mass_after = grid.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Step mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d2q9_stream_periodic() {
        let nx = 4;
        let ny = 4;
        let mut f: Vec<[f64; 9]> = vec![[0.0; 9]; nx * ny];
        f[2 * nx + 2][1] = 1.0;
        let mass_before: f64 = f.iter().flat_map(|fi| fi.iter()).sum();
        D2Q9::stream(&mut f, nx, ny);
        let mass_after: f64 = f.iter().flat_map(|fi| fi.iter()).sum();
        assert!((mass_before - mass_after).abs() < 1e-14, "Stream lost mass");
        assert!(f[2 * nx + 3][1] > 0.9, "Packet should have moved East");
    }
}
#[cfg(test)]
mod freefunction_tests {
    use super::*;
    #[test]
    fn test_feq_d2q9_sums_to_rho() {
        let rho = 1.25;
        let feq = feq_d2q9(rho, 0.07, -0.03);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-13,
            "feq_d2q9 sum={sum}, expected {rho}"
        );
    }
    #[test]
    fn test_feq_d2q9_zero_velocity() {
        let rho = 0.9;
        let feq = feq_d2q9(rho, 0.0, 0.0);
        for i in 0..9 {
            let expected = D2Q9_WEIGHTS[i] * rho;
            assert!(
                (feq[i] - expected).abs() < 1e-14,
                "feq_d2q9[{i}]={}, expected {expected}",
                feq[i]
            );
        }
    }
    #[test]
    fn test_feq_d3q19_sums_to_rho() {
        let rho = 2.0;
        let feq = feq_d3q19(rho, 0.05, -0.02, 0.01);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-12,
            "feq_d3q19 sum={sum}, expected {rho}"
        );
    }
    #[test]
    fn test_feq_d3q19_zero_velocity() {
        let rho = 1.5;
        let feq = feq_d3q19(rho, 0.0, 0.0, 0.0);
        for i in 0..19 {
            let expected = D3Q19_WEIGHTS[i] * rho;
            assert!(
                (feq[i] - expected).abs() < 1e-14,
                "feq_d3q19[{i}]={}, expected {expected}",
                feq[i]
            );
        }
    }
    #[test]
    fn test_macros_d2q9_at_rest() {
        let rho0 = 1.3;
        let feq = feq_d2q9(rho0, 0.0, 0.0);
        let (rho, ux, uy) = macros_d2q9(&feq);
        assert!((rho - rho0).abs() < 1e-13, "rho={rho}");
        assert!(ux.abs() < 1e-14, "ux={ux}");
        assert!(uy.abs() < 1e-14, "uy={uy}");
    }
    #[test]
    fn test_macros_d2q9_nonzero_velocity() {
        let rho0 = 1.0;
        let ux0 = 0.06;
        let uy0 = -0.04;
        let feq = feq_d2q9(rho0, ux0, uy0);
        let (rho, ux, uy) = macros_d2q9(&feq);
        assert!((rho - rho0).abs() < 1e-13, "rho={rho}");
        assert!((ux - ux0).abs() < 1e-13, "ux={ux}");
        assert!((uy - uy0).abs() < 1e-13, "uy={uy}");
    }
    #[test]
    fn test_stream_d2q9_conserves_mass() {
        let nx = 8;
        let ny = 6;
        let feq = feq_d2q9(1.0, 0.0, 0.0);
        let mut f: Vec<f64> = (0..nx * ny).flat_map(|_| feq.iter().copied()).collect();
        f[(3 * nx + 2) * 9] += 0.3;
        f[(3 * nx + 2) * 9 + 1] -= 0.15;
        f[(3 * nx + 2) * 9 + 2] -= 0.15;
        let mass_before: f64 = f.iter().sum();
        stream_d2q9(&mut f, nx, ny);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "stream_d2q9 mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_stream_d2q9_east_packet() {
        let nx = 6;
        let ny = 4;
        let mut f = vec![0.0_f64; nx * ny * 9];
        f[(2 * nx + 2) * 9 + 1] = 1.0;
        stream_d2q9(&mut f, nx, ny);
        assert!(
            f[(2 * nx + 3) * 9 + 1] > 0.9,
            "East packet should be at (3,2) after streaming"
        );
    }
    #[test]
    fn test_stream_d3q19_conserves_mass() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let feq = feq_d3q19(1.0, 0.0, 0.0, 0.0);
        let mut f: Vec<f64> = (0..nx * ny * nz)
            .flat_map(|_| feq.iter().copied())
            .collect();
        f[0] += 0.1;
        let mass_before: f64 = f.iter().sum();
        stream_d3q19(&mut f, nx, ny, nz);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-11,
            "stream_d3q19 mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_bounce_back_d2q9_reverses_directions() {
        let mut f = [0.0_f64; 9];
        f[1] = 1.0;
        bounce_back_d2q9(&mut f);
        assert!(
            (f[3] - 1.0).abs() < 1e-15,
            "East should bounce to West: f[3]={}",
            f[3]
        );
        assert!(f[1].abs() < 1e-15, "East should be zero after bounce");
    }
    #[test]
    fn test_bounce_back_d2q9_conserves_mass() {
        let feq = feq_d2q9(1.2, 0.05, -0.02);
        let mut f = feq;
        let mass_before: f64 = f.iter().sum();
        bounce_back_d2q9(&mut f);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-14,
            "bounce_back_d2q9 should conserve mass"
        );
    }
}
/// Compute the zeroth-order Hermite expansion coefficient `a^(0)` = rho.
///
/// The zeroth-order coefficient is simply the macroscopic density.
pub fn hermite_a0(f: &[f64; 9]) -> f64 {
    f.iter().sum()
}
/// Compute the first-order Hermite coefficient `a^(1)_alpha` = rho * u_alpha.
///
/// Returns `[rho*ux, rho*uy]` — the momentum vector.
pub fn hermite_a1(f: &[f64; 9]) -> [f64; 2] {
    let mut mx = 0.0_f64;
    let mut my = 0.0_f64;
    for i in 0..9 {
        mx += f[i] * D2Q9_VELOCITIES[i][0] as f64;
        my += f[i] * D2Q9_VELOCITIES[i][1] as f64;
    }
    [mx, my]
}
/// Compute the second-order symmetric Hermite coefficient `a^(2)_{alpha,beta}`.
///
/// Returns `[Pxx, Pxy, Pyy]` where:
/// - `Pxx = sum_i f_i * cx_i * cx_i`
/// - `Pxy = sum_i f_i * cx_i * cy_i`
/// - `Pyy = sum_i f_i * cy_i * cy_i`
pub fn hermite_a2(f: &[f64; 9]) -> [f64; 3] {
    let mut pxx = 0.0_f64;
    let mut pxy = 0.0_f64;
    let mut pyy = 0.0_f64;
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        pxx += f[i] * cx * cx;
        pxy += f[i] * cx * cy;
        pyy += f[i] * cy * cy;
    }
    [pxx, pxy, pyy]
}
/// Compute second-order Hermite coefficients for a D3Q19 distribution.
///
/// Returns `[Pxx, Pxy, Pxz, Pyy, Pyz, Pzz]`.
pub fn hermite_a2_3d(f: &[f64; 19]) -> [f64; 6] {
    let mut pxx = 0.0_f64;
    let mut pxy = 0.0_f64;
    let mut pxz = 0.0_f64;
    let mut pyy = 0.0_f64;
    let mut pyz = 0.0_f64;
    let mut pzz = 0.0_f64;
    for i in 0..19 {
        let cx = D3Q19_VELOCITIES[i][0] as f64;
        let cy = D3Q19_VELOCITIES[i][1] as f64;
        let cz = D3Q19_VELOCITIES[i][2] as f64;
        pxx += f[i] * cx * cx;
        pxy += f[i] * cx * cy;
        pxz += f[i] * cx * cz;
        pyy += f[i] * cy * cy;
        pyz += f[i] * cy * cz;
        pzz += f[i] * cz * cz;
    }
    [pxx, pxy, pxz, pyy, pyz, pzz]
}
/// Reconstruct the D2Q9 equilibrium distribution from Hermite coefficients.
///
/// Uses the expansion up to second order:
/// `f_eq_i = w_i * (a0 + cx*a1x/cs2 + cy*a1y/cs2
///           + (cx*cx - cs2)*a2xx/(2*cs4) + cx*cy*a2xy/cs4
///           + (cy*cy - cs2)*a2yy/(2*cs4))`
pub fn equilibrium_from_hermite(_a0: f64, _a1: [f64; 2], _a2: [f64; 3]) -> [f64; 9] {
    let rho = _a0;
    let ux = if rho.abs() > 1e-15 { _a1[0] / rho } else { 0.0 };
    let uy = if rho.abs() > 1e-15 { _a1[1] / rho } else { 0.0 };
    D2Q9::equilibrium(rho, ux, uy)
}
/// Scale a physical velocity to lattice units given a reference Mach number.
///
/// `u_lattice = u_physical / u_ref * Ma_target * cs`
/// where `cs = 1/sqrt(3)` in lattice units.
pub fn scale_to_lattice_velocity(u_physical: f64, u_ref: f64, ma_target: f64) -> f64 {
    let cs = CS2.sqrt();
    (u_physical / u_ref) * ma_target * cs
}
/// Convert lattice viscosity to relaxation time tau.
///
/// `tau = nu_lattice / cs^2 + 0.5`
pub fn viscosity_to_tau(nu_lattice: f64) -> f64 {
    nu_lattice / CS2 + 0.5
}
/// Convert relaxation time tau to lattice kinematic viscosity.
///
/// `nu_lattice = cs^2 * (tau - 0.5)`
pub fn tau_to_viscosity(tau: f64) -> f64 {
    CS2 * (tau - 0.5)
}
/// Convert relaxation frequency omega to Reynolds number given grid resolution and velocity.
///
/// `Re = u_lb * L / nu_lb`
pub fn omega_to_reynolds(omega: f64, u_lb: f64, l: usize) -> f64 {
    let nu = tau_to_viscosity(1.0 / omega);
    u_lb * l as f64 / nu
}
/// Compute the target omega for a desired Reynolds number.
///
/// `omega = 1 / (Re * nu / (u_lb * L) + 0.5)`
pub fn reynolds_to_omega(re: f64, u_lb: f64, l: usize) -> f64 {
    let nu = u_lb * l as f64 / re;
    1.0 / viscosity_to_tau(nu)
}
/// Compute the lattice Mach number for a given lattice velocity.
///
/// `Ma = u / cs`
pub fn mach_number(u: f64) -> f64 {
    u / CS2.sqrt()
}
/// Compute acoustic pressure from density using equation of state `p = cs^2 * rho`.
pub fn acoustic_pressure(rho: f64) -> f64 {
    CS2 * rho
}
/// Compute the non-equilibrium stress tensor component P_xy from f and f_eq.
///
/// `Pi_xy = sum_i (f_i - f_eq_i) * cx_i * cy_i`
pub fn non_equilibrium_stress_xy(f: &[f64; 9], feq: &[f64; 9]) -> f64 {
    let mut pi_xy = 0.0_f64;
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        pi_xy += (f[i] - feq[i]) * cx * cy;
    }
    pi_xy
}
/// Compute the effective viscosity from the non-equilibrium stress using Smagorinsky model.
///
/// `nu_eff = nu_0 + (sqrt(nu_0^2 + 18 * C_s^2 * |S|) - nu_0) / 2`
/// where `|S|` is the strain-rate magnitude estimated from Pi_neq.
pub fn effective_viscosity_smagorinsky(nu0: f64, cs_sgs: f64, pi_neq_mag: f64) -> f64 {
    let discriminant = nu0 * nu0 + 18.0 * cs_sgs * cs_sgs * pi_neq_mag;
    (discriminant.sqrt() - nu0) * 0.5 + nu0
}
/// Compute the strain rate tensor magnitude from non-equilibrium distributions (D2Q9).
///
/// `|S| = sqrt(2 * Pi_neq : Pi_neq) / (2 * rho * cs^2)`
pub fn strain_rate_magnitude_d2q9(f: &[f64; 9], feq: &[f64; 9], rho: f64, tau: f64) -> f64 {
    let mut pi_sq = 0.0_f64;
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        let fneq = f[i] - feq[i];
        pi_sq += fneq * fneq * (cx * cx + cy * cy);
    }
    let pi_mag = pi_sq.sqrt();
    pi_mag / (2.0 * rho * CS2 * tau)
}
/// Compute D2Q9 equilibrium using the full Hermite polynomial expansion up to order 2.
///
/// Equivalent to the standard BGK equilibrium but written explicitly in terms
/// of the Hermite basis functions for clarity and validation purposes.
pub fn d2q9_hermite_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u_sq = ux * ux + uy * uy;
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        let eu = cx * ux + cy * uy;
        feq[i] =
            D2Q9_WEIGHTS[i] * rho * (1.0 + eu / CS2 + (eu * eu - CS2 * u_sq) / (2.0 * CS2 * CS2));
    }
    feq
}
/// Compute D3Q19 equilibrium using explicit Hermite expansion up to order 2.
pub fn d3q19_hermite_equilibrium(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; 19] {
    let u_sq = ux * ux + uy * uy + uz * uz;
    let mut feq = [0.0_f64; 19];
    for i in 0..19 {
        let cx = D3Q19_VELOCITIES[i][0] as f64;
        let cy = D3Q19_VELOCITIES[i][1] as f64;
        let cz = D3Q19_VELOCITIES[i][2] as f64;
        let eu = cx * ux + cy * uy + cz * uz;
        feq[i] =
            D3Q19_WEIGHTS[i] * rho * (1.0 + eu / CS2 + (eu * eu - CS2 * u_sq) / (2.0 * CS2 * CS2));
    }
    feq
}
/// Compute the full pressure tensor (momentum flux) for D2Q9.
///
/// Returns `[[Pxx, Pxy\], [Pyx, Pyy]]` as a flat `[f64; 4]` in row-major order.
pub fn pressure_tensor_d2q9(f: &[f64; 9]) -> [f64; 4] {
    let mut pxx = 0.0_f64;
    let mut pxy = 0.0_f64;
    let mut pyy = 0.0_f64;
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        pxx += f[i] * cx * cx;
        pxy += f[i] * cx * cy;
        pyy += f[i] * cy * cy;
    }
    [pxx, pxy, pxy, pyy]
}
/// Compute the energy flux vector `q_alpha = sum_i f_i * c_i^2 * c_alpha` for D2Q9.
///
/// Returns `[qx, qy]`.
pub fn energy_flux_d2q9(f: &[f64; 9]) -> [f64; 2] {
    let mut qx = 0.0_f64;
    let mut qy = 0.0_f64;
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        let c_sq = cx * cx + cy * cy;
        qx += f[i] * c_sq * cx;
        qy += f[i] * c_sq * cy;
    }
    [qx, qy]
}
#[cfg(test)]
mod lattice_physics_tests {
    use super::*;
    #[test]
    fn test_d3q27_equilibrium_sums_to_rho() {
        let rho = 1.4;
        let u = [0.03, -0.02, 0.01];
        let feq = D3Q27::equilibrium(rho, u);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-12,
            "D3Q27 feq sum={sum}, expected {rho}"
        );
    }
    #[test]
    fn test_d3q27_equilibrium_zero_velocity() {
        let rho = 0.8;
        let feq = D3Q27::equilibrium(rho, [0.0, 0.0, 0.0]);
        for i in 0..27 {
            let expected = D3Q27_WEIGHTS[i] * rho;
            assert!(
                (feq[i] - expected).abs() < 1e-14,
                "D3Q27 feq[{i}]={}, expected {expected}",
                feq[i]
            );
        }
    }
    #[test]
    fn test_d3q27_compute_rho() {
        let rho0 = 1.2;
        let feq = D3Q27::equilibrium(rho0, [0.05, -0.03, 0.01]);
        let rho = D3Q27::compute_rho(&feq);
        assert!((rho - rho0).abs() < 1e-12, "D3Q27 compute_rho={rho}");
    }
    #[test]
    fn test_d3q27_compute_velocity_at_rest() {
        let rho = 1.0;
        let feq = D3Q27::equilibrium(rho, [0.0, 0.0, 0.0]);
        let vel = D3Q27::compute_velocity(&feq, rho);
        assert!(
            vel[0].abs() < 1e-14 && vel[1].abs() < 1e-14 && vel[2].abs() < 1e-14,
            "D3Q27 velocity at rest should be zero: {vel:?}"
        );
    }
    #[test]
    fn test_d3q27_compute_velocity_nonzero() {
        let rho = 1.0;
        let u0 = [0.05_f64, -0.03, 0.02];
        let feq = D3Q27::equilibrium(rho, u0);
        let vel = D3Q27::compute_velocity(&feq, rho);
        for d in 0..3 {
            assert!(
                (vel[d] - u0[d]).abs() < 1e-12,
                "D3Q27 velocity[{d}]={}, expected {}",
                vel[d],
                u0[d]
            );
        }
    }
    #[test]
    fn test_d3q27_stream_conserves_mass() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let feq = D3Q27::equilibrium(1.0, [0.0, 0.0, 0.0]);
        let mut f = vec![feq; nx * ny * nz];
        f[0][1] += 0.1;
        let mass_before: f64 = f.iter().flat_map(|fi| fi.iter()).sum();
        D3Q27::stream(&mut f, nx, ny, nz);
        let mass_after: f64 = f.iter().flat_map(|fi| fi.iter()).sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-11,
            "D3Q27 stream mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d3q27_bounce_back_conserves_mass() {
        let rho = 1.1;
        let feq = D3Q27::equilibrium(rho, [0.04, -0.02, 0.01]);
        let mut f = feq;
        let mass_before: f64 = f.iter().sum();
        D3Q27::bounce_back(&mut f);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-14,
            "D3Q27 bounce_back mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d3q27_bgk_conserves_mass() {
        let rho = 1.0;
        let u = [0.03, 0.01, -0.02];
        let mut f = D3Q27::equilibrium(rho, u);
        f[0] += 0.05;
        f[1] -= 0.05;
        let mass_before: f64 = f.iter().sum();
        D3Q27::bgk(&mut f, rho, u, 1.2);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-13,
            "D3Q27 BGK mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_hermite_a0_equals_rho() {
        let rho = 1.3;
        let f = D2Q9::equilibrium(rho, 0.05, -0.03);
        let a0 = hermite_a0(&f);
        assert!((a0 - rho).abs() < 1e-13, "hermite_a0={a0}, expected {rho}");
    }
    #[test]
    fn test_hermite_a1_equals_momentum_at_rest() {
        let rho = 1.0;
        let f = D2Q9::equilibrium(rho, 0.0, 0.0);
        let a1 = hermite_a1(&f);
        assert!(
            a1[0].abs() < 1e-14 && a1[1].abs() < 1e-14,
            "hermite_a1 at rest should be zero: {a1:?}"
        );
    }
    #[test]
    fn test_hermite_a1_equals_momentum_nonzero() {
        let rho = 1.2;
        let ux = 0.07;
        let uy = -0.04;
        let f = D2Q9::equilibrium(rho, ux, uy);
        let a1 = hermite_a1(&f);
        assert!(
            (a1[0] - rho * ux).abs() < 1e-13,
            "hermite_a1[0]={}, expected {}",
            a1[0],
            rho * ux
        );
        assert!(
            (a1[1] - rho * uy).abs() < 1e-13,
            "hermite_a1[1]={}, expected {}",
            a1[1],
            rho * uy
        );
    }
    #[test]
    fn test_hermite_a2_diagonal_at_rest() {
        let rho = 1.5;
        let f = D2Q9::equilibrium(rho, 0.0, 0.0);
        let a2 = hermite_a2(&f);
        let expected_diag = rho * CS2;
        assert!(
            (a2[0] - expected_diag).abs() < 1e-12,
            "Pxx={}, expected {}",
            a2[0],
            expected_diag
        );
        assert!(a2[1].abs() < 1e-14, "Pxy at rest should be zero: {}", a2[1]);
        assert!(
            (a2[2] - expected_diag).abs() < 1e-12,
            "Pyy={}, expected {}",
            a2[2],
            expected_diag
        );
    }
    #[test]
    fn test_hermite_a2_3d_diagonal_at_rest() {
        let rho = 1.0;
        let f = D3Q19::equilibrium(rho, [0.0, 0.0, 0.0]);
        let a2 = hermite_a2_3d(&f);
        let expected = rho * CS2;
        assert!((a2[0] - expected).abs() < 1e-12, "Pxx={}", a2[0]);
        assert!(a2[1].abs() < 1e-14, "Pxy={}", a2[1]);
        assert!(a2[2].abs() < 1e-14, "Pxz={}", a2[2]);
        assert!((a2[3] - expected).abs() < 1e-12, "Pyy={}", a2[3]);
        assert!(a2[4].abs() < 1e-14, "Pyz={}", a2[4]);
        assert!((a2[5] - expected).abs() < 1e-12, "Pzz={}", a2[5]);
    }
    #[test]
    fn test_equilibrium_from_hermite_matches_d2q9() {
        let rho = 1.1;
        let ux = 0.06;
        let uy = -0.04;
        let f = D2Q9::equilibrium(rho, ux, uy);
        let a0 = hermite_a0(&f);
        let a1 = hermite_a1(&f);
        let a2 = hermite_a2(&f);
        let feq_reconstructed = equilibrium_from_hermite(a0, a1, a2);
        for i in 0..9 {
            assert!(
                (feq_reconstructed[i] - f[i]).abs() < 1e-12,
                "feq_reconstructed[{i}]={}, original={}",
                feq_reconstructed[i],
                f[i]
            );
        }
    }
    #[test]
    fn test_viscosity_tau_roundtrip() {
        let nu = 0.01_f64;
        let tau = viscosity_to_tau(nu);
        let nu_back = tau_to_viscosity(tau);
        assert!(
            (nu - nu_back).abs() < 1e-14,
            "viscosity roundtrip: nu={nu}, tau={tau}, nu_back={nu_back}"
        );
    }
    #[test]
    fn test_tau_to_viscosity_at_half() {
        let nu = tau_to_viscosity(0.5);
        assert!(nu.abs() < 1e-15, "nu at tau=0.5 should be 0, got {nu}");
    }
    #[test]
    fn test_mach_number_at_cs() {
        let cs = CS2.sqrt();
        let ma = mach_number(cs);
        assert!((ma - 1.0).abs() < 1e-14, "Ma at cs should be 1, got {ma}");
    }
    #[test]
    fn test_acoustic_pressure_linear() {
        let rho = 2.0;
        let p = acoustic_pressure(rho);
        assert!(
            (p - CS2 * rho).abs() < 1e-14,
            "acoustic pressure={p}, expected {}",
            CS2 * rho
        );
    }
    #[test]
    fn test_reynolds_omega_roundtrip() {
        let re = 100.0_f64;
        let u_lb = 0.05_f64;
        let l = 64_usize;
        let omega = reynolds_to_omega(re, u_lb, l);
        let re_back = omega_to_reynolds(omega, u_lb, l);
        assert!(
            (re - re_back).abs() < 1e-10,
            "Reynolds roundtrip: Re={re}, back={re_back}"
        );
    }
    #[test]
    fn test_pressure_tensor_d2q9_trace() {
        let rho = 1.0;
        let p = pressure_tensor_d2q9(&D2Q9::equilibrium(rho, 0.0, 0.0));
        let trace = p[0] + p[3];
        assert!(
            (trace - 2.0 * rho * CS2).abs() < 1e-12,
            "pressure tensor trace={trace}, expected {}",
            2.0 * rho * CS2
        );
    }
    #[test]
    fn test_pressure_tensor_d2q9_symmetry() {
        let rho = 1.2;
        let p = pressure_tensor_d2q9(&D2Q9::equilibrium(rho, 0.05, 0.03));
        assert!(
            (p[1] - p[2]).abs() < 1e-14,
            "Pressure tensor not symmetric: Pxy={}, Pyx={}",
            p[1],
            p[2]
        );
    }
    #[test]
    fn test_energy_flux_d2q9_at_rest_zero() {
        let rho = 1.0;
        let q = energy_flux_d2q9(&D2Q9::equilibrium(rho, 0.0, 0.0));
        assert!(
            q[0].abs() < 1e-14 && q[1].abs() < 1e-14,
            "energy flux at rest should be zero: {q:?}"
        );
    }
    #[test]
    fn test_non_equilibrium_stress_at_equilibrium_zero() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = -0.02;
        let feq = D2Q9::equilibrium(rho, ux, uy);
        let pi_xy = non_equilibrium_stress_xy(&feq, &feq);
        assert!(
            pi_xy.abs() < 1e-14,
            "Non-eq stress at equilibrium should be zero: {pi_xy}"
        );
    }
    #[test]
    fn test_effective_viscosity_smagorinsky_no_stress() {
        let nu0 = 0.01_f64;
        let cs_sgs = 0.1_f64;
        let nu_eff = effective_viscosity_smagorinsky(nu0, cs_sgs, 0.0);
        assert!(
            (nu_eff - nu0).abs() < 1e-13,
            "With zero stress, nu_eff should equal nu0: {nu_eff}"
        );
    }
    #[test]
    fn test_d2q9_hermite_equilibrium_matches_standard() {
        let rho = 1.3;
        let ux = 0.04;
        let uy = -0.03;
        let feq_std = D2Q9::equilibrium(rho, ux, uy);
        let feq_herm = d2q9_hermite_equilibrium(rho, ux, uy);
        for i in 0..9 {
            assert!(
                (feq_std[i] - feq_herm[i]).abs() < 1e-14,
                "feq mismatch at [{i}]: std={}, hermite={}",
                feq_std[i],
                feq_herm[i]
            );
        }
    }
    #[test]
    fn test_d3q19_hermite_equilibrium_matches_standard() {
        let rho = 1.1;
        let ux = 0.03;
        let uy = -0.02;
        let uz = 0.01;
        let feq_std = D3Q19::equilibrium(rho, [ux, uy, uz]);
        let feq_herm = d3q19_hermite_equilibrium(rho, ux, uy, uz);
        for i in 0..19 {
            assert!(
                (feq_std[i] - feq_herm[i]).abs() < 1e-14,
                "D3Q19 feq mismatch at [{i}]: std={}, hermite={}",
                feq_std[i],
                feq_herm[i]
            );
        }
    }
    #[test]
    fn test_strain_rate_magnitude_at_equilibrium_near_zero() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.0;
        let tau = 1.0;
        let feq = D2Q9::equilibrium(rho, ux, uy);
        let s_mag = strain_rate_magnitude_d2q9(&feq, &feq, rho, tau);
        assert!(
            s_mag.abs() < 1e-14,
            "Strain rate at equilibrium should be ~0: {s_mag}"
        );
    }
    #[test]
    fn test_scale_to_lattice_velocity() {
        let u_lb = scale_to_lattice_velocity(1.0, 1.0, 0.1);
        let expected = 0.1 * CS2.sqrt();
        assert!(
            (u_lb - expected).abs() < 1e-14,
            "u_lb={u_lb}, expected {expected}"
        );
    }
    #[test]
    fn test_d3q27_bounce_back_reverses_east() {
        let mut f = [0.0_f64; 27];
        let mut east_idx = 0;
        let mut west_idx = 0;
        for (i, &vel) in D3Q27_VELOCITIES.iter().enumerate() {
            if vel == [1, 0, 0] {
                east_idx = i;
            }
            if vel == [-1, 0, 0] {
                west_idx = i;
            }
        }
        f[east_idx] = 1.0;
        D3Q27::bounce_back(&mut f);
        assert!(
            (f[west_idx] - 1.0).abs() < 1e-15,
            "East should bounce to West in D3Q27: f[west]={}",
            f[west_idx]
        );
    }
}
#[cfg(test)]
mod lattice_extended_tests {
    use super::*;
    #[test]
    fn test_lattice_type_q_d2q9() {
        assert_eq!(LatticeType::D2Q9.q(), 9);
    }
    #[test]
    fn test_lattice_type_q_d3q19() {
        assert_eq!(LatticeType::D3Q19.q(), 19);
    }
    #[test]
    fn test_lattice_type_q_d3q27() {
        assert_eq!(LatticeType::D3Q27.q(), 27);
    }
    #[test]
    fn test_lattice_type_dim_d2q9() {
        assert_eq!(LatticeType::D2Q9.dim(), 2);
    }
    #[test]
    fn test_lattice_type_dim_d3q19() {
        assert_eq!(LatticeType::D3Q19.dim(), 3);
    }
    #[test]
    fn test_lattice_type_dim_d3q27() {
        assert_eq!(LatticeType::D3Q27.dim(), 3);
    }
    #[test]
    fn test_d2q9_weights_sum_to_one() {
        let sum: f64 = D2Q9_WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14, "D2Q9 weights sum = {sum}");
    }
    #[test]
    fn test_d3q19_weights_sum_to_one() {
        let sum: f64 = D3Q19_WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14, "D3Q19 weights sum = {sum}");
    }
    #[test]
    fn test_d3q27_weights_sum_to_one() {
        let sum: f64 = D3Q27_WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14, "D3Q27 weights sum = {sum}");
    }
    #[test]
    fn test_d2q9_opposites_are_involutory() {
        for (i, &j) in D2Q9_OPPOSITES.iter().enumerate() {
            assert_eq!(
                D2Q9_OPPOSITES[j], i,
                "D2Q9_OPPOSITES involutory failed at {i}"
            );
        }
    }
    #[test]
    fn test_d3q19_opposites_are_involutory() {
        for (i, &j) in D3Q19_OPPOSITES.iter().enumerate() {
            assert_eq!(
                D3Q19_OPPOSITES[j], i,
                "D3Q19_OPPOSITES involutory failed at {i}"
            );
        }
    }
    #[test]
    fn test_d3q27_opposites_are_involutory() {
        for (i, &j) in D3Q27_OPPOSITES.iter().enumerate() {
            assert_eq!(
                D3Q27_OPPOSITES[j], i,
                "D3Q27_OPPOSITES involutory failed at {i}"
            );
        }
    }
    #[test]
    fn test_d2q9_opposite_velocity_negated() {
        for i in 0..9 {
            let j = D2Q9_OPPOSITES[i];
            assert_eq!(D2Q9_VELOCITIES[i][0], -D2Q9_VELOCITIES[j][0]);
            assert_eq!(D2Q9_VELOCITIES[i][1], -D2Q9_VELOCITIES[j][1]);
        }
    }
    #[test]
    fn test_d3q19_opposite_velocity_negated() {
        for i in 0..19 {
            let j = D3Q19_OPPOSITES[i];
            assert_eq!(D3Q19_VELOCITIES[i][0], -D3Q19_VELOCITIES[j][0]);
            assert_eq!(D3Q19_VELOCITIES[i][1], -D3Q19_VELOCITIES[j][1]);
            assert_eq!(D3Q19_VELOCITIES[i][2], -D3Q19_VELOCITIES[j][2]);
        }
    }
    #[test]
    fn test_d3q27_opposite_velocity_negated() {
        for i in 0..27 {
            let j = D3Q27_OPPOSITES[i];
            assert_eq!(D3Q27_VELOCITIES[i][0], -D3Q27_VELOCITIES[j][0]);
            assert_eq!(D3Q27_VELOCITIES[i][1], -D3Q27_VELOCITIES[j][1]);
            assert_eq!(D3Q27_VELOCITIES[i][2], -D3Q27_VELOCITIES[j][2]);
        }
    }
    #[test]
    fn test_d3q27_equilibrium_sums_to_rho_at_rest() {
        let rho = 1.5;
        let feq = D3Q27::equilibrium(rho, [0.0, 0.0, 0.0]);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-14, "D3Q27 feq sum at rest = {sum}");
    }
    #[test]
    fn test_d3q27_equilibrium_sums_to_rho_with_velocity() {
        let rho = 0.9;
        let feq = D3Q27::equilibrium(rho, [0.05, -0.03, 0.02]);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-13,
            "D3Q27 feq sum with velocity = {sum}"
        );
    }
    #[test]
    fn test_d3q27_equilibrium_at_zero_velocity_matches_weights() {
        let rho = 2.0;
        let feq = D3Q27::equilibrium(rho, [0.0, 0.0, 0.0]);
        for i in 0..27 {
            let expected = D3Q27_WEIGHTS[i] * rho;
            assert!(
                (feq[i] - expected).abs() < 1e-14,
                "D3Q27 feq[{i}]={}, expected {expected}",
                feq[i]
            );
        }
    }
    #[test]
    fn test_d3q27_compute_rho_recovery() {
        let rho = 1.3;
        let feq = D3Q27::equilibrium(rho, [0.01, 0.02, 0.03]);
        let rho_out = D3Q27::compute_rho(&feq);
        assert!(
            (rho_out - rho).abs() < 1e-13,
            "D3Q27 compute_rho = {rho_out}"
        );
    }
    #[test]
    fn test_d3q27_compute_velocity_recovery() {
        let rho = 1.0;
        let u_in = [0.04, -0.02, 0.01];
        let feq = D3Q27::equilibrium(rho, u_in);
        let u_out = D3Q27::compute_velocity(&feq, rho);
        for k in 0..3 {
            assert!(
                (u_out[k] - u_in[k]).abs() < 1e-13,
                "D3Q27 velocity[{k}] = {}, expected {}",
                u_out[k],
                u_in[k]
            );
        }
    }
    #[test]
    fn test_d3q27_stream_preserves_mass() {
        let nx = 3;
        let ny = 3;
        let nz = 3;
        let rho = 1.0;
        let mut f: Vec<[f64; 27]> = vec![[0.0; 27]; nx * ny * nz];
        for node in f.iter_mut() {
            *node = D3Q27::equilibrium(rho, [0.0, 0.0, 0.0]);
        }
        f[4][1] += 0.1;
        f[4][2] -= 0.1;
        let mass_before: f64 = f.iter().flat_map(|n| n.iter()).sum();
        D3Q27::stream(&mut f, nx, ny, nz);
        let mass_after: f64 = f.iter().flat_map(|n| n.iter()).sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "D3Q27 stream mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d3q27_bgk_conserves_mass() {
        let rho = 1.0;
        let u = [0.05, 0.0, 0.0];
        let mut f = D3Q27::equilibrium(rho, u);
        f[0] += 0.1;
        f[1] -= 0.05;
        f[2] -= 0.05;
        let mass_before: f64 = f.iter().sum();
        D3Q27::bgk(&mut f, rho, u, 1.0);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-13,
            "D3Q27 BGK mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_d3q27_bounce_back_involutory() {
        let mut f = [0.0_f64; 27];
        for (i, fi) in f.iter_mut().enumerate() {
            *fi = (i + 1) as f64 * 0.01;
        }
        let f_orig = f;
        D3Q27::bounce_back(&mut f);
        D3Q27::bounce_back(&mut f);
        for i in 0..27 {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-15,
                "D3Q27 bounce-back not involutory at [{i}]"
            );
        }
    }
    #[test]
    fn test_d3q19_equilibrium_momentum_consistency() {
        let rho = 1.1;
        let u = [0.03, -0.02, 0.01];
        let feq = D3Q19::equilibrium(rho, u);
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for i in 0..19 {
            mx += feq[i] * D3Q19_VELOCITIES[i][0] as f64;
            my += feq[i] * D3Q19_VELOCITIES[i][1] as f64;
            mz += feq[i] * D3Q19_VELOCITIES[i][2] as f64;
        }
        assert!((mx - rho * u[0]).abs() < 1e-13, "D3Q19 momentum x = {mx}");
        assert!((my - rho * u[1]).abs() < 1e-13, "D3Q19 momentum y = {my}");
        assert!((mz - rho * u[2]).abs() < 1e-13, "D3Q19 momentum z = {mz}");
    }
    #[test]
    fn test_d3q27_equilibrium_momentum_consistency() {
        let rho = 1.2;
        let u = [0.04, 0.01, -0.02];
        let feq = D3Q27::equilibrium(rho, u);
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for i in 0..27 {
            mx += feq[i] * D3Q27_VELOCITIES[i][0] as f64;
            my += feq[i] * D3Q27_VELOCITIES[i][1] as f64;
            mz += feq[i] * D3Q27_VELOCITIES[i][2] as f64;
        }
        assert!((mx - rho * u[0]).abs() < 1e-13, "D3Q27 momentum x = {mx}");
        assert!((my - rho * u[1]).abs() < 1e-13, "D3Q27 momentum y = {my}");
        assert!((mz - rho * u[2]).abs() < 1e-13, "D3Q27 momentum z = {mz}");
    }
    #[test]
    fn test_d2q9_equilibrium_momentum_consistency() {
        let rho = 1.4;
        let ux = 0.06;
        let uy = -0.03;
        let feq = D2Q9::equilibrium(rho, ux, uy);
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for i in 0..9 {
            mx += feq[i] * D2Q9_VELOCITIES[i][0] as f64;
            my += feq[i] * D2Q9_VELOCITIES[i][1] as f64;
        }
        assert!((mx - rho * ux).abs() < 1e-13, "D2Q9 momentum x = {mx}");
        assert!((my - rho * uy).abs() < 1e-13, "D2Q9 momentum y = {my}");
    }
    #[test]
    fn test_d2q9_equilibrium_second_moment_xx() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.0;
        let feq = D2Q9::equilibrium(rho, ux, uy);
        let expected_pxx = rho * CS2 + rho * ux * ux;
        let mut pxx = 0.0_f64;
        for i in 0..9 {
            let cx = D2Q9_VELOCITIES[i][0] as f64;
            pxx += feq[i] * cx * cx;
        }
        assert!(
            (pxx - expected_pxx).abs() < 1e-13,
            "D2Q9 P_xx = {pxx}, expected {expected_pxx}"
        );
    }
    #[test]
    fn test_lattice_struct_d2q9_weight_sum() {
        let lattice = Lattice::new(LatticeType::D2Q9);
        let sum: f64 = (0..9).map(|i| lattice.weight(i)).sum();
        assert!((sum - 1.0).abs() < 1e-14, "Lattice D2Q9 weight sum = {sum}");
    }
    #[test]
    fn test_lattice_struct_d3q19_weight_sum() {
        let lattice = Lattice::new(LatticeType::D3Q19);
        let sum: f64 = (0..19).map(|i| lattice.weight(i)).sum();
        assert!(
            (sum - 1.0).abs() < 1e-14,
            "Lattice D3Q19 weight sum = {sum}"
        );
    }
    #[test]
    fn test_lattice_struct_d3q27_weight_sum() {
        let lattice = Lattice::new(LatticeType::D3Q27);
        let sum: f64 = (0..27).map(|i| lattice.weight(i)).sum();
        assert!(
            (sum - 1.0).abs() < 1e-14,
            "Lattice D3Q27 weight sum = {sum}"
        );
    }
    #[test]
    fn test_lattice_struct_cs2() {
        let lattice = Lattice::new(LatticeType::D2Q9);
        assert!((lattice.cs2() - 1.0 / 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_lattice_struct_d2q9_opposite_consistent() {
        let lattice = Lattice::new(LatticeType::D2Q9);
        for i in 0..9 {
            let j = lattice.opposite(i);
            assert_eq!(lattice.opposite(j), i);
        }
    }
    #[test]
    fn test_lattice_struct_d3q19_opposite_consistent() {
        let lattice = Lattice::new(LatticeType::D3Q19);
        for i in 0..19 {
            let j = lattice.opposite(i);
            assert_eq!(lattice.opposite(j), i);
        }
    }
    #[test]
    fn test_lattice_dimensions_single_cell() {
        let dims = LatticeDimensions::new(1, 1, 1);
        assert_eq!(dims.total_cells(), 1);
        assert_eq!(dims.idx_3d(0, 0, 0), 0);
        assert_eq!(dims.coords(0), (0, 0, 0));
    }
    #[test]
    fn test_lattice_dimensions_large() {
        let dims = LatticeDimensions::new(100, 100, 100);
        assert_eq!(dims.total_cells(), 1_000_000);
    }
    #[test]
    fn test_lattice_dimensions_idx_boundary() {
        let dims = LatticeDimensions::new(5, 5, 5);
        let idx = dims.idx_3d(4, 4, 4);
        let (x, y, z) = dims.coords(idx);
        assert_eq!((x, y, z), (4, 4, 4));
    }
    #[test]
    fn test_feq_d2q9_sum_to_rho() {
        let rho = 1.8;
        let ux = 0.07;
        let uy = -0.05;
        let feq = feq_d2q9(rho, ux, uy);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "feq_d2q9 sum = {sum}");
    }
    #[test]
    fn test_feq_d3q19_sum_to_rho() {
        let rho = 0.8;
        let feq = feq_d3q19(rho, 0.02, -0.01, 0.03);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "feq_d3q19 sum = {sum}");
    }
    #[test]
    fn test_macros_d2q9_equilibrium_recovers_input() {
        let rho_in = 1.2;
        let ux_in = 0.04;
        let uy_in = -0.02;
        let feq = feq_d2q9(rho_in, ux_in, uy_in);
        let mut flat = [0.0_f64; 9];
        flat.copy_from_slice(&feq);
        let (rho_out, ux_out, uy_out) = macros_d2q9(&flat);
        assert!(
            (rho_out - rho_in).abs() < 1e-13,
            "macros_d2q9 rho = {rho_out}"
        );
        assert!((ux_out - ux_in).abs() < 1e-13, "macros_d2q9 ux = {ux_out}");
        assert!((uy_out - uy_in).abs() < 1e-13, "macros_d2q9 uy = {uy_out}");
    }
    #[test]
    fn test_stream_d2q9_mass_conservation() {
        let nx = 5;
        let ny = 4;
        let n = nx * ny;
        let rho = 1.0;
        let mut f = vec![0.0_f64; n * 9];
        for idx in 0..n {
            let feq = feq_d2q9(rho, 0.0, 0.0);
            f[idx * 9..idx * 9 + 9].copy_from_slice(&feq);
        }
        f[0] += 0.1;
        f[1] -= 0.1;
        let mass_before: f64 = f.iter().sum();
        stream_d2q9(&mut f, nx, ny);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "stream_d2q9 mass: {mass_before} vs {mass_after}"
        );
    }
    #[test]
    fn test_stream_d3q19_mass_conservation() {
        let nx = 3;
        let ny = 3;
        let nz = 3;
        let n = nx * ny * nz;
        let rho = 1.0;
        let mut f = vec![0.0_f64; n * 19];
        for idx in 0..n {
            let feq = feq_d3q19(rho, 0.0, 0.0, 0.0);
            f[idx * 19..idx * 19 + 19].copy_from_slice(&feq);
        }
        f[10] += 0.2;
        f[11] -= 0.2;
        let mass_before: f64 = f.iter().sum();
        stream_d3q19(&mut f, nx, ny, nz);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "stream_d3q19 mass: {mass_before} vs {mass_after}"
        );
    }
    #[test]
    fn test_compute_rho_equilibrium() {
        let rho_in = 1.7;
        let feq = feq_d2q9(rho_in, 0.02, -0.01);
        let mut f9 = [0.0_f64; 9];
        f9.copy_from_slice(&feq);
        let rho_out = compute_rho(&f9);
        assert!((rho_out - rho_in).abs() < 1e-13, "compute_rho = {rho_out}");
    }
    #[test]
    fn test_compute_velocity_equilibrium() {
        let rho = 1.0;
        let ux_in = 0.05;
        let uy_in = -0.03;
        let feq = feq_d2q9(rho, ux_in, uy_in);
        let mut f9 = [0.0_f64; 9];
        f9.copy_from_slice(&feq);
        let vel = compute_velocity(&f9, rho);
        assert!((vel[0] - ux_in).abs() < 1e-13, "ux = {}", vel[0]);
        assert!((vel[1] - uy_in).abs() < 1e-13, "uy = {}", vel[1]);
    }
    #[test]
    fn test_bgk_collision_leaves_equilibrium_unchanged() {
        let rho = 1.0;
        let ux = 0.04;
        let uy = -0.02;
        let feq = feq_d2q9(rho, ux, uy);
        let mut f9 = [0.0_f64; 9];
        f9.copy_from_slice(&feq);
        let f_before = f9;
        bgk_collision(&mut f9, rho, ux, uy, 1.5);
        for i in 0..9 {
            assert!(
                (f9[i] - f_before[i]).abs() < 1e-13,
                "BGK changed equilibrium at [{i}]"
            );
        }
    }
    #[test]
    fn test_bounce_back_d2q9_east_to_west() {
        let mut f = [0.0_f64; 9];
        f[1] = 1.0;
        bounce_back_d2q9(&mut f);
        assert!(
            (f[3] - 1.0).abs() < 1e-15,
            "bounce_back_d2q9: East->West, f[3]={}",
            f[3]
        );
        assert!(
            f[1].abs() < 1e-15,
            "Original dir should be zero after bounce"
        );
    }
    #[test]
    fn test_hermite_a0_sums_distribution() {
        let rho = 1.3;
        let feq = feq_d2q9(rho, 0.02, -0.01);
        let mut f9 = [0.0_f64; 9];
        f9.copy_from_slice(&feq);
        let a0 = hermite_a0(&f9);
        assert!((a0 - rho).abs() < 1e-13, "hermite_a0 = {a0}");
    }
    #[test]
    fn test_hermite_a1_gives_momentum() {
        let rho = 1.0;
        let ux = 0.06;
        let uy = -0.04;
        let feq = feq_d2q9(rho, ux, uy);
        let mut f9 = [0.0_f64; 9];
        f9.copy_from_slice(&feq);
        let a1 = hermite_a1(&f9);
        assert!((a1[0] - rho * ux).abs() < 1e-13, "hermite_a1 x = {}", a1[0]);
        assert!((a1[1] - rho * uy).abs() < 1e-13, "hermite_a1 y = {}", a1[1]);
    }
    #[test]
    fn test_d3q27_bgk_approaches_equilibrium() {
        let rho = 1.0;
        let u = [0.03, 0.0, 0.0];
        let feq = D3Q27::equilibrium(rho, u);
        let mut f = feq;
        f[0] += 0.5;
        f[1] -= 0.25;
        f[2] -= 0.25;
        for _ in 0..200 {
            D3Q27::bgk(&mut f, rho, u, 1.9);
        }
        for i in 0..27 {
            assert!(
                (f[i] - feq[i]).abs() < 1e-3,
                "D3Q27 BGK did not converge to eq at [{i}]: f={}, feq={}",
                f[i],
                feq[i]
            );
        }
    }
    #[test]
    fn test_cs2_constant_value() {
        assert!((CS2 - 1.0 / 3.0).abs() < 1e-15, "CS2 = {CS2}");
    }
    #[test]
    fn test_d3q19_struct_stream_conserves_mass() {
        let nx = 4;
        let ny = 4;
        let rho = 1.0;
        let mut f: Vec<[f64; 9]> = vec![[0.0; 9]; nx * ny];
        for node in f.iter_mut() {
            let feq = feq_d2q9(rho, 0.01, 0.0);
            node.copy_from_slice(&feq);
        }
        f[5][0] += 0.3;
        f[5][1] -= 0.3;
        let mass_before: f64 = f.iter().flat_map(|n| n.iter()).sum();
        D2Q9::stream(&mut f, nx, ny);
        let mass_after: f64 = f.iter().flat_map(|n| n.iter()).sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "D2Q9::stream mass: {mass_before} vs {mass_after}"
        );
    }
}
