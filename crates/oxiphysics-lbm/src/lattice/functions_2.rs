//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions::{
    CS2, D2Q9_OPPOSITES, D2Q9_VELOCITIES, D2Q9_WEIGHTS, D3Q19_VELOCITIES, D3Q19_WEIGHTS,
    D3Q27_VELOCITIES, D3Q27_WEIGHTS,
};
#[allow(unused_imports)]
use super::streaming::{equilibrium_d2q9, equilibrium_d3q19, non_equilibrium_d2q9};
#[allow(unused_imports)]
use super::types::{D2Q9, D3Q19, D3Q27};

/// Compute the 2×2 deviatoric stress tensor components `[Pxx, Pyy, Pxy]`
/// from the non-equilibrium part of a D2Q9 distribution.
#[allow(dead_code)]
pub fn stress_tensor_d2q9(fneq: &[f64; 9]) -> [f64; 3] {
    let mut pxx = 0.0f64;
    let mut pyy = 0.0f64;
    let mut pxy = 0.0f64;
    for q in 0..9 {
        let cx = D2Q9_VELOCITIES[q][0] as f64;
        let cy = D2Q9_VELOCITIES[q][1] as f64;
        pxx += cx * cx * fneq[q];
        pyy += cy * cy * fneq[q];
        pxy += cx * cy * fneq[q];
    }
    [pxx, pyy, pxy]
}
/// Compute the shear-rate magnitude from the non-equilibrium stress tensor (2D).
///
/// `|S| = 1/(2 * rho * cs^2 * tau) * sqrt(2 * Pi_neq : Pi_neq)`
#[allow(dead_code)]
pub fn shear_rate_magnitude_d2q9(fneq: &[f64; 9], rho: f64, tau: f64) -> f64 {
    let [pxx, pyy, pxy] = stress_tensor_d2q9(fneq);
    let pi_sq = pxx * pxx + pyy * pyy + 2.0 * pxy * pxy;
    let factor = 1.0 / (2.0 * rho * CS2 * tau);
    factor * (2.0 * pi_sq).sqrt()
}
/// Compute TRT (two-relaxation-time) equilibrium split for D2Q9.
///
/// Returns symmetric (`feq_sym`) and anti-symmetric (`feq_anti`) parts.
#[allow(dead_code)]
pub fn trt_equilibrium_split_d2q9(rho: f64, ux: f64, uy: f64) -> ([f64; 9], [f64; 9]) {
    let feq = equilibrium_d2q9(rho, ux, uy);
    let mut sym = [0.0f64; 9];
    let mut anti = [0.0f64; 9];
    for q in 0..9 {
        let opp = D2Q9_OPPOSITES[q];
        sym[q] = 0.5 * (feq[q] + feq[opp]);
        anti[q] = 0.5 * (feq[q] - feq[opp]);
    }
    (sym, anti)
}
/// Apply TRT collision to a D2Q9 node.
///
/// * `omega_sym`  – relaxation rate for symmetric part (= 1/tau)
/// * `omega_anti` – relaxation rate for anti-symmetric part (controls numerical diffusion).
///   The "magic" parameter Λ = (1/omega_sym - 0.5)*(1/omega_anti - 0.5) = 3/16
///   gives exact Poiseuille flow with `omega_anti = 8*(2-omega_sym)/(8-omega_sym)`.
#[allow(dead_code)]
pub fn trt_collision_d2q9(
    f: &mut [f64; 9],
    rho: f64,
    ux: f64,
    uy: f64,
    omega_sym: f64,
    omega_anti: f64,
) {
    let (feq_sym, feq_anti) = trt_equilibrium_split_d2q9(rho, ux, uy);
    let mut f_sym = [0.0f64; 9];
    let mut f_anti = [0.0f64; 9];
    for q in 0..9 {
        let opp = D2Q9_OPPOSITES[q];
        f_sym[q] = 0.5 * (f[q] + f[opp]);
        f_anti[q] = 0.5 * (f[q] - f[opp]);
    }
    for q in 0..9 {
        f[q] = f[q] - omega_sym * (f_sym[q] - feq_sym[q]) - omega_anti * (f_anti[q] - feq_anti[q]);
    }
}
/// Compute the "magic" TRT anti-symmetric relaxation rate for Poiseuille accuracy.
///
/// Λ = (τ_sym - 0.5)(τ_anti - 0.5) = 3/16  →  τ_anti = 0.5 + 3/(16*(τ_sym-0.5))
#[allow(dead_code)]
pub fn trt_magic_omega_anti(omega_sym: f64) -> f64 {
    let tau_sym = 1.0 / omega_sym;
    let tau_anti = 0.5 + 3.0 / (16.0 * (tau_sym - 0.5));
    1.0 / tau_anti
}
/// Zou–He inlet BC: prescribe velocity `(ux, uy)` on the **west** (x=0) boundary.
///
/// Assumes the missing populations at the inlet are `q = 1, 5, 8` (eastward).
/// Computes `rho` from the known and wall populations, then fills the missing ones.
#[allow(dead_code)]
pub fn zou_he_west_velocity(f: &mut [f64; 9], ux: f64, uy: f64) {
    let rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - ux);
    let ru = rho * ux;
    f[1] = f[3] + (2.0 / 3.0) * ru;
    f[5] = f[7] - 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * ru + 0.5 * rho * uy;
    f[8] = f[6] + 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * ru - 0.5 * rho * uy;
}
/// Zou–He outlet BC: prescribe velocity `(ux, uy)` on the **east** (x=nx-1) boundary.
#[allow(dead_code)]
pub fn zou_he_east_velocity(f: &mut [f64; 9], ux: f64, uy: f64) {
    let rho = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + ux);
    let ru = rho * ux;
    f[3] = f[1] - (2.0 / 3.0) * ru;
    f[7] = f[5] + 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * ru - 0.5 * rho * uy;
    f[6] = f[8] - 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * ru + 0.5 * rho * uy;
}
/// Zou–He south-wall BC: prescribe velocity on the **south** (y=0) boundary.
#[allow(dead_code)]
pub fn zou_he_south_velocity(f: &mut [f64; 9], ux: f64, uy: f64) {
    let rho = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / (1.0 - uy);
    let rv = rho * uy;
    f[2] = f[4] + (2.0 / 3.0) * rv;
    f[5] = f[7] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux + (1.0 / 6.0) * rv;
    f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux + (1.0 / 6.0) * rv;
}
/// Convert physical kinematic viscosity to lattice relaxation time τ.
///
/// τ = 0.5 + ν_phys / (cs² * dx² / dt)
#[allow(dead_code)]
pub fn physical_nu_to_tau(nu_phys: f64, dx: f64, dt: f64) -> f64 {
    let nu_lb = nu_phys * dt / (dx * dx);
    0.5 + nu_lb / CS2
}
/// Convert τ back to lattice kinematic viscosity.
#[allow(dead_code)]
pub fn tau_to_nu_lb(tau: f64) -> f64 {
    CS2 * (tau - 0.5)
}
/// Convert lattice velocity to physical velocity.
#[allow(dead_code)]
pub fn lb_to_physical_velocity(u_lb: f64, dx: f64, dt: f64) -> f64 {
    u_lb * dx / dt
}
/// Convert lattice density to physical pressure (isothermal EOS: p = rho * cs²).
#[allow(dead_code)]
pub fn lb_to_physical_pressure(rho_lb: f64, rho0: f64, cs_phys: f64) -> f64 {
    (rho_lb - rho0) * cs_phys * cs_phys
}
/// Return the 0th moment (density) of a D2Q9 population.
#[allow(dead_code)]
pub fn moment0_d2q9(f: &[f64; 9]) -> f64 {
    f.iter().sum()
}
/// Return the 1st-order moments (momentum) of a D2Q9 population as `[jx, jy]`.
#[allow(dead_code)]
pub fn moment1_d2q9(f: &[f64; 9]) -> [f64; 2] {
    let mut jx = 0.0f64;
    let mut jy = 0.0f64;
    for q in 0..9 {
        jx += D2Q9_VELOCITIES[q][0] as f64 * f[q];
        jy += D2Q9_VELOCITIES[q][1] as f64 * f[q];
    }
    [jx, jy]
}
/// Return the 2nd-order moments (stress) as `[Pxx, Pyy, Pxy]`.
#[allow(dead_code)]
pub fn moment2_d2q9(f: &[f64; 9]) -> [f64; 3] {
    let mut pxx = 0.0f64;
    let mut pyy = 0.0f64;
    let mut pxy = 0.0f64;
    for q in 0..9 {
        let cx = D2Q9_VELOCITIES[q][0] as f64;
        let cy = D2Q9_VELOCITIES[q][1] as f64;
        pxx += cx * cx * f[q];
        pyy += cy * cy * f[q];
        pxy += cx * cy * f[q];
    }
    [pxx, pyy, pxy]
}
/// Compute the H-function (Boltzmann entropy) for a D2Q9 population.
///
/// H = Σ_i f_i * ln(f_i / w_i)
#[allow(dead_code)]
pub fn entropy_d2q9(f: &[f64; 9]) -> f64 {
    f.iter()
        .enumerate()
        .map(|(q, &fi)| {
            let wi = D2Q9_WEIGHTS[q];
            if fi > 0.0 && wi > 0.0 {
                fi * (fi / wi).ln()
            } else {
                0.0
            }
        })
        .sum()
}
/// Initialize a uniform equilibrium population array at rest.
///
/// Returns `Vec<[f64; 9]>` of length `nx * ny` with `feq(rho=1, u=0)`.
#[allow(dead_code)]
pub fn init_equilibrium_rest_2d(nx: usize, ny: usize) -> Vec<[f64; 9]> {
    let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
    vec![feq; nx * ny]
}
/// Initialize a uniform equilibrium population array for D3Q19.
#[allow(dead_code)]
pub fn init_equilibrium_rest_3d(nx: usize, ny: usize, nz: usize) -> Vec<[f64; 19]> {
    let feq = equilibrium_d3q19(1.0, 0.0, 0.0, 0.0);
    vec![feq; nx * ny * nz]
}
/// Set inlet velocity profile (parabolic Poiseuille) on the west face of a 2D grid.
///
/// `u_max` is the centreline velocity; the profile is `u_x(y) = u_max * 4y(H-y)/H²`.
#[allow(dead_code)]
pub fn set_poiseuille_inlet_west(f: &mut Vec<[f64; 9]>, nx: usize, ny: usize, u_max: f64) {
    for j in 0..ny {
        let y = j as f64 + 0.5;
        let h = ny as f64;
        let ux = u_max * 4.0 * y * (h - y) / (h * h);
        zou_he_west_velocity(&mut f[j * nx], ux, 0.0);
    }
}
/// Perform a full periodic streaming step for D3Q27 populations.
///
/// `f` is stored as `f[z*ny*nx + y*nx + x][q]`.
#[allow(dead_code)]
pub fn stream_d3q27_periodic(f: &mut Vec<[f64; 27]>, nx: usize, ny: usize, nz: usize) {
    use crate::lattice::D3Q27_VELOCITIES;
    let src = f.clone();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let dst = z * nx * ny + y * nx + x;
                for q in 0..27usize {
                    let sx = ((x as i32 - D3Q27_VELOCITIES[q][0]).rem_euclid(nx as i32)) as usize;
                    let sy = ((y as i32 - D3Q27_VELOCITIES[q][1]).rem_euclid(ny as i32)) as usize;
                    let sz = ((z as i32 - D3Q27_VELOCITIES[q][2]).rem_euclid(nz as i32)) as usize;
                    f[dst][q] = src[sz * nx * ny + sy * nx + sx][q];
                }
            }
        }
    }
}
#[cfg(test)]
mod extended_lattice_tests {
    use super::*;
    use crate::lattice::streaming::*;

    #[test]
    fn test_eq_d2q9_sums_to_rho() {
        let (rho, ux, uy) = (1.2, 0.05, -0.03);
        let feq = equilibrium_d2q9(rho, ux, uy);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "feq sum = {sum} ≠ {rho}");
    }
    #[test]
    fn test_eq_d2q9_zero_velocity_equals_weights() {
        let rho = 1.0;
        let feq = equilibrium_d2q9(rho, 0.0, 0.0);
        for (q, &fi) in feq.iter().enumerate() {
            assert!((fi - D2Q9_WEIGHTS[q]).abs() < 1e-14, "q={q}: fi={fi}");
        }
    }
    #[test]
    fn test_eq_d2q9_momentum_matches() {
        let (rho, ux, uy) = (0.98, 0.1, 0.04);
        let feq = equilibrium_d2q9(rho, ux, uy);
        let [jx, jy] = moment1_d2q9(&feq);
        assert!((jx - rho * ux).abs() < 1e-12, "jx mismatch");
        assert!((jy - rho * uy).abs() < 1e-12, "jy mismatch");
    }
    #[test]
    fn test_eq_d2q9_all_positive_at_small_velocity() {
        let feq = equilibrium_d2q9(1.0, 0.05, 0.02);
        for (q, &fi) in feq.iter().enumerate() {
            assert!(fi > 0.0, "feq[{q}] = {fi} ≤ 0");
        }
    }
    #[test]
    fn test_eq_d3q19_sums_to_rho() {
        let feq = equilibrium_d3q19(1.1, 0.04, -0.02, 0.01);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.1).abs() < 1e-13, "D3Q19 feq sum = {sum}");
    }
    #[test]
    fn test_eq_d3q19_zero_velocity_equals_weights() {
        let feq = equilibrium_d3q19(1.0, 0.0, 0.0, 0.0);
        for (q, &fi) in feq.iter().enumerate() {
            assert!((fi - D3Q19_WEIGHTS[q]).abs() < 1e-14, "D3Q19 q={q}");
        }
    }
    #[test]
    fn test_eq_d3q19_momentum_matches() {
        let (rho, ux, uy, uz) = (1.0, 0.07, 0.03, -0.02);
        let feq = equilibrium_d3q19(rho, ux, uy, uz);
        let jx: f64 = (0..19)
            .map(|q| D3Q19_VELOCITIES[q][0] as f64 * feq[q])
            .sum();
        let jy: f64 = (0..19)
            .map(|q| D3Q19_VELOCITIES[q][1] as f64 * feq[q])
            .sum();
        let jz: f64 = (0..19)
            .map(|q| D3Q19_VELOCITIES[q][2] as f64 * feq[q])
            .sum();
        assert!((jx - rho * ux).abs() < 1e-12, "D3Q19 jx mismatch");
        assert!((jy - rho * uy).abs() < 1e-12, "D3Q19 jy mismatch");
        assert!((jz - rho * uz).abs() < 1e-12, "D3Q19 jz mismatch");
    }
    #[test]
    fn test_macros_from_d2q9_at_rest() {
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        let (rho, ux, uy) = macros_from_d2q9(&feq);
        assert!((rho - 1.0).abs() < 1e-13);
        assert!(ux.abs() < 1e-13);
        assert!(uy.abs() < 1e-13);
    }
    #[test]
    fn test_macros_from_d2q9_recovers_velocity() {
        let feq = equilibrium_d2q9(1.05, 0.08, -0.04);
        let (rho, ux, uy) = macros_from_d2q9(&feq);
        assert!((rho - 1.05).abs() < 1e-12);
        assert!((ux - 0.08).abs() < 1e-12);
        assert!((uy + 0.04).abs() < 1e-12);
    }
    #[test]
    fn test_macros_from_d3q19_at_rest() {
        let feq = equilibrium_d3q19(1.0, 0.0, 0.0, 0.0);
        let (rho, ux, uy, uz) = macros_from_d3q19(&feq);
        assert!((rho - 1.0).abs() < 1e-13);
        assert!(ux.abs() < 1e-13);
        assert!(uy.abs() < 1e-13);
        assert!(uz.abs() < 1e-13);
    }
    #[test]
    fn test_stream_d2q9_periodic_conserves_mass() {
        let nx = 6usize;
        let ny = 4usize;
        let mut f = init_equilibrium_rest_2d(nx, ny);
        f[2][1] += 0.05;
        f[2][3] -= 0.05;
        let mass_before: f64 = f.iter().flat_map(|n| n.iter()).sum();
        stream_d2q9_periodic(&mut f, nx, ny);
        let mass_after: f64 = f.iter().flat_map(|n| n.iter()).sum();
        assert!((mass_before - mass_after).abs() < 1e-12);
    }
    #[test]
    fn test_stream_d2q9_periodic_moves_population_east() {
        let nx = 4usize;
        let ny = 1usize;
        let mut f = vec![[0.0f64; 9]; nx];
        f[0][1] = 1.0;
        stream_d2q9_periodic(&mut f, nx, ny);
        assert!(
            (f[1][1] - 1.0).abs() < 1e-14,
            "east streaming failed: {:?}",
            f
        );
        assert!(f[0][1].abs() < 1e-14);
    }
    #[test]
    fn test_stream_d3q19_periodic_conserves_mass() {
        let (nx, ny, nz) = (4, 4, 4);
        let mut f = init_equilibrium_rest_3d(nx, ny, nz);
        f[0][1] += 0.1;
        f[1][1] -= 0.1;
        let mass_before: f64 = f.iter().flat_map(|n| n.iter()).sum();
        stream_d3q19_periodic(&mut f, nx, ny, nz);
        let mass_after: f64 = f.iter().flat_map(|n| n.iter()).sum();
        assert!((mass_before - mass_after).abs() < 1e-11);
    }
    #[test]
    fn test_bgk_d2q9_converges_to_equilibrium() {
        let (rho, ux, uy) = (1.0, 0.05, 0.02);
        let feq = equilibrium_d2q9(rho, ux, uy);
        let mut f = feq;
        f[1] += 0.02;
        f[3] -= 0.02;
        for _ in 0..2000 {
            let (r, u, v) = macros_from_d2q9(&f);
            bgk_d2q9(&mut f, r, u, v, 1.8);
        }
        let (r_final, u_final, v_final) = macros_from_d2q9(&f);
        let feq_final = equilibrium_d2q9(r_final, u_final, v_final);
        for q in 0..9 {
            assert!((f[q] - feq_final[q]).abs() < 1e-6, "BGK q={q}: {}", f[q]);
        }
    }
    #[test]
    fn test_bgk_d2q9_conserves_mass() {
        let feq = equilibrium_d2q9(1.0, 0.04, 0.0);
        let mut f = feq;
        f[5] += 0.02;
        f[7] -= 0.02;
        let mass_before: f64 = f.iter().sum();
        let (r, u, v) = macros_from_d2q9(&f);
        bgk_d2q9(&mut f, r, u, v, 1.5);
        let mass_after: f64 = f.iter().sum();
        assert!((mass_before - mass_after).abs() < 1e-14);
    }
    #[test]
    fn test_bgk_d3q19_conserves_mass() {
        let mut f = equilibrium_d3q19(1.0, 0.03, 0.01, 0.0);
        f[1] += 0.02;
        f[2] -= 0.02;
        let mass_before: f64 = f.iter().sum();
        let (r, u, v, w) = macros_from_d3q19(&f);
        bgk_d3q19(&mut f, r, u, v, w, 1.5);
        let mass_after: f64 = f.iter().sum();
        assert!((mass_before - mass_after).abs() < 1e-14);
    }
    #[test]
    fn test_bounce_back_node_d2q9_conserves_mass() {
        let mut f = equilibrium_d2q9(1.0, 0.05, 0.0);
        let mass_before: f64 = f.iter().sum();
        bounce_back_node_d2q9(&mut f);
        let mass_after: f64 = f.iter().sum();
        assert!((mass_before - mass_after).abs() < 1e-14);
    }
    #[test]
    fn test_bounce_back_node_d2q9_reverses_momentum() {
        let mut f = equilibrium_d2q9(1.0, 0.05, 0.0);
        let [jx_before, jy_before] = moment1_d2q9(&f);
        bounce_back_node_d2q9(&mut f);
        let [jx_after, jy_after] = moment1_d2q9(&f);
        assert!(
            (jx_after + jx_before).abs() < 1e-13,
            "momentum not reversed"
        );
        assert!((jy_after + jy_before).abs() < 1e-13);
    }
    #[test]
    fn test_bounce_back_node_d3q19_conserves_mass() {
        let mut f = equilibrium_d3q19(1.0, 0.03, 0.01, 0.02);
        let mass_before: f64 = f.iter().sum();
        bounce_back_node_d3q19(&mut f);
        let mass_after: f64 = f.iter().sum();
        assert!((mass_before - mass_after).abs() < 1e-14);
    }
    #[test]
    fn test_trt_magic_omega_anti_poiseuille() {
        let omega_sym = 1.0;
        let omega_anti = trt_magic_omega_anti(omega_sym);
        let tau_sym = 1.0 / omega_sym;
        let tau_anti = 1.0 / omega_anti;
        let lambda = (tau_sym - 0.5) * (tau_anti - 0.5);
        assert!((lambda - 3.0 / 16.0).abs() < 1e-12, "Λ = {lambda}");
    }
    #[test]
    fn test_trt_collision_conserves_mass() {
        let (rho, ux, uy) = (1.0, 0.04, 0.02);
        let mut f = equilibrium_d2q9(rho, ux, uy);
        f[0] += 0.03;
        f[3] -= 0.03;
        let mass_before: f64 = f.iter().sum();
        trt_collision_d2q9(&mut f, rho, ux, uy, 1.5, trt_magic_omega_anti(1.5));
        let mass_after: f64 = f.iter().sum();
        assert!((mass_before - mass_after).abs() < 1e-13);
    }
    #[test]
    fn test_trt_split_symmetric_part() {
        let (sym, _anti) = trt_equilibrium_split_d2q9(1.0, 0.05, 0.0);
        for q in 0..9 {
            let opp = D2Q9_OPPOSITES[q];
            assert!((sym[q] - sym[opp]).abs() < 1e-14, "sym not symmetric q={q}");
        }
    }
    #[test]
    fn test_zou_he_west_preserves_prescribed_velocity() {
        let mut f = equilibrium_d2q9(1.0, 0.0, 0.0);
        zou_he_west_velocity(&mut f, 0.1, 0.0);
        let (_, ux, _) = macros_from_d2q9(&f);
        assert!((ux - 0.1).abs() < 1e-12, "ux = {ux}");
    }
    #[test]
    fn test_zou_he_east_preserves_prescribed_velocity() {
        let mut f = equilibrium_d2q9(1.0, 0.1, 0.0);
        zou_he_east_velocity(&mut f, 0.1, 0.0);
        let (_, ux, _) = macros_from_d2q9(&f);
        assert!((ux - 0.1).abs() < 1e-12, "ux = {ux}");
    }
    #[test]
    fn test_zou_he_south_preserves_prescribed_vy() {
        let mut f = equilibrium_d2q9(1.0, 0.0, 0.05);
        zou_he_south_velocity(&mut f, 0.0, 0.05);
        let (_, _, uy) = macros_from_d2q9(&f);
        assert!((uy - 0.05).abs() < 1e-12, "uy = {uy}");
    }
    #[test]
    fn test_stress_tensor_d2q9_zero_at_equilibrium() {
        let feq = equilibrium_d2q9(1.0, 0.04, 0.02);
        let fneq = non_equilibrium_d2q9(&feq, &feq);
        let [pxx, pyy, pxy] = stress_tensor_d2q9(&fneq);
        assert!(pxx.abs() < 1e-14);
        assert!(pyy.abs() < 1e-14);
        assert!(pxy.abs() < 1e-14);
    }
    #[test]
    fn test_shear_rate_magnitude_positive_for_non_newtonian() {
        let f = equilibrium_d2q9(1.0, 0.08, 0.0);
        let mut fpert = f;
        fpert[1] += 0.01;
        fpert[3] -= 0.01;
        let feq = equilibrium_d2q9(1.0, 0.08, 0.0);
        let fneq = non_equilibrium_d2q9(&fpert, &feq);
        let gamma = shear_rate_magnitude_d2q9(&fneq, 1.0, 1.0);
        assert!(gamma >= 0.0, "shear rate = {gamma}");
    }
    #[test]
    fn test_entropy_d2q9_nonnegative_at_equilibrium() {
        let feq = equilibrium_d2q9(1.0, 0.05, 0.02);
        let h_eq = entropy_d2q9(&feq);
        assert!(h_eq.is_finite(), "H_eq should be finite: {h_eq}");
    }
    #[test]
    fn test_physical_nu_to_tau_roundtrip() {
        let nu_phys = 1e-6f64;
        let dx = 1e-3f64;
        let dt = 1e-6f64;
        let tau = physical_nu_to_tau(nu_phys, dx, dt);
        let nu_back = tau_to_nu_lb(tau) * dx * dx / dt;
        assert!((nu_back - nu_phys).abs() / nu_phys < 1e-10, "nu roundtrip");
    }
    #[test]
    fn test_tau_to_nu_lb_positive() {
        assert!(tau_to_nu_lb(0.8) > 0.0);
    }
    #[test]
    fn test_lb_to_physical_velocity_scaling() {
        let u = lb_to_physical_velocity(0.1, 1e-3, 1e-6);
        assert!((u - 100.0).abs() < 1e-10, "u_phys = {u}");
    }
    #[test]
    fn test_moment0_d2q9_equals_sum() {
        let f = equilibrium_d2q9(1.3, 0.04, 0.0);
        let rho = moment0_d2q9(&f);
        assert!((rho - 1.3).abs() < 1e-13);
    }
    #[test]
    fn test_moment1_d2q9_gives_momentum() {
        let feq = equilibrium_d2q9(1.0, 0.06, -0.03);
        let [jx, jy] = moment1_d2q9(&feq);
        assert!((jx - 0.06).abs() < 1e-12, "jx={jx}");
        assert!((jy + 0.03).abs() < 1e-12, "jy={jy}");
    }
    #[test]
    fn test_moment2_d2q9_pxx_plus_pyy_equals_rho_cs2_at_rest() {
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        let [pxx, pyy, _] = moment2_d2q9(&feq);
        assert!((pxx - 1.0 / 3.0).abs() < 1e-13, "pxx={pxx}");
        assert!((pyy - 1.0 / 3.0).abs() < 1e-13, "pyy={pyy}");
    }
    #[test]
    fn test_set_poiseuille_inlet_west_centreline_max() {
        let (nx, ny) = (8usize, 10usize);
        let u_max = 0.1;
        let mut f = init_equilibrium_rest_2d(nx, ny);
        set_poiseuille_inlet_west(&mut f, nx, ny, u_max);
        let j = ny / 2;
        let (_, ux, _) = macros_from_d2q9(&f[j * nx]);
        assert!(ux > 0.8 * u_max, "centreline ux = {ux}");
    }
    #[test]
    fn test_stream_d3q27_periodic_conserves_mass() {
        use crate::lattice::D3Q27_WEIGHTS;
        let (nx, ny, nz) = (3, 3, 3);
        let feq_node: [f64; 27] = {
            let mut arr = [0.0f64; 27];
            for (q, &w) in D3Q27_WEIGHTS.iter().enumerate() {
                arr[q] = w;
            }
            arr
        };
        let mut f = vec![feq_node; nx * ny * nz];
        f[0][1] += 0.05;
        f[1][1] -= 0.05;
        let mass_before: f64 = f.iter().flat_map(|n| n.iter()).sum();
        stream_d3q27_periodic(&mut f, nx, ny, nz);
        let mass_after: f64 = f.iter().flat_map(|n| n.iter()).sum();
        assert!((mass_before - mass_after).abs() < 1e-12);
    }
}
